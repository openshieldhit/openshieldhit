#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "beam/osh_beam.h"
#include "beam/runtime/osh_beam_runtime.h"
#include "common/osh_diag.h"
#include "common/osh_particle_pool.h"
#include "common/osh_time.h"
#include "gemca/osh_gemca2.h"
#include "gemca/osh_gemca2_defines.h"
#include "gemca/runtime/osh_gemca_runtime.h"
#include "material/osh_material.h"
#include "material/runtime/osh_material_compile.h"
#include "openshieldhit/beam_defs.h"
#include "openshieldhit/simulation.h"
#include "physics/nuclear/osh_nuclear_handler.h"
#include "scoring/runtime/osh_scoring_compile.h"
#include "scoring/runtime/osh_scoring_postprocess.h"
#include "scoring/runtime/osh_scoring_shadow.h"
#include "scoring/save/osh_scoring_save.h"
#include "scoring/save/osh_scoring_sink.h"
#include "transport/osh_checkpoint_policy.h"
#include "transport/osh_fragment_pool.h"
#include "transport/osh_neutron_pool.h"
#include "transport/osh_run_control.h"
#include "transport/osh_transport.h"

/* ---- Private definitions of the opaque handles --------------------------- */

struct osh_results {
    unsigned long long requested_nstat;
    unsigned long long completed_nstat;
    char has_completed_run;
    struct osh_scoring_runtime const *scoring;
};

struct osh_simulation {
    struct osh_beam_workspace const *beam;
    struct osh_scoring_workspace const *scoring;
    struct osh_diag_sink const *diag;

    struct osh_beam_runtime *beam_rt;
    struct osh_gemca_runtime geom_rt;
    struct osh_material_runtime transport_tables;
    struct osh_scoring_runtime scoring_runtime;
    struct osh_transport_context transport_ctx;
    struct osh_nuclear_handler nuclear_handler;
    struct osh_fragment_pool fragment_pool;
    struct osh_particle_pool ion_pool;
    struct osh_neutron_pool neutron_pool;
    struct osh_transport_profile profile;
    struct osh_run_control run_control;
    struct osh_checkpoint_policy checkpoint_policy;

    /* Partial-result dump plumbing (issue #193).  The shadow is the reusable
     * out-of-place postprocess scratch (bound to scoring_runtime); the file sink
     * writes selected outputs to the workspace's resolved paths.  Bound once on
     * the first run that has dumping armed; freed with the simulation. */
    struct osh_scoring_shadow dump_shadow;
    struct osh_scoring_file_sink dump_file_sink;
    struct osh_scoring_sink dump_sink;
    int dump_bound; /* non-zero once dump_shadow + dump_sink are initialised */

    unsigned long long requested_nstat;
    unsigned long long completed_nstat;
    struct osh_results results;
};

/* True when any early-stop or dump trigger is armed, i.e. the transport should
 * be handed the run-control block at all.  Keeps the "is a policy active?"
 * decision in one place so the two independent setters cannot disagree. */
static int simulation_run_control_active(struct osh_simulation const *sim) {
    struct osh_run_control const *ctl = &sim->run_control;

    /* Check each trigger on its own line so the reason the block is (or is not)
     * needed is self-documenting, rather than a single five-term OR. */
    if (ctl->wall_budget_s > 0.0) {
        return 1; /* wall-time budget (--max-time / MAXTIME) */
    }
    if (ctl->should_stop != NULL) {
        return 1; /* graceful-stop callback (e.g. Ctrl-C) */
    }
    if (ctl->dump_every_s > 0.0) {
        return 1; /* wall-time dump cadence */
    }
    if (ctl->dump_every_primaries != 0u) {
        return 1; /* primary-count dump cadence */
    }
    if (ctl->should_dump != NULL) {
        return 1; /* on-demand dump callback (e.g. SIGUSR1) */
    }
    return 0;
}

static int prepared_has_voxel_body(struct osh_gemca_prepared const *gemca) {
    size_t ib;

    if (!gemca || !gemca->bodies) {
        return 0;
    }

    for (ib = 0u; ib < gemca->nbodies; ++ib) {
        if (gemca->bodies[ib] && gemca->bodies[ib]->type == OSH_GEMCA_BODY_VOX) {
            return 1;
        }
    }

    return 0;
}

/* ---- Transport pool + scratch helpers ------------------------------------ */

/* Ion wavefront width: primaries injected per pool fill — the user override or
 * compiled default, never more than nstat.  This is the cache/parallelism perf
 * knob ("primaries in flight"); it does NOT bound the pool allocation. */
static size_t simulation_ion_wavefront_width(struct osh_transport_params const *p) {
    size_t width = (p->pool_capacity != 0u) ? p->pool_capacity : (size_t) OSH_TRANSPORT_POOL_CAPACITY;
    if (p->nstat > 0u && width > p->nstat) {
        width = p->nstat;
    }
    return width;
}

/* Allocated ion-pool capacity: the wavefront width PLUS headroom for the nuclear
 * secondaries a full wavefront can inject before the next compaction.
 *
 * Clamping capacity to the wavefront width (the pre-#213 behaviour) left a full
 * primary wavefront with zero room for its recoils/abrasion nucleons/break-up
 * fragments, so they were silently dropped — and, because the drop count tracked
 * the checkpoint batch schedule, scored results depended on the batch size
 * (issue #213).  Doubling the width gives 100 % headroom: secondaries drain each
 * wavefront pass, so peak occupancy stays far below capacity in practice and the
 * default configuration drops nothing, so scored output is batch- and
 * capacity-invariant.  Any residual overflow on an aggressively small capacity is
 * counted (osh_particle_pool::n_dropped) and WARNed; it cannot perturb a surviving
 * history's stream (osh_rng_split is lineage+ordinal-keyed), so the only cost is
 * the dropped secondaries' own deposited energy. */
static size_t simulation_ion_pool_capacity(struct osh_transport_params const *p) {
    size_t const width = simulation_ion_wavefront_width(p);
    size_t const headroom = width;
    if (width > ((size_t) -1) - headroom) {
        return (size_t) -1; /* saturate rather than overflow the sum */
    }
    return width + headroom;
}

/* Neutron accumulation buffer size: must span the full ion run (all nstat batches),
 * so use nstat directly rather than the wavefront batch size. */
static size_t simulation_neutron_pool_capacity(struct osh_transport_params const *p) {
    return (p->nstat > 0u) ? p->nstat : (size_t) OSH_TRANSPORT_POOL_CAPACITY;
}

/* Allocate (or re-allocate) both pools and the shared geometry scratch.
 * On failure the caller is expected to call osh_simulation_free(). */
static enum osh_status simulation_alloc_pools(struct osh_simulation *sim) {
    size_t ion_cap;
    size_t npool_cap;
    enum osh_status rc;
    struct osh_zone_ref *zr; /* locals avoid CodeQL double-free false positive */
    double *db;              /* (sim->transport_ctx fields freed above then reused) */

    ion_cap = simulation_ion_pool_capacity(&sim->transport_ctx.params);
    npool_cap = simulation_neutron_pool_capacity(&sim->transport_ctx.params);

    /* Free existing slabs before reallocating (safe on first call: pointers are NULL). */
    osh_particle_pool_free(&sim->ion_pool);
    osh_neutron_pool_free(&sim->neutron_pool);
    free(sim->transport_ctx.zone_refs);
    free(sim->transport_ctx.dist_batch);
    sim->transport_ctx.zone_refs = NULL;
    sim->transport_ctx.dist_batch = NULL;
    sim->transport_ctx.scratch_capacity = 0u;

    rc = osh_particle_pool_init(&sim->ion_pool, ion_cap);
    if (rc != OSH_OK) {
        return rc;
    }
    rc = osh_neutron_pool_init(&sim->neutron_pool, npool_cap);
    if (rc != OSH_OK) {
        return rc;
    }

    zr = (struct osh_zone_ref *) malloc(ion_cap * sizeof(struct osh_zone_ref));
    db = (double *) malloc(ion_cap * sizeof(double));
    if (!zr || !db) {
        free(zr);
        free(db);
        return OSH_ENOMEM;
    }
    sim->transport_ctx.zone_refs = zr;
    sim->transport_ctx.dist_batch = db;
    sim->transport_ctx.scratch_capacity = ion_cap;

    /* Fill injects only a wavefront's worth of primaries; the remaining capacity
     * (ion_cap - width) is the secondary headroom that keeps issue #213's drop
     * from happening in the default configuration. */
    sim->transport_ctx.ion_wavefront_width = simulation_ion_wavefront_width(&sim->transport_ctx.params);

    sim->transport_ctx.ion_pool = &sim->ion_pool;
    sim->transport_ctx.neutron_pool = &sim->neutron_pool;
    return OSH_OK;
}

/* ---- Lifecycle ----------------------------------------------------------- */

enum osh_status osh_simulation_create(struct osh_beam_workspace *beam,
                                      struct osh_geometry_workspace *geo,
                                      struct osh_material_workspace *mat,
                                      struct osh_scoring_workspace *scoring,
                                      struct osh_diag_sink const *diag,
                                      struct osh_simulation **sim_out) {
    struct osh_simulation *sim;
    struct osh_gemca_prepared *gemca;
    unsigned int z_max;
    size_t iz;
    int has_voxel_body;
    int geometry_hu_table_type;
    enum osh_status rc;

    if (!beam || !geo || !mat || !scoring || !sim_out) {
        return OSH_EINVAL;
    }
    *sim_out = NULL;
    if (!geo->prepared) {
        OSH_DIAG_ERRORF(diag, "%s", "simulation: geometry workspace has not been prepared");
        return OSH_ESTATE;
    }
    if (!beam->prepared) {
        OSH_DIAG_ERRORF(diag, "%s", "simulation: beam workspace has not been prepared");
        return OSH_ESTATE;
    }
    if (!beam->has_primary) {
        OSH_DIAG_ERRORF(diag, "%s", "simulation: beam workspace has no primary particle defined");
        return OSH_EINVAL;
    }

    sim = (struct osh_simulation *) calloc(1, sizeof(*sim));
    if (!sim) {
        return OSH_ENOMEM;
    }
    sim->beam = beam;
    sim->scoring = scoring;
    sim->diag = diag;
    sim->requested_nstat = (unsigned long long) beam->nstat;
    sim->results.requested_nstat = sim->requested_nstat;
    sim->results.completed_nstat = 0ull;
    sim->results.has_completed_run = 0;
    sim->results.scoring = &sim->scoring_runtime;

    /* ---- 1. Zone → material index resolution ----------------------------- */

    gemca = geo->prepared;
    has_voxel_body = prepared_has_voxel_body(gemca);
    if (has_voxel_body && mat->hu_table_type == OSH_HU_TABLE_NONE) {
        OSH_DIAG_ERRORF(
            diag, "%s", "simulation: voxel geometry requires an explicit HUTABLE card in the material file");
        rc = OSH_EPARSE;
        goto fail;
    }
    geometry_hu_table_type = has_voxel_body ? mat->hu_table_type : OSH_HU_TABLE_NONE;

    for (iz = 0u; iz < gemca->nzones; ++iz) {
        struct zone *z = gemca->zones[iz];
        struct osh_material const *m;

        if (!z->material_name) {
            OSH_DIAG_ERRORF(diag, "simulation: zone '%s' has no material assigned", z->name ? z->name : "(unnamed)");
            rc = OSH_EPARSE;
            goto fail;
        }
        m = osh_material_by_name(mat, z->material_name);
        if (!m) {
            OSH_DIAG_ERRORF(diag,
                            "simulation: zone '%s': unknown material '%s'",
                            z->name ? z->name : "(unnamed)",
                            z->material_name);
            rc = OSH_EPARSE;
            goto fail;
        }
        z->material_idx = m->index;
    }

    /* ---- 2. Geometry runtime -------------------------------------------- */

    rc = osh_gemca_compile(
        gemca, geometry_hu_table_type, has_voxel_body ? mat->hu_first_material_idx : 0u, diag, &sim->geom_rt);
    if (rc != OSH_OK) {
        OSH_DIAG_ERRORF(diag, "%s", "simulation: failed to compile geometry runtime");
        goto fail;
    }
    OSH_DIAG_INFOF(diag,
                   "simulation: geometry runtime zone batch dispatcher: %s",
                   osh_gemca_runtime_zone_batch_dispatch_name(&sim->geom_rt));

    /* ---- 3. Transport tables --------------------------------------------- */

    z_max = (beam->primary.z > 0u) ? (unsigned int) beam->primary.z : 1u;
    if (beam->nuclear_inelastic && z_max < OSH_FERMI_BREAKUP_ZMAX) {
        /* Fermi break-up emits ions up to Z = OSH_FERMI_BREAKUP_ZMAX; without
         * their stopping-power columns the transport loop would silently kill
         * every charged break-up product (find_projectile_index). */
        z_max = OSH_FERMI_BREAKUP_ZMAX;
    }
    rc = osh_material_compile(mat, z_max, diag, &sim->transport_tables);
    if (rc != OSH_OK) {
        OSH_DIAG_ERRORF(diag, "%s", "simulation: failed to compile transport tables");
        goto fail;
    }

    /* ---- 4. Scoring runtime --------------------------------------------- */

    rc = osh_scoring_compile(scoring, diag, &sim->scoring_runtime);
    if (rc != OSH_OK) {
        OSH_DIAG_ERRORF(diag, "%s", "simulation: failed to compile scoring runtime");
        goto fail;
    }

    /* ---- 4b. Resolve scoring settings material names and wire tables ------ */

    for (iz = 0u; iz < sim->scoring_runtime.nsettings; ++iz) {
        struct osh_scoring_settings_runtime *sset = &sim->scoring_runtime.settings[iz];
        char const *mat_name = scoring->settings[iz].material_name;
        if (mat_name) {
            struct osh_material const *m = osh_material_by_name(mat, mat_name);
            if (!m) {
                OSH_DIAG_ERRORF(diag,
                                "scoring settings '%s': unknown material '%s'",
                                sset->name ? sset->name : "(unnamed)",
                                mat_name);
                rc = OSH_EPARSE;
                goto fail;
            }
            sset->medium = (int) m->index;
            sset->has_medium = 1;
        }
    }
    /* Validate all medium indices that were set numerically (Material name path
     * already resolves via osh_material_by_name, so it is always in-range). */
    for (iz = 0u; iz < sim->scoring_runtime.nsettings; ++iz) {
        struct osh_scoring_settings_runtime const *sset = &sim->scoring_runtime.settings[iz];
        if (sset->has_medium) {
            if (sset->medium < 0 || (size_t) sset->medium >= sim->transport_tables.nmaterials) {
                OSH_DIAG_ERRORF(diag,
                                "scoring settings '%s': medium index %d out of range [0, %zu)",
                                sset->name ? sset->name : "(unnamed)",
                                sset->medium,
                                sim->transport_tables.nmaterials);
                rc = OSH_EPARSE;
                goto fail;
            }
        }
    }
    sim->scoring_runtime.mat_tables = &sim->transport_tables;
    osh_scoring_runtime_finalize_ssets(&sim->scoring_runtime);

    /* ---- 5. Transport parameters from beam ------------------------------- */

    sim->transport_ctx.params.nstat = beam->nstat;
    sim->transport_ctx.params.deltae = beam->deltae;
    sim->transport_ctx.params.demin = beam->demin;
    sim->transport_ctx.params.tcut = beam->tcut;
    sim->transport_ctx.params.ncut = beam->ncut;
    sim->transport_ctx.params.rndseed = beam->rndseed;
    sim->transport_ctx.params.rndoffset = beam->rndoffset;

    switch (beam->scatter) {
    case OSH_BEAM_MSCAT_GAUSS:
        sim->transport_ctx.params.mcs_mode = OSH_TRANSPORT_MCS_GAUSSIAN;
        break;
    case OSH_BEAM_MSCAT_MOLIERE:
        sim->transport_ctx.params.mcs_mode = OSH_TRANSPORT_MCS_MOLIERE;
        break;
    case OSH_BEAM_MSCAT_WENTZEL:
        sim->transport_ctx.params.mcs_mode = OSH_TRANSPORT_MCS_WENTZEL;
        break;
    default:
        sim->transport_ctx.params.mcs_mode = OSH_TRANSPORT_MCS_OFF;
        break;
    }
    switch (beam->straggl) {
    case OSH_BEAM_STRAGG_GAUSS:
        sim->transport_ctx.params.straggling_mode = OSH_TRANSPORT_STRAGGLING_GAUSSIAN;
        break;
    case OSH_BEAM_STRAGG_VAVILOV:
        sim->transport_ctx.params.straggling_mode = OSH_TRANSPORT_STRAGGLING_VAVILOV;
        break;
    case OSH_BEAM_STRAGG_URBAN:
        sim->transport_ctx.params.straggling_mode = OSH_TRANSPORT_STRAGGLING_URBAN;
        break;
    default:
        sim->transport_ctx.params.straggling_mode = OSH_TRANSPORT_STRAGGLING_OFF;
        break;
    }
    sim->transport_ctx.params.nuclear_inelastic = beam->nuclear_inelastic ? 1 : 0;
    sim->transport_ctx.params.nuclear_elastic = beam->nuclear_elastic ? 1 : 0;
    sim->transport_ctx.diag = diag;
    sim->fragment_pool.n_created = 0u;
    sim->fragment_pool.n_sent_breakup = 0u;
    sim->fragment_pool.n_breakup = 0u;
    sim->transport_ctx.fragment_pool = &sim->fragment_pool;

    rc = simulation_alloc_pools(sim);
    if (rc != OSH_OK) {
        OSH_DIAG_ERRORF(diag, "%s", "simulation: failed to allocate transport pools");
        goto fail;
    }

    /* ---- 5b. Nuclear handler --------------------------------------------- */

    sim->transport_ctx.nuclear_handler = NULL;
    if (sim->transport_ctx.params.nuclear_inelastic || sim->transport_ctx.params.nuclear_elastic) {
        rc = osh_nuclear_handler_compile(mat, &sim->nuclear_handler);
        if (rc != OSH_OK) {
            OSH_DIAG_ERRORF(diag, "%s", "simulation: failed to compile nuclear handler");
            goto fail;
        }
        sim->transport_ctx.nuclear_handler = &sim->nuclear_handler;
    }

    /* ---- 6. Beam runtime ------------------------------------------------- */

    rc = osh_beam_compile(beam, &sim->beam_rt);
    if (rc != OSH_OK) {
        OSH_DIAG_ERRORF(diag, "%s", "simulation: failed to initialise beam runtime");
        goto fail;
    }

    *sim_out = sim;
    return OSH_OK;

fail:
    osh_simulation_free(sim);
    return rc;
}

enum osh_status osh_simulation_set_profiling(struct osh_simulation *sim, int enable) {
    if (!sim) {
        return OSH_EINVAL;
    }
    memset(&sim->profile, 0, sizeof(sim->profile));
    sim->transport_ctx.profile = enable ? &sim->profile : NULL;
    return OSH_OK;
}

enum osh_status osh_simulation_set_run_control(struct osh_simulation *sim,
                                               double wall_budget_s,
                                               int (*should_stop)(void *user),
                                               void *user) {
    if (!sim) {
        return OSH_EINVAL;
    }
    /* Set only the stop/budget fields — never a full osh_run_control_init(), which
     * would also wipe any dump policy configured by osh_simulation_set_dump_control().
     * The block starts calloc-zeroed with the simulation, so unset fields are
     * already at their "off" defaults.  The two setters are order-independent. */
    /* Normalise a non-positive budget to "off" (0 = unlimited). */
    if (wall_budget_s > 0.0) {
        sim->run_control.wall_budget_s = wall_budget_s;
    } else {
        sim->run_control.wall_budget_s = 0.0;
    }
    sim->run_control.should_stop = should_stop;
    sim->run_control.should_stop_user = user;
    /* Hand the block to transport only when some policy (stop or dump) is armed;
     * otherwise leave the pointer NULL so the hot path stays exactly as before. */
    if (simulation_run_control_active(sim)) {
        sim->transport_ctx.run_control = &sim->run_control;
    } else {
        sim->transport_ctx.run_control = NULL;
    }
    return OSH_OK;
}

enum osh_status osh_simulation_set_dump_control(struct osh_simulation *sim,
                                                double dump_every_s,
                                                unsigned long long dump_every_primaries,
                                                int (*should_dump)(void *user),
                                                void *user) {
    if (!sim) {
        return OSH_EINVAL;
    }
#if SIZE_MAX < ULLONG_MAX
    /* The count cadence is stored as size_t; reject a value that would wrap on the
     * narrowing cast (only reachable where size_t is narrower than ULL). */
    if (dump_every_primaries > (unsigned long long) SIZE_MAX) {
        return OSH_EINVAL;
    }
#endif
    /* Normalise a non-positive time cadence to "off" (0). */
    if (dump_every_s > 0.0) {
        sim->run_control.dump_every_s = dump_every_s;
    } else {
        sim->run_control.dump_every_s = 0.0;
    }
    sim->run_control.dump_every_primaries = (size_t) dump_every_primaries;
    sim->run_control.should_dump = should_dump;
    sim->run_control.should_dump_user = user;

    /* The dump cadence *is* the checkpoint cadence: a scheduled dump needs
     * intermediate family-complete checkpoints to fire at, so put the run in LIVE
     * batching whenever a cadence is set (overriding any prior
     * osh_simulation_set_checkpoint_policy()).  An on-demand-only trigger adds no
     * cadence, so the run stays final-only — and since the transport dump hook
     * skips the final boundary, an on-demand request with no cadence never fires.
     * That is intentional (see the API note on osh_simulation_set_dump_control):
     * SIGUSR1 is meaningful only alongside a cadence, and the CLI passes the
     * callback through only when a cadence exists. */
    osh_checkpoint_policy_init(&sim->checkpoint_policy);
    if (run_ctl_has_scheduled_dump(&sim->run_control)) {
        sim->checkpoint_policy.mode = OSH_PARTIAL_LIVE;
        sim->checkpoint_policy.completeness = OSH_PARTIAL_EXACT;
        sim->checkpoint_policy.every_s = sim->run_control.dump_every_s;
        sim->checkpoint_policy.every_primaries = sim->run_control.dump_every_primaries;
    }
    /* Only hand a LIVE policy to transport; a final-only policy leaves the pointer
     * NULL so the one-batch fast path is byte-for-byte identical to not calling this. */
    if (sim->checkpoint_policy.mode != OSH_PARTIAL_NONE) {
        sim->transport_ctx.checkpoint_policy = &sim->checkpoint_policy;
    } else {
        sim->transport_ctx.checkpoint_policy = NULL;
    }

    if (simulation_run_control_active(sim)) {
        sim->transport_ctx.run_control = &sim->run_control;
    } else {
        sim->transport_ctx.run_control = NULL;
    }
    return OSH_OK;
}

enum osh_status osh_simulation_set_checkpoint_policy(struct osh_simulation *sim, unsigned long long every_primaries) {
    if (!sim) {
        return OSH_EINVAL;
    }
#if SIZE_MAX < ULLONG_MAX
    /* The cadence is stored as size_t.  On platforms where size_t is narrower
     * than unsigned long long (e.g. 32-bit builds), reject a value that would
     * wrap on the narrowing cast so the requested cadence is never silently
     * corrupted into a smaller one.  Compiled out where size_t == ULL (the usual
     * 64-bit case), where the comparison would be trivially false. */
    if (every_primaries > (unsigned long long) SIZE_MAX) {
        return OSH_EINVAL;
    }
#endif
    osh_checkpoint_policy_init(&sim->checkpoint_policy);
    if (every_primaries > 0ull) {
        sim->checkpoint_policy.mode = OSH_PARTIAL_LIVE;
        sim->checkpoint_policy.completeness = OSH_PARTIAL_EXACT;
        sim->checkpoint_policy.every_primaries = (size_t) every_primaries;
    }
    /* Wire the policy in only when it changes behaviour (LIVE batching); a
     * final-only policy leaves the transport pointer NULL so the fast path — one
     * batch of K = nstat — is byte-for-byte identical to not calling this. */
    sim->transport_ctx.checkpoint_policy =
        (sim->checkpoint_policy.mode != OSH_PARTIAL_NONE) ? &sim->checkpoint_policy : NULL;
    return OSH_OK;
}

enum osh_status osh_simulation_set_score_replicas(struct osh_simulation *sim, size_t replicas) {
    if (!sim) {
        return OSH_EINVAL;
    }
    /* A contiguous integer partition of [0, nstat) into N parts only tiles the
     * range with no empty part when N <= nstat; reject the misuse here, at the
     * earliest point nstat is known, with a clear diagnostic. */
    if (replicas > (size_t) sim->transport_ctx.params.nstat) {
        OSH_DIAG_ERRORF(sim->diag,
                        "simulation: --score-replicas %zu exceeds nstat %zu (each replica needs at least one history)",
                        replicas,
                        (size_t) sim->transport_ctx.params.nstat);
        return OSH_EINVAL;
    }
    sim->transport_ctx.params.score_replicas = replicas;
    return OSH_OK;
}

enum osh_status osh_simulation_set_pool_capacity(struct osh_simulation *sim, size_t capacity) {
    if (!sim) {
        return OSH_EINVAL;
    }
    sim->transport_ctx.params.pool_capacity = capacity;
    return simulation_alloc_pools(sim);
}

enum osh_status osh_simulation_get_profile(struct osh_simulation const *sim, struct osh_simulation_profile *out) {
    if (!sim || !out) {
        return OSH_EINVAL;
    }
    if (!sim->transport_ctx.profile) {
        return OSH_ESTATE;
    }
    out->transport_s = sim->profile.total_s;
    out->phase_fill_s = sim->profile.fill_s;
    out->phase_zone_ref_s = sim->profile.zone_ref_s;
    out->phase_distance_s = sim->profile.distance_s;
    out->phase_step_s = sim->profile.step_s;
    out->phase_compact_s = sim->profile.compact_s;
    out->steps = sim->profile.steps;
    out->iterations = sim->profile.iterations;
    out->nuclear_events = sim->profile.nuclear_events;
    out->secondaries = sim->profile.secondaries;
    out->neutrons_banked = (unsigned long long) sim->neutron_pool.n_created;
    out->fragments_banked = (unsigned long long) sim->fragment_pool.n_created;
    out->ion_secondaries_dropped = (unsigned long long) sim->ion_pool.n_dropped;
    return OSH_OK;
}

enum osh_status osh_simulation_run(struct osh_simulation *sim) {
    enum osh_status rc;
    int dump_armed; /* set below: a dump destination is needed this run */

    if (!sim) {
        return OSH_EINVAL;
    }

    /* Variance batching (issue #209): a VARIANCE run needs >= 2 checkpoint batches
     * to have any degrees of freedom.  If the user set no other cadence and no
     * score-replica split, derive a count cadence that yields the requested number
     * of batches (ws->variance), so a plain single-threaded run produces error bars
     * out of the box.  An explicit dump/checkpoint cadence or --score-replicas
     * already provides batches, so those are left untouched (the variance fold rides
     * whatever cadence is active). */
    if (osh_scoring_runtime_tracks_variance(&sim->scoring_runtime) && sim->transport_ctx.params.score_replicas == 0u
        && sim->checkpoint_policy.mode == OSH_PARTIAL_NONE) {
        int const nbatches = sim->scoring ? sim->scoring->variance : 0;
        unsigned long long const nstat = (unsigned long long) sim->transport_ctx.params.nstat;
        if (nbatches >= 2 && nstat > 0ull) {
            /* ceil(nstat / nbatches): the count cadence that splits [0, nstat) into
             * about `nbatches` family-complete batches. */
            unsigned long long const every =
                (nstat + (unsigned long long) nbatches - 1ull) / (unsigned long long) nbatches;
            osh_simulation_set_checkpoint_policy(sim, every);
        }
    }

    /* Bind the dump destination when any dump trigger is armed.  Done here, not in
     * the setter, because the shadow aliases the compiled scoring runtime and the
     * file sink needs the workspace's resolved output paths — both stable across
     * runs, so bind once and reuse.  Fail-soft: if setup fails, warn and run
     * without dumps rather than abort (a dump is a preview, not the product).
     *
     * dump_armed: a dump trigger is present if either a cadence is set
     * (run_ctl_has_scheduled_dump) or the on-demand callback is wired.  Named so the
     * guard reads plainly instead of a four-term condition. */
    dump_armed = run_ctl_has_scheduled_dump(&sim->run_control) || sim->run_control.should_dump != NULL;
    if (sim->transport_ctx.run_control && dump_armed) {
        if (!sim->dump_bound) {
            rc = osh_scoring_shadow_init(&sim->dump_shadow, &sim->scoring_runtime);
            if (rc == OSH_OK) {
                rc = osh_scoring_file_sink_init(&sim->dump_file_sink, sim->scoring, &sim->dump_sink);
            }
            if (rc == OSH_OK) {
                sim->dump_bound = 1;
            } else {
                OSH_DIAG_WARNF(
                    sim->diag, "%s", "simulation: could not initialise partial-result dumps; continuing without them");
                osh_scoring_shadow_free(&sim->dump_shadow);
            }
        }
        if (sim->dump_bound) {
            sim->run_control.dump_sink = &sim->dump_sink;
            sim->run_control.dump_shadow = &sim->dump_shadow;
            sim->run_control.dump_outputs = NULL; /* all outputs; G2 selector reserved for the web target */
            sim->run_control.dump_noutputs = 0u;
        }
    }

    /* Arm the wall-time budget from the run's start instant, so "elapsed" in the
     * transport loop is measured against now (no-op when run control is off). */
    if (sim->transport_ctx.run_control) {
        osh_run_control_start(&sim->run_control, osh_monotonic_seconds());
    }
    sim->transport_ctx.completed_primaries = (size_t) sim->transport_ctx.params.nstat;

    rc = osh_transport_run_minimal(
        &sim->transport_ctx, sim->beam_rt, &sim->geom_rt, &sim->transport_tables, &sim->scoring_runtime);
    if (rc != OSH_OK) {
        OSH_DIAG_ERRORF(sim->diag, "%s", "simulation: transport failed");
        return rc;
    }

    /* Real number of primaries whose histories finished.  A clean stop drains
     * the pool, so this is exact; all banked secondaries from these primaries
     * are drained by the family scheduler, so the snapshot is family-exact and
     * normalising every output by it is correct (issue #192 / #195). */
    sim->completed_nstat = (unsigned long long) sim->transport_ctx.completed_primaries;

    if (sim->neutron_pool.n_created > 0u) {
        OSH_DIAG_INFOF(sim->diag,
                       "simulation: %zu neutron(s) created, %zu dropped (pool overflow)",
                       sim->neutron_pool.n_created,
                       sim->neutron_pool.n_dropped);
    }
    if (sim->ion_pool.n_dropped > 0u) {
        /* A drop means deposited energy was lost; the default configuration
         * reserves secondary headroom so this never fires (issue #213).  Warn
         * (not info) and point at the fix so it can never pass unnoticed. */
        OSH_DIAG_WARNF(sim->diag,
                       "simulation: %zu ion secondary/secondaries dropped (pool overflow); "
                       "raise the pool capacity to avoid losing their deposited energy",
                       sim->ion_pool.n_dropped);
    }
    if (sim->fragment_pool.n_sent_breakup > 0u) {
        OSH_DIAG_INFOF(sim->diag,
                       "simulation: %zu prefragment(s) sent to Fermi break-up, %zu de-excited",
                       sim->fragment_pool.n_sent_breakup,
                       sim->fragment_pool.n_breakup);
    }
    if (sim->fragment_pool.n_created > 0u) {
        OSH_DIAG_INFOF(sim->diag,
                       "simulation: %zu heavy nuclear fragment(s) emitted (Fermi break-up residues and "
                       "prefragments outside the break-up domain), transported as recoil ions or "
                       "point-deposited when below the transport threshold",
                       sim->fragment_pool.n_created);
    }

    /* Finalise per-bin Monte-Carlo standard error (issue #209) BEFORE postprocess:
     * it reads the raw sums (and, for LET/Qeff, the two-pass weight) that
     * postprocess is about to rescale/collapse.  A no-op unless the VARIANCE card
     * allocated the M2 arrays. */
    rc = osh_scoring_finalize_errors(&sim->scoring_runtime);
    if (rc != OSH_OK) {
        OSH_DIAG_ERRORF(sim->diag, "%s", "simulation: scoring error finalize failed");
        return rc;
    }

    rc = osh_scoring_postprocess(&sim->scoring_runtime);
    if (rc != OSH_OK) {
        OSH_DIAG_ERRORF(sim->diag, "%s", "simulation: scoring postprocess failed");
        return rc;
    }

    sim->results.requested_nstat = sim->requested_nstat;
    sim->results.completed_nstat = sim->completed_nstat;
    sim->results.has_completed_run = 1;
    sim->results.scoring = &sim->scoring_runtime;

    return OSH_OK;
}

enum osh_status osh_simulation_get_results(struct osh_simulation const *sim, struct osh_results const **out) {
    if (!sim || !out) {
        return OSH_EINVAL;
    }

    *out = &sim->results;
    return OSH_OK;
}

unsigned long long osh_results_requested_nstat(struct osh_results const *results) {
    if (!results) {
        return 0ull;
    }
    return results->requested_nstat;
}

unsigned long long osh_results_completed_nstat(struct osh_results const *results) {
    if (!results) {
        return 0ull;
    }
    return results->completed_nstat;
}

int osh_results_has_completed_run(struct osh_results const *results) {
    if (!results) {
        return 0;
    }
    return results->has_completed_run ? 1 : 0;
}

enum osh_status osh_simulation_save(struct osh_simulation const *sim) {
    enum osh_status rc;

    if (!sim) {
        return OSH_EINVAL;
    }
    if (!sim->results.has_completed_run) {
        OSH_DIAG_ERRORF(sim->diag, "%s", "simulation: save requested before a completed run");
        return OSH_ESTATE;
    }

    rc = osh_scoring_save(sim->scoring, &sim->scoring_runtime, sim->completed_nstat);
    if (rc != OSH_OK) {
        OSH_DIAG_ERRORF(sim->diag, "%s", "simulation: scoring save failed");
    }
    return rc;
}

enum osh_status osh_simulation_free(struct osh_simulation *sim) {
    if (!sim) {
        return OSH_OK;
    }
    osh_scoring_shadow_free(&sim->dump_shadow); /* no-op when dumps were never armed (calloc-zeroed) */
    osh_nuclear_handler_free(&sim->nuclear_handler);
    osh_scoring_runtime_free(&sim->scoring_runtime);
    osh_gemca_runtime_free(&sim->geom_rt);
    osh_material_runtime_free(&sim->transport_tables);
    osh_particle_pool_free(&sim->ion_pool);
    osh_neutron_pool_free(&sim->neutron_pool);
    free(sim->transport_ctx.zone_refs);
    free(sim->transport_ctx.dist_batch);
    osh_beam_runtime_free(&sim->beam_rt);
    free(sim);
    return OSH_OK;
}
