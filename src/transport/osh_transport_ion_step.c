#include "transport/osh_transport_ion_step.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "transport/osh_neutron_pool.h"

#include "common/osh_coord.h"
#include "common/osh_diag.h"
#include "common/osh_particle_pool.h"
#include "common/osh_ray.h"
#include "common/osh_step.h"
#include "common/osh_step_segment.h"
#include "gemca/runtime/osh_gemca_runtime.h"
#include "material/osh_material.h"
#include "material/runtime/osh_material_runtime.h"
#include "physics/atomic/osh_physics_bethe.h"
#include "physics/atomic/osh_physics_moliere.h"
#include "physics/atomic/osh_physics_straggling.h"
#include "physics/nuclear/osh_nuclear_handler.h"
#include "physics/nuclear/osh_nuclear_tripathi.h"
#include "random/osh_rng.h"
#include "scoring/runtime/osh_scoring_step.h"
#include "transport/osh_transport.h"

/*
 * Spatial boundary epsilon [cm]: particles closer than this to a zone boundary
 * are nudged forward rather than stepping into it.
 */
#define OSH_TRANSPORT_BOUNDARY_EPS 1e-8

/*
 * Relative tolerance for step-length equality comparisons [dimensionless].
 * Used to detect that the actual step equals the requested CSDA step after a
 * hinge clip.  Not a spatial quantity — do not confuse with BOUNDARY_EPS.
 */
#define OSH_TRANSPORT_STEP_LEN_REL_TOL 1.0e-12

/*
 * Maximum projected RMS scattering angle per substep [rad].
 * Steps are split so that θ₀ ≤ this threshold before MCS is sampled,
 * keeping the Highland Gaussian approximation valid.
 */
#define OSH_TRANSPORT_THETA_MAX_RAD 0.1

/*
 * Hard minimum ion kinetic energy [MeV/nucleon].
 * Particles below this threshold are killed regardless of tcut or table emin.
 * Mirrors the OSH_BEAM_TMIN floor previously imported from beam headers.
 */
#define OSH_TRANSPORT_ION_EMIN_MEV_PER_U 0.1

/* ---- Per-step computation context --------------------------------------- */

/*
 * Holds all intermediate state for one transport substep, populated
 * progressively by the five step-phase functions below.  The pool SoA arrays
 * (pool->x, pool->e, ...) are NOT cached here; each phase accesses them
 * directly via (pool, slot).
 *
 * Phase ownership:
 *   ion_step_setup()                — particle/material fields + done flag
 *   ion_step_vacuum()               — step geometry (straight flight)
 *   ion_step_length()               — step-length limits (CSDA/theta/boundary)
 *   ion_step_hinge_and_scatter()    — hinge geometry + MCS direction
 *   ion_step_energy_and_straggling() — exit energy + straggling + boundary MCS
 *   ion_step_nuclear()              — stochastic nuclear reaction kill
 *   ion_step_commit()               — scoring + pool update + nudge
 */
struct ion_step_ctx {
    /* --- Set by ion_step_setup() ----------------------------------------- */
    struct particle const *part; /* shortcut: pool->species[slot]        */
    size_t projectile_idx;       /* index into material_rt->projectile_* */
    size_t zone_idx;             /* caller-provided zone index           */
    size_t zone_material_idx;    /* zone->material_idx                   */
    size_t voxel_idx;            /* flat CT voxel index (from zone_ref)  */
    char has_voxel;              /* non-zero when voxel_idx is valid     */
    double e0;                   /* entry total kinetic energy [MeV]     */
    double a_proj;               /* mass number (float cast, ≥ 1)        */
    double cutoff;               /* energy cutoff for this particle [MeV]*/
    double demin_total;          /* minimum energy loss per step [MeV]   */
    double boundary_ds;          /* distance to zone boundary [cm]       */
    double rho;                  /* segment density [g/cm³]              */
    double mat_z_mean;
    double mat_z_over_a;
    double mat_x0_gcm2;
    double proj_mass_mev;
    char enable_mcs;
    char enable_straggling;
    char is_vacuum;          /* 1 if vacuum zone or ρ ≤ 0            */
    char done;               /* 1 if step already handled (kill/nudge/error) */
    enum osh_status done_rc; /* return value when done == 1           */

    /* --- Set by ion_step_length() ---------------------------------------- */
    double r0;                    /* CSDA range at e0 [g/cm²]             */
    double e1_target;             /* CSDA energy target (DELTAE/DEMIN)    */
    double preclip_step_len;      /* step length before hinge clip [cm]   */
    char preclip_hits_boundary;   /* 1 if preclip step is boundary-limited */
    char preclip_is_csda_limited; /* 1 if preclip step is set by CSDA     */
    char demin_limited;

    /* --- Set by ion_step_hinge_and_scatter() / ion_step_vacuum() --------- */
    double step_len; /* actual track length [cm]              */
    double h;        /* hinge position along entry dir [cm]   */
    double tail_len; /* post-hinge leg [cm]                   */
    char hit_boundary;
    double nudge_dir[3]; /* direction for post-boundary nudge     */
    double w_scat[3];    /* exit direction after MCS              */

    /* --- Set by ion_step_energy_and_straggling() ------------------------- */
    double exit_energy; /* exit total kinetic energy [MeV]       */
    double ds_gcm2;     /* areal density of actual step [g/cm²] */

    /* --- Set by ion_step_nuclear() --------------------------------------- */
    double incident_dir[3];                 /* entry direction before MCS         */
    struct osh_nuclear_event nuclear_event; /* kind=NONE until handler fires   */
};

/* ---- Forward declarations ------------------------------------------------ */

static void ion_step_setup(struct ion_step_ctx *ctx,
                           struct osh_particle_pool *pool,
                           size_t slot,
                           struct osh_zone_ref const *zone_ref,
                           struct osh_step_segment const *step_segments,
                           size_t n_step_segments,
                           struct osh_gemca_runtime const *geom_rt,
                           struct osh_transport_context const *transport_ctx,
                           struct osh_material_runtime const *material_rt);

static void ion_step_vacuum(struct ion_step_ctx *ctx);

static void ion_step_length(struct ion_step_ctx *ctx,
                            struct osh_particle_pool *pool,
                            size_t slot,
                            struct osh_transport_context *transport_ctx,
                            struct osh_material_runtime const *material_rt);

static enum osh_status ion_step_hinge_and_scatter(struct ion_step_ctx *ctx,
                                                  struct osh_particle_pool *pool,
                                                  size_t slot,
                                                  struct osh_gemca_runtime const *geom_rt,
                                                  struct osh_transport_context const *transport_ctx,
                                                  struct osh_material_runtime const *material_rt,
                                                  struct osh_rng *rng);

static void ion_step_energy_and_straggling(struct ion_step_ctx *ctx,
                                           struct osh_particle_pool *pool,
                                           size_t slot,
                                           struct osh_material_runtime const *material_rt,
                                           struct osh_rng *rng);

static void ion_step_nuclear(struct ion_step_ctx *ctx,
                             struct osh_transport_context const *transport_ctx,
                             struct osh_material_runtime const *material_rt,
                             struct osh_rng *rng);

static enum osh_status ion_step_commit(struct ion_step_ctx const *ctx,
                                       struct osh_particle_pool *pool,
                                       size_t slot,
                                       struct osh_transport_context const *transport_ctx,
                                       struct osh_scoring_runtime *score_rt);

/* Table helpers */
static double cutoff_total_energy(struct osh_transport_params const *params,
                                  struct osh_material_runtime const *material_rt,
                                  struct particle const *part);
static double energy_from_residual_range(struct osh_material_runtime const *material_rt,
                                         size_t material_idx,
                                         size_t projectile_idx,
                                         double residual_range);
static double energy_grid_value(struct osh_material_runtime const *material_rt, size_t energy_idx);
static enum osh_status find_projectile_index(struct osh_material_runtime const *material_rt,
                                             struct particle const *part,
                                             size_t *projectile_idx_out);
static int is_blackhole_material(size_t material_idx);
static int is_vacuum_material(size_t material_idx);
static void step_from_pool(struct step *st,
                           struct osh_particle_pool const *pool,
                           size_t slot,
                           double step_len,
                           double exit_energy,
                           double rho,
                           int medium,
                           int zone);

/* ---- Public API ---------------------------------------------------------- */

/**
 * Top-level ion step function: sequences the five physics phases.
 *
 * Transport is entirely encapsulated here; the wavefront loop only calls this
 * function and does not touch the physics.
 */
enum osh_status osh_transport_ion_step(struct osh_particle_pool *pool,
                                       size_t slot,
                                       struct osh_zone_ref const *zone_ref,
                                       struct osh_step_segment const *step_segments,
                                       size_t n_step_segments,
                                       struct osh_gemca_runtime const *geom_rt,
                                       struct osh_transport_context *transport_ctx,
                                       struct osh_material_runtime const *material_rt,
                                       struct osh_scoring_runtime *score_rt,
                                       struct osh_rng *rng) {
    struct ion_step_ctx ctx;
    enum osh_status rc;

    /* Phase 1 — identify particle, load material, handle early exits */
    ion_step_setup(&ctx, pool, slot, zone_ref, step_segments, n_step_segments, geom_rt, transport_ctx, material_rt);
    if (ctx.done)
        return ctx.done_rc;

    if (ctx.is_vacuum) {
        /* Phase 2a — vacuum: straight shot, no energy loss or scatter */
        ion_step_vacuum(&ctx);
    } else {
        /* Phase 2b — determine step length (boundary / CSDA / θ limits) */
        ion_step_length(&ctx, pool, slot, transport_ctx, material_rt);
        if (ctx.done)
            return ctx.done_rc;

        /* Phase 3 — sample random hinge; MCS scatter for physics-limited steps;
         *           clip post-hinge leg against zone boundary */
        rc = ion_step_hinge_and_scatter(&ctx, pool, slot, geom_rt, transport_ctx, material_rt, rng);
        if (rc != OSH_OK)
            return rc;
        if (ctx.done)
            return ctx.done_rc;

        /* Phase 4 — exit energy from CSDA residual range; Bohr straggling;
         *           end-of-step MCS for boundary-limited steps */
        ion_step_energy_and_straggling(&ctx, pool, slot, material_rt, rng);
        if (ctx.done)
            return ctx.done_rc;

        /* Phase 5 — nuclear interaction sampling (no-op when disabled) */
        ion_step_nuclear(&ctx, transport_ctx, material_rt, rng);
    }

    /* Phase 6 — score step, update pool position/energy/direction, nudge */
    rc = ion_step_commit(&ctx, pool, slot, transport_ctx, score_rt);
    if (rc != OSH_OK) {
        return rc;
    }

    /* Inject secondaries produced by a nuclear event (pp elastic recoil, future
     * fragments).  Secondaries are appended past the current wavefront and are
     * processed on the next pass.  Silently skip if the pool is full. */
    if (ctx.nuclear_event.n_secondaries > 0u) {
        size_t si;
        struct osh_nuclear_event const *ev = &ctx.nuclear_event;
        for (si = 0u; si < ev->n_secondaries; ++si) {
            size_t s;
            /* Neutrons are not transported through the CSDA charged-particle
             * loop — route them to the neutron pool instead. */
            if (transport_ctx->neutron_pool != NULL &&
                ev->secondaries[si].species != NULL &&
                ev->secondaries[si].species->pdg == 2112) {
                transport_ctx->neutron_pool->n_created++;
                continue;
            }
            if (pool->n >= pool->capacity) {
                break;
            }
            s = pool->n;
            pool->x[s] = pool->x[slot];
            pool->y[s] = pool->y[slot];
            pool->z[s] = pool->z[slot];
            pool->ux[s] = ev->secondaries[si].dir[0];
            pool->uy[s] = ev->secondaries[si].dir[1];
            pool->uz[s] = ev->secondaries[si].dir[2];
            pool->e[s] = ev->secondaries[si].energy;
            pool->wt[s] = pool->wt[slot];
            pool->prim_idx[s] = pool->prim_idx[slot];
            pool->gen[s] = (pool->gen[slot] < 255u) ? (uint8_t) (pool->gen[slot] + 1u) : 255u;
            pool->species[s] = ev->secondaries[si].species;
            pool->n++;
        }
    }

    return OSH_OK;
}

/* ---- Phase 1: setup ------------------------------------------------------ */

/**
 * Identify the particle species, load material properties, and handle all
 * early-exit conditions (unknown species, energy cutoff, invalid zone,
 * blackhole, boundary nudge).
 *
 * Sets ctx->done = 1 if the particle is handled here; the caller should
 * return ctx->done_rc without proceeding to later phases.
 */
static void ion_step_setup(struct ion_step_ctx *ctx,
                           struct osh_particle_pool *pool,
                           size_t slot,
                           struct osh_zone_ref const *zone_ref,
                           struct osh_step_segment const *step_segments,
                           size_t n_step_segments,
                           struct osh_gemca_runtime const *geom_rt,
                           struct osh_transport_context const *transport_ctx,
                           struct osh_material_runtime const *material_rt) {
    struct osh_transport_params const *params;
    enum osh_status rc;

    (void) geom_rt;
    params = transport_ctx ? &transport_ctx->params : NULL;
    ctx->done = 0;
    ctx->done_rc = OSH_OK;
    ctx->nuclear_event.kind = OSH_NUCLEAR_EVENT_NONE;
    ctx->nuclear_event.n_secondaries = 0u;
    ctx->part = pool->species[slot];
    ctx->a_proj = (ctx->part->a > 0u) ? (double) ctx->part->a : 1.0;
    ctx->e0 = pool->e[slot];
    ctx->zone_idx = zone_ref ? zone_ref->zone_idx : OSH_GEMCA_ZONE_INDEX_INVALID;
    ctx->has_voxel = zone_ref ? zone_ref->has_hu : 0;
    ctx->voxel_idx = (zone_ref && zone_ref->has_hu) ? zone_ref->voxel_idx : 0u;
    ctx->demin_total = 0.0;

    ctx->boundary_ds = (n_step_segments > 0u) ? step_segments[0].ds : 0.0;
    if (params && params->demin > 0.0f) {
        ctx->demin_total = (double) params->demin * ctx->a_proj;
    }

    rc = find_projectile_index(material_rt, ctx->part, &ctx->projectile_idx);
    if (rc != OSH_OK) {
        pool->e[slot] = 0.0; /* unknown species — kill silently */
        ctx->done = 1;
        return;
    }

    ctx->cutoff = cutoff_total_energy(params, material_rt, ctx->part);
    if (ctx->e0 <= ctx->cutoff) {
        pool->e[slot] = 0.0;
        ctx->done = 1;
        return;
    }

    if (!zone_ref || zone_ref->zone_idx == OSH_GEMCA_ZONE_INDEX_INVALID) {
        pool->e[slot] = 0.0; /* escaped geometry */
        ctx->done = 1;
        return;
    }

    ctx->zone_material_idx = zone_ref->material_idx;
    ctx->rho = osh_material_runtime_get_rho(material_rt, zone_ref);

    if (is_blackhole_material(zone_ref->material_idx)) {
        pool->e[slot] = 0.0;
        ctx->done = 1;
        return;
    }

    if (zone_ref->material_idx >= material_rt->nmaterials) {
        ctx->done = 1;
        ctx->done_rc = OSH_ESTATE;
        return;
    }

    if (ctx->boundary_ds < 0.0) {
        ctx->done = 1;
        ctx->done_rc = OSH_ESTATE;
        return;
    }
    if (ctx->boundary_ds <= OSH_TRANSPORT_BOUNDARY_EPS) {
        /* Particle is sitting on a boundary: nudge it forward and re-step
         * in the next wavefront round without scoring a zero-length step. */
        pool->x[slot] += pool->ux[slot] * OSH_TRANSPORT_BOUNDARY_EPS;
        pool->y[slot] += pool->uy[slot] * OSH_TRANSPORT_BOUNDARY_EPS;
        pool->z[slot] += pool->uz[slot] * OSH_TRANSPORT_BOUNDARY_EPS;
        ctx->done = 1;
        return;
    }

    /* Load material scalars used by later phases. */
    ctx->mat_z_mean = (double) material_rt->z_mean[zone_ref->material_idx];
    ctx->mat_z_over_a = (double) material_rt->z_over_a[zone_ref->material_idx];
    ctx->mat_x0_gcm2 = (double) material_rt->rad_length[zone_ref->material_idx];
    ctx->proj_mass_mev = material_rt->projectile_mass_mev[ctx->projectile_idx];
    ctx->enable_mcs = (params && params->mcs_mode != OSH_TRANSPORT_MCS_OFF);
    ctx->enable_straggling = (params && params->straggling_mode != OSH_TRANSPORT_STRAGGLING_OFF);
    ctx->is_vacuum = (is_vacuum_material(zone_ref->material_idx) || ctx->rho <= 0.0);

    /* Default exit direction: incident direction (overwritten by MCS phases) */
    ctx->incident_dir[0] = pool->ux[slot];
    ctx->incident_dir[1] = pool->uy[slot];
    ctx->incident_dir[2] = pool->uz[slot];
    ctx->w_scat[0] = pool->ux[slot];
    ctx->w_scat[1] = pool->uy[slot];
    ctx->w_scat[2] = pool->uz[slot];
    ctx->nudge_dir[0] = pool->ux[slot];
    ctx->nudge_dir[1] = pool->uy[slot];
    ctx->nudge_dir[2] = pool->uz[slot];
    ctx->h = 0.0;
    ctx->tail_len = 0.0;
    ctx->step_len = 0.0;
    ctx->hit_boundary = 0;
    ctx->exit_energy = ctx->e0;
    ctx->ds_gcm2 = 0.0;
    ctx->demin_limited = 0;
}

/* ---- Phase 2a: vacuum step ----------------------------------------------- */

/**
 * Vacuum: fly straight to the next boundary at constant energy.
 * No energy loss, no scatter, no straggling.
 */
static void ion_step_vacuum(struct ion_step_ctx *ctx) {
    ctx->step_len = ctx->boundary_ds;
    ctx->h = ctx->boundary_ds;
    ctx->tail_len = 0.0;
    ctx->hit_boundary = 1;
    ctx->preclip_hits_boundary = 1;
    ctx->preclip_is_csda_limited = 0;
    ctx->preclip_step_len = ctx->boundary_ds;
    ctx->exit_energy = ctx->e0;
}

/* ---- Phase 2b: step-length determination --------------------------------- */

/**
 * Determine the substep length as the minimum of three criteria:
 *   1. boundary_ds     — next geometry boundary in the incident direction
 *   2. ds_csda         — DELTAE fraction of the CSDA range, but never less
 *                        than DEMIN energy loss unless the boundary clips it
 *   3. ds_theta        — maximum step keeping θ₀ ≤ OSH_TRANSPORT_THETA_MAX_RAD
 *
 * Sets ctx->preclip_step_len, ctx->preclip_hits_boundary,
 * ctx->preclip_is_csda_limited, ctx->r0, and ctx->e1_target.
 * Kills the particle (ctx->done = 1) if the CSDA step limit is degenerate.
 */
static void ion_step_length(struct ion_step_ctx *ctx,
                            struct osh_particle_pool *pool,
                            size_t slot,
                            struct osh_transport_context *transport_ctx,
                            struct osh_material_runtime const *material_rt) {
    struct osh_transport_params const *params;
    double target_energy_loss;
    double r1_csda;
    double ds_csda;
    double ds_theta;
    double z_eff_0;

    if (!transport_ctx) {
        pool->e[slot] = 0.0;
        ctx->done = 1;
        ctx->done_rc = OSH_EINVAL;
        return;
    }
    params = &transport_ctx->params;
    target_energy_loss = ctx->e0 * (double) params->deltae;
    ctx->demin_limited = 0;
    if (ctx->demin_total > target_energy_loss) {
        target_energy_loss = ctx->demin_total;
        ctx->demin_limited = 1;
    }

    ctx->e1_target = ctx->e0 - target_energy_loss;
    if (ctx->e1_target < ctx->cutoff) {
        ctx->e1_target = ctx->cutoff;
    }
    if (ctx->e1_target >= ctx->e0) {
        pool->e[slot] = 0.0;
        ctx->done = 1;
        return;
    }

    ctx->r0 = osh_material_runtime_range_lookup(
        material_rt, ctx->zone_material_idx, ctx->projectile_idx, ctx->e0 / ctx->a_proj);
    r1_csda = osh_material_runtime_range_lookup(
        material_rt, ctx->zone_material_idx, ctx->projectile_idx, ctx->e1_target / ctx->a_proj);

    ds_csda = (ctx->r0 - r1_csda) / ctx->rho;
    if (ds_csda <= 0.0) {
        pool->e[slot] = 0.0;
        ctx->done = 1;
        return;
    }

    if (ctx->enable_mcs) {
        z_eff_0 = osh_physics_bethe_z_eff(ctx->e0 / ctx->a_proj, (double) ctx->part->z, ctx->a_proj, ctx->mat_z_mean);
        ds_theta = osh_physics_moliere_s_theta(
            ctx->e0, ctx->proj_mass_mev, z_eff_0, ctx->rho, ctx->mat_x0_gcm2, OSH_TRANSPORT_THETA_MAX_RAD);
    } else {
        ds_theta = 0.0;
    }

    ctx->preclip_step_len = ctx->boundary_ds;
    if (ds_csda < ctx->preclip_step_len) {
        ctx->preclip_step_len = ds_csda;
    }
    if (ds_theta > 0.0 && ds_theta < ctx->preclip_step_len) {
        ctx->preclip_step_len = ds_theta;
    }

    ctx->preclip_hits_boundary = (ctx->boundary_ds - ctx->preclip_step_len <= OSH_TRANSPORT_BOUNDARY_EPS);
    ctx->preclip_is_csda_limited =
        (fabs(ctx->preclip_step_len - ds_csda) <= OSH_TRANSPORT_STEP_LEN_REL_TOL * fmax(1.0, ctx->preclip_step_len));
    if (ctx->demin_limited && ctx->boundary_ds + OSH_TRANSPORT_BOUNDARY_EPS < ds_csda && transport_ctx
        && !transport_ctx->warned_boundary_demin_override) {
        OSH_DIAG_WARNF(transport_ctx->diag,
                       "transport: boundary-limited step shorter than DEMIN; allowing sub-DEMIN step near boundary "
                       "(boundary_ds=%.17g cm, demin_loss=%.17g MeV, e0=%.17g MeV, zone=%zu)",
                       ctx->boundary_ds,
                       ctx->demin_total,
                       ctx->e0,
                       ctx->zone_idx);
        transport_ctx->warned_boundary_demin_override = 1;
    }
}

/* ---- Phase 3: random hinge + MCS ---------------------------------------- */

/**
 * For boundary-limited steps: record the straight path to the boundary;
 * MCS is deferred to ion_step_energy_and_straggling() once the exit energy is
 * known.
 *
 * For physics-limited steps: sample a random hinge point uniformly along the
 * preclip step, compute the Highland θ₀ at mid-step kinematics, rotate the
 * direction there, then clip the post-hinge leg against the zone boundary in
 * the scattered direction.  The actual track length is h + tail_len.
 *
 * The random-hinge method avoids the parallel-boundary pathology of applying
 * MCS only at the step endpoint and distributes the deflection along the step.
 *
 * @par References
 * Fippel M, Soukup M. A Monte Carlo dose calculation algorithm for proton
 * therapy. Med Phys. 2004;31(8):2263-2273. doi:10.1118/1.1769631.
 */
static enum osh_status ion_step_hinge_and_scatter(struct ion_step_ctx *ctx,
                                                  struct osh_particle_pool *pool,
                                                  size_t slot,
                                                  struct osh_gemca_runtime const *geom_rt,
                                                  struct osh_transport_context const *transport_ctx,
                                                  struct osh_material_runtime const *material_rt,
                                                  struct osh_rng *rng) {
    double residual_range;
    double proposed_exit_energy;
    double proposed_ds_gcm2;
    double e_mid;
    double z_eff;
    double theta0;
    double boundary_tail_ds;
    double v_in[3];
    struct ray hinge_ray;

    if (ctx->preclip_hits_boundary) {
        /* Boundary-limited: fly straight to the boundary. */
        ctx->h = ctx->preclip_step_len;
        ctx->tail_len = 0.0;
        ctx->step_len = ctx->preclip_step_len;
        ctx->hit_boundary = 1;
        return OSH_OK;
    }

    /* Physics-limited: sample hinge, compute scatter, clip tail. */
    ctx->h = ctx->preclip_step_len * osh_rng_double(rng);
    ctx->tail_len = ctx->preclip_step_len - ctx->h;

    proposed_ds_gcm2 = ctx->rho * ctx->preclip_step_len;
    residual_range = ctx->r0 - proposed_ds_gcm2;
    if (residual_range <= 0.0) {
        pool->e[slot] = 0.0;
        ctx->done = 1;
        return OSH_OK;
    }

    /* Proposed exit energy before hinge clipping, used for mid-step kinematics. */
    if (ctx->preclip_is_csda_limited) {
        proposed_exit_energy = ctx->e1_target;
    } else {
        proposed_exit_energy =
            energy_from_residual_range(material_rt, ctx->zone_material_idx, ctx->projectile_idx, residual_range)
            * ctx->a_proj;
        if (proposed_exit_energy > ctx->e0)
            proposed_exit_energy = ctx->e0;
        if (proposed_exit_energy < ctx->cutoff)
            proposed_exit_energy = ctx->cutoff;
    }

    e_mid = 0.5 * (ctx->e0 + proposed_exit_energy);
    z_eff = osh_physics_bethe_z_eff(e_mid / ctx->a_proj, (double) ctx->part->z, ctx->a_proj, ctx->mat_z_mean);

    if (ctx->enable_mcs && ctx->mat_x0_gcm2 > 0.0) {
        theta0 = osh_physics_moliere_theta0(e_mid, ctx->proj_mass_mev, z_eff, proposed_ds_gcm2, ctx->mat_x0_gcm2);
        if (theta0 > 0.0) {
            v_in[0] = pool->ux[slot];
            v_in[1] = pool->uy[slot];
            v_in[2] = pool->uz[slot];
            osh_physics_moliere_scatter(v_in, ctx->w_scat, theta0, rng);
        }
    }

    /* Build hinge ray and query boundary distance in the scattered direction */
    hinge_ray.p[0] = pool->x[slot] + pool->ux[slot] * ctx->h;
    hinge_ray.p[1] = pool->y[slot] + pool->uy[slot] * ctx->h;
    hinge_ray.p[2] = pool->z[slot] + pool->uz[slot] * ctx->h;
    hinge_ray.cp[0] = ctx->w_scat[0];
    hinge_ray.cp[1] = ctx->w_scat[1];
    hinge_ray.cp[2] = ctx->w_scat[2];
    hinge_ray.system = OSH_COORD_UNIVERSE;

    boundary_tail_ds = osh_gemca_runtime_get_distance(geom_rt, ctx->zone_idx, &hinge_ray);
    if (boundary_tail_ds < 0.0) {
        OSH_DIAG_ERRORF(transport_ctx->diag,
                        "transport: negative hinge boundary distance zone=%zu h=%.17g tail_req=%.17g hinge=(%.17g, "
                        "%.17g, %.17g) w=(%.17g, %.17g, %.17g)",
                        ctx->zone_idx,
                        ctx->h,
                        ctx->tail_len,
                        hinge_ray.p[0],
                        hinge_ray.p[1],
                        hinge_ray.p[2],
                        hinge_ray.cp[0],
                        hinge_ray.cp[1],
                        hinge_ray.cp[2]);
        return OSH_ESTATE;
    }
    if (boundary_tail_ds < ctx->tail_len) {
        ctx->tail_len = boundary_tail_ds;
        ctx->hit_boundary = 1;
        ctx->nudge_dir[0] = ctx->w_scat[0];
        ctx->nudge_dir[1] = ctx->w_scat[1];
        ctx->nudge_dir[2] = ctx->w_scat[2];
    }
    ctx->step_len = ctx->h + ctx->tail_len;
    return OSH_OK;
}

/* ---- Phase 4: exit energy + straggling ----------------------------------- */

/**
 * Compute the exit energy from the CSDA residual range for the actual
 * (possibly hinge-clipped) step length.  Apply Bohr Gaussian straggling.
 * For boundary-limited steps, also apply end-of-step Highland/Molière MCS
 * using the post-straggling mid-step kinematics.
 */
static void ion_step_energy_and_straggling(struct ion_step_ctx *ctx,
                                           struct osh_particle_pool *pool,
                                           size_t slot,
                                           struct osh_material_runtime const *material_rt,
                                           struct osh_rng *rng) {
    double residual_range;
    double e_mid;
    double z_eff;
    double sigma_strag;
    double theta0;
    double v_in[3];

    ctx->ds_gcm2 = ctx->rho * ctx->step_len;

    residual_range = ctx->r0 - ctx->ds_gcm2;
    if (residual_range <= 0.0) {
        pool->e[slot] = 0.0;
        ctx->done = 1;
        return;
    }

    /* Exit energy: shortcut to e1_target when the actual step equals the preclip CSDA step. */
    if (ctx->preclip_is_csda_limited
        && fabs(ctx->step_len - ctx->preclip_step_len)
               <= OSH_TRANSPORT_STEP_LEN_REL_TOL * fmax(1.0, ctx->preclip_step_len)) {
        ctx->exit_energy = ctx->e1_target;
    } else {
        ctx->exit_energy =
            energy_from_residual_range(material_rt, ctx->zone_material_idx, ctx->projectile_idx, residual_range)
            * ctx->a_proj;
        if (ctx->exit_energy > ctx->e0)
            ctx->exit_energy = ctx->e0;
        if (ctx->exit_energy < ctx->cutoff)
            ctx->exit_energy = ctx->cutoff;
    }

    /* Bohr Gaussian energy straggling */
    if (ctx->enable_straggling) {
        e_mid = 0.5 * (ctx->e0 + ctx->exit_energy);
        z_eff = osh_physics_bethe_z_eff(e_mid / ctx->a_proj, (double) ctx->part->z, ctx->a_proj, ctx->mat_z_mean);
        sigma_strag = osh_physics_straggling_sigma(z_eff, ctx->mat_z_over_a, ctx->ds_gcm2);
        if (sigma_strag > 0.0) {
            ctx->exit_energy += osh_rng_gauss(rng, 0.0, sigma_strag);
            if (ctx->exit_energy > ctx->e0)
                ctx->exit_energy = ctx->e0;
            if (ctx->exit_energy < ctx->cutoff)
                ctx->exit_energy = ctx->cutoff;
        }
    }

    /* End-of-step MCS for boundary-limited steps.
     * Physics-limited steps already received their scatter at the hinge.
     * Recompute e_mid from the post-straggling exit energy for consistency. */
    if (ctx->preclip_hits_boundary && ctx->enable_mcs && ctx->mat_x0_gcm2 > 0.0) {
        /* e_mid and z_eff are recomputed here (not reused from phase 3) so that
         * MCS uses the post-straggling exit energy for consistency. */
        e_mid = 0.5 * (ctx->e0 + ctx->exit_energy);
        z_eff = osh_physics_bethe_z_eff(e_mid / ctx->a_proj, (double) ctx->part->z, ctx->a_proj, ctx->mat_z_mean);
        theta0 = osh_physics_moliere_theta0(e_mid, ctx->proj_mass_mev, z_eff, ctx->ds_gcm2, ctx->mat_x0_gcm2);
        if (theta0 > 0.0) {
            v_in[0] = pool->ux[slot];
            v_in[1] = pool->uy[slot];
            v_in[2] = pool->uz[slot];
            osh_physics_moliere_scatter(v_in, ctx->w_scat, theta0, rng);
        }
    }
}

/* ---- Phase 5: nuclear interaction sampling ------------------------------- */

/*
 * Sample nuclear interactions over the current step.  Dispatches to
 * osh_nuclear_handler_step() when a compiled handler is present; otherwise
 * falls back to the legacy scalar Tripathi path.
 *
 * When pp elastic fires, w_scat is overwritten with incident_dir so that the
 * bent-path commit formula  q = p + h·u0 + tail·w_scat  reduces to a straight
 * endpoint  q = p + step_len·u0.  This is correct because elastic scattering
 * replaces MCS for this step.
 */
static void ion_step_nuclear(struct ion_step_ctx *ctx,
                             struct osh_transport_context const *transport_ctx,
                             struct osh_material_runtime const *material_rt,
                             struct osh_rng *rng) {
    (void) material_rt;

    if (!transport_ctx->params.nuclear_inelastic && !transport_ctx->params.nuclear_elastic)
        return;
    if (ctx->ds_gcm2 <= 0.0)
        return;

    if (transport_ctx->nuclear_handler) {
        osh_nuclear_handler_step(transport_ctx->nuclear_handler,
                                 ctx->e0,
                                 ctx->exit_energy,
                                 ctx->incident_dir,
                                 ctx->zone_material_idx,
                                 ctx->ds_gcm2,
                                 ctx->part,
                                 &transport_ctx->params,
                                 rng,
                                 &ctx->nuclear_event);

        /* Elastic overrides MCS: straight-line endpoint. */
        if (ctx->nuclear_event.kind == OSH_NUCLEAR_EVENT_ELASTIC_PP) {
            ctx->w_scat[0] = ctx->incident_dir[0];
            ctx->w_scat[1] = ctx->incident_dir[1];
            ctx->w_scat[2] = ctx->incident_dir[2];
        }
    } else if (transport_ctx->params.nuclear_inelastic) {
        /* Legacy scalar Tripathi path (no handler compiled). */
        double at;
        double e_per_nucleon;
        double sigma;
        double lambda;
        double p_survive;

        if (ctx->mat_z_over_a <= 0.0)
            return;
        at = ctx->mat_z_mean / ctx->mat_z_over_a;
        e_per_nucleon = ctx->e0 / ctx->a_proj;
        sigma = osh_nuclear_tripathi_sigma(ctx->part->z, ctx->part->a, ctx->mat_z_mean, at, e_per_nucleon);
        if (sigma <= 0.0)
            return;
        lambda = osh_nuclear_lambda_gcm2(at, sigma);
        p_survive = osh_nuclear_survival_prob(ctx->ds_gcm2, lambda);
        if (osh_rng_double(rng) > p_survive) {
            ctx->nuclear_event.kind = OSH_NUCLEAR_EVENT_ABSORB;
            ctx->nuclear_event.primary_energy = 0.0;
        }
    }
}

/* ---- Phase 6: commit ----------------------------------------------------- */

/**
 * Finalise the step: compute the bent-path exit position, score the step,
 * write the updated state back to the pool, and nudge past any boundary.
 *
 * Exit position (two-segment bent path):
 *   q = p + h·u0 + tail·w_scat
 *
 * For vacuum (h == boundary_ds, w_scat == u0) this reduces to q = p + ds·u0.
 * For material steps the hinge bends the path; the scoring chord p→q has an
 * O(θ₀² × step_len) ≈ 10⁻⁵ cm error in voxel DDA scoring, negligible
 * compared to voxel resolution.
 */
static enum osh_status ion_step_commit(struct ion_step_ctx const *ctx,
                                       struct osh_particle_pool *pool,
                                       size_t slot,
                                       struct osh_transport_context const *transport_ctx,
                                       struct osh_scoring_runtime *score_rt) {
    double qx;
    double qy;
    double qz;
    struct step st;
    enum osh_status rc;

    if (ctx->step_len <= 0.0) {
        /* Zero-length step: update direction only, then nudge if at a boundary */
        pool->ux[slot] = ctx->w_scat[0];
        pool->uy[slot] = ctx->w_scat[1];
        pool->uz[slot] = ctx->w_scat[2];
        if (ctx->hit_boundary) {
            pool->x[slot] += ctx->nudge_dir[0] * OSH_TRANSPORT_BOUNDARY_EPS;
            pool->y[slot] += ctx->nudge_dir[1] * OSH_TRANSPORT_BOUNDARY_EPS;
            pool->z[slot] += ctx->nudge_dir[2] * OSH_TRANSPORT_BOUNDARY_EPS;
        }
        return OSH_OK;
    }

    /* Bent-path exit position: p + h·u0 + tail·w_scat */
    qx = pool->x[slot] + pool->ux[slot] * ctx->h + ctx->w_scat[0] * ctx->tail_len;
    qy = pool->y[slot] + pool->uy[slot] * ctx->h + ctx->w_scat[1] * ctx->tail_len;
    qz = pool->z[slot] + pool->uz[slot] * ctx->h + ctx->w_scat[2] * ctx->tail_len;

    step_from_pool(
        &st, pool, slot, ctx->step_len, ctx->exit_energy, ctx->rho, (int) ctx->zone_material_idx, (int) ctx->zone_idx);
    st.q[0] = qx;
    st.q[1] = qy;
    st.q[2] = qz;
    st.q[3] = ctx->exit_energy;
    st.w[0] = ctx->w_scat[0];
    st.w[1] = ctx->w_scat[1];
    st.w[2] = ctx->w_scat[2];
    st.voxel_idx = ctx->voxel_idx;
    st.has_voxel = ctx->has_voxel;

    if (ctx->nuclear_event.kind == OSH_NUCLEAR_EVENT_ABSORB ||
        ctx->nuclear_event.kind == OSH_NUCLEAR_EVENT_ABRASION) {
        /* Primary destroyed by inelastic nuclear reaction.  st.de already
         * holds the ionisation energy deposited — do not modify it.
         * st.q[3] = 0 signals that the primary exits dead. */
        st.q[3] = 0.0;
    } else if (ctx->nuclear_event.kind == OSH_NUCLEAR_EVENT_ELASTIC_PP) {
        /* Primary survives with a new direction and energy from elastic scatter. */
        st.q[3] = ctx->nuclear_event.primary_energy;
        st.w[0] = ctx->nuclear_event.primary_dir[0];
        st.w[1] = ctx->nuclear_event.primary_dir[1];
        st.w[2] = ctx->nuclear_event.primary_dir[2];
    }

    rc = osh_scoring_score_step(score_rt, ctx->part, &st);
    if (rc != OSH_OK) {
        OSH_DIAG_ERRORF(transport_ctx->diag,
                        "transport: scoring rejected step rc=%d ds=%.17g p=(%.17g, %.17g, %.17g) "
                        "q=(%.17g, %.17g, %.17g)",
                        (int) rc,
                        st.ds,
                        st.p[0],
                        st.p[1],
                        st.p[2],
                        st.q[0],
                        st.q[1],
                        st.q[2]);
        return rc;
    }

    /* Write the updated particle state — always use the bent-path exit position */
    pool->x[slot] = qx;
    pool->y[slot] = qy;
    pool->z[slot] = qz;
    if (ctx->nuclear_event.kind == OSH_NUCLEAR_EVENT_ABSORB ||
        ctx->nuclear_event.kind == OSH_NUCLEAR_EVENT_ABRASION) {
        pool->e[slot] = 0.0;
        pool->ux[slot] = ctx->w_scat[0];
        pool->uy[slot] = ctx->w_scat[1];
        pool->uz[slot] = ctx->w_scat[2];
    } else if (ctx->nuclear_event.kind == OSH_NUCLEAR_EVENT_ELASTIC_PP) {
        pool->e[slot] = ctx->nuclear_event.primary_energy;
        pool->ux[slot] = ctx->nuclear_event.primary_dir[0];
        pool->uy[slot] = ctx->nuclear_event.primary_dir[1];
        pool->uz[slot] = ctx->nuclear_event.primary_dir[2];
    } else {
        pool->e[slot] = ctx->exit_energy;
        pool->ux[slot] = ctx->w_scat[0];
        pool->uy[slot] = ctx->w_scat[1];
        pool->uz[slot] = ctx->w_scat[2];
    }

    if (ctx->hit_boundary) {
        /* Nudge past the boundary in the leg direction that reached it:
         *   straight/boundary steps:   u0 (nudge_dir == initial direction)
         *   hinge-clipped scatter:     w_scat (nudge_dir set in phase 3) */
        pool->x[slot] += ctx->nudge_dir[0] * OSH_TRANSPORT_BOUNDARY_EPS;
        pool->y[slot] += ctx->nudge_dir[1] * OSH_TRANSPORT_BOUNDARY_EPS;
        pool->z[slot] += ctx->nudge_dir[2] * OSH_TRANSPORT_BOUNDARY_EPS;
    }

    return OSH_OK;
}

/* ---- Table helpers ------------------------------------------------------- */

static double cutoff_total_energy(struct osh_transport_params const *params,
                                  struct osh_material_runtime const *material_rt,
                                  struct particle const *part) {
    double a_proj;
    double cutoff_total;
    double cutoff_from_params;
    double cutoff_from_transport;

    a_proj = (part->a > 0u) ? (double) part->a : 1.0;
    cutoff_total = OSH_TRANSPORT_ION_EMIN_MEV_PER_U * a_proj;
    if (params) {
        cutoff_from_params = (double) params->tcut * a_proj;
        if (cutoff_from_params > cutoff_total) {
            cutoff_total = cutoff_from_params;
        }
    }
    cutoff_from_transport = material_rt->emin * a_proj;
    if (cutoff_from_transport > cutoff_total) {
        cutoff_total = cutoff_from_transport;
    }
    return cutoff_total;
}

static double energy_from_residual_range(struct osh_material_runtime const *material_rt,
                                         size_t material_idx,
                                         size_t projectile_idx,
                                         double residual_range) {
    /* TODO(simd): add a batch form alongside osh_material_runtime_range_lookup()
     * so the transport stepper can process residual-range -> energy inversion
     * for multiple lanes without scalar binary-searching each slot separately. */
    float const *range_col;
    size_t lo;
    size_t hi;
    size_t mid;
    double r_lo;
    double r_hi;
    double e_lo;
    double e_hi;
    double frac;

    range_col =
        material_rt->range_csda + (material_idx * material_rt->nprojectiles + projectile_idx) * material_rt->nenergy;

    if (residual_range <= (double) range_col[0]) {
        return material_rt->emin;
    }
    if (residual_range >= (double) range_col[material_rt->nenergy - 1u]) {
        return material_rt->emax;
    }

    lo = 0u;
    hi = material_rt->nenergy - 1u;
    while (hi - lo > 1u) {
        mid = lo + (hi - lo) / 2u;
        if ((double) range_col[mid] <= residual_range) {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    r_lo = (double) range_col[lo];
    r_hi = (double) range_col[hi];
    e_lo = energy_grid_value(material_rt, lo);
    e_hi = energy_grid_value(material_rt, hi);
    if (r_hi <= r_lo) {
        return e_lo;
    }
    frac = (residual_range - r_lo) / (r_hi - r_lo);
    return e_lo * (1.0 - frac) + e_hi * frac;
}

static double energy_grid_value(struct osh_material_runtime const *material_rt, size_t energy_idx) {
    return exp(material_rt->log_emin + (double) energy_idx / material_rt->inv_dlog);
}

static enum osh_status find_projectile_index(struct osh_material_runtime const *material_rt,
                                             struct particle const *part,
                                             size_t *projectile_idx_out) {
    unsigned int z_match;
    unsigned int a_match;
    size_t projectile_idx;

    if (!material_rt || !part || !projectile_idx_out) {
        return OSH_EINVAL;
    }
    z_match = part->z;
    a_match = part->a;
    if (part->pdg == 2212 && z_match == 0u && a_match == 0u) {
        z_match = 1u;
        a_match = 1u;
    }
    if (z_match == 0u) {
        return OSH_ENOTSUP;
    }
    projectile_idx = (size_t) (z_match - 1u);
    if (projectile_idx >= material_rt->nprojectiles) {
        return OSH_ENOTSUP;
    }
    if (material_rt->projectile_z[projectile_idx] != z_match) {
        return OSH_ESTATE;
    }
    if (a_match != 0u && material_rt->projectile_a[projectile_idx] != a_match) {
        /* Runtime columns are keyed primarily by Z; differing isotopes share the
         * representative projectile for that Z for now. */
    }
    *projectile_idx_out = projectile_idx;
    return OSH_OK;
}

static int is_blackhole_material(size_t material_idx) {
    return material_idx == OSH_MATERIAL_INDEX_BLACKHOLE;
}

static int is_vacuum_material(size_t material_idx) {
    return material_idx == OSH_MATERIAL_INDEX_VACUUM;
}

/**
 * @brief Build a scoring step from pool slot @p slot and step kinematics.
 *
 * @details
 * Populates all fields of @p st that can be derived from the pool SoA state
 * and the computed step kinematics.  The caller must still override st.q and
 * st.w with the actual bent-path endpoint and exit direction.
 */
static void step_from_pool(struct step *st,
                           struct osh_particle_pool const *pool,
                           size_t slot,
                           double step_len,
                           double exit_energy,
                           double rho,
                           int medium,
                           int zone) {
    st->p[0] = pool->x[slot];
    st->p[1] = pool->y[slot];
    st->p[2] = pool->z[slot];
    st->p[3] = pool->e[slot];
    /* st->q and st->w are overridden by the caller with bent-path values */
    st->q[0] = pool->x[slot] + pool->ux[slot] * step_len;
    st->q[1] = pool->y[slot] + pool->uy[slot] * step_len;
    st->q[2] = pool->z[slot] + pool->uz[slot] * step_len;
    st->q[3] = exit_energy;
    st->v[0] = pool->ux[slot];
    st->v[1] = pool->uy[slot];
    st->v[2] = pool->uz[slot];
    st->w[0] = pool->ux[slot];
    st->w[1] = pool->uy[slot];
    st->w[2] = pool->uz[slot];
    st->ds = step_len;
    st->de = pool->e[slot] - exit_energy;
    if (st->de < 0.0) {
        st->de = 0.0;
    }
    st->rho = rho;
    st->wt = pool->wt[slot];
    st->medium = medium;
    st->zone = zone;
    st->system = OSH_COORD_UNIVERSE;
    st->prim_idx = pool->prim_idx[slot];
    st->gen = pool->gen[slot];
}
