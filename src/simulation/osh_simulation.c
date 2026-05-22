#include <stdlib.h>
#include <string.h>

#include "beam/osh_beam.h"
#include "beam/runtime/osh_beam_runtime.h"
#include "common/osh_diag.h"
#include "gemca/osh_gemca2.h"
#include "gemca/osh_gemca2_defines.h"
#include "gemca/runtime/osh_gemca_runtime.h"
#include "material/osh_material.h"
#include "material/runtime/osh_material_compile.h"
#include "openshieldhit/beam_defs.h"
#include "openshieldhit/simulation.h"
#include "scoring/runtime/osh_scoring_compile.h"
#include "scoring/runtime/osh_scoring_postprocess.h"
#include "scoring/save/osh_scoring_save.h"
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

    unsigned long long requested_nstat;
    unsigned long long completed_nstat;
    struct osh_results results;
};

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
    sim->transport_ctx.params.rndseed = beam->rndseed;
    sim->transport_ctx.params.rndoffset = beam->rndoffset;

    switch (beam->scatter) {
    case OSH_BEAM_MSCAT_GAUSS:
        sim->transport_ctx.params.mcs_mode = OSH_TRANSPORT_MCS_GAUSSIAN;
        break;
    case OSH_BEAM_MSCAT_MOLIERE:
        sim->transport_ctx.params.mcs_mode = OSH_TRANSPORT_MCS_MOLIERE;
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
    default:
        sim->transport_ctx.params.straggling_mode = OSH_TRANSPORT_STRAGGLING_OFF;
        break;
    }
    sim->transport_ctx.params.nuclear = beam->nuclear ? 1 : 0;
    sim->transport_ctx.diag = diag;

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

enum osh_status osh_simulation_run(struct osh_simulation *sim) {
    enum osh_status rc;

    if (!sim) {
        return OSH_EINVAL;
    }

    rc = osh_transport_run_minimal(
        &sim->transport_ctx, sim->beam_rt, &sim->geom_rt, &sim->transport_tables, &sim->scoring_runtime);
    if (rc != OSH_OK) {
        OSH_DIAG_ERRORF(sim->diag, "%s", "simulation: transport failed");
        return rc;
    }

    sim->completed_nstat = (unsigned long long) sim->transport_ctx.params.nstat;

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
    osh_scoring_runtime_free(&sim->scoring_runtime);
    osh_gemca_runtime_free(&sim->geom_rt);
    osh_material_runtime_free(&sim->transport_tables);
    osh_beam_runtime_free(sim->beam_rt);
    free(sim);
    return OSH_OK;
}
