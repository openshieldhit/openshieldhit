#include <stdlib.h>
#include <string.h>

#include "beam/osh_beam.h"
#include "beam/runtime/osh_beam_runtime.h"
#include "gemca/osh_gemca2.h"
#include "gemca/runtime/osh_gemca_runtime.h"
#include "material/runtime/osh_material_compile.h"
#include "openshieldhit/beam_defs.h"
#include "openshieldhit/logger.h"
#include "openshieldhit/simulation.h"
#include "scoring/runtime/osh_scoring_compile.h"
#include "scoring/runtime/osh_scoring_postprocess.h"
#include "scoring/save/osh_scoring_save.h"
#include "transport/osh_transport.h"

/* ---- Private definition of the opaque handle ----------------------------- */

struct osh_simulation {
    struct osh_beam_workspace const *beam;
    struct osh_scoring_workspace const *scoring;

    struct osh_beam_runtime *beam_rt;
    struct osh_gemca_runtime geom_rt;
    struct osh_material_runtime transport_tables;
    struct osh_scoring_runtime scoring_runtime;
    struct osh_transport_context transport_ctx;
};

/* ---- Lifecycle ----------------------------------------------------------- */

enum osh_status osh_simulation_create(struct osh_beam_workspace *beam,
                                      struct osh_geometry_workspace *geo,
                                      struct osh_material_workspace *mat,
                                      struct osh_scoring_workspace *scoring,
                                      struct osh_simulation **sim_out) {
    struct osh_simulation *sim;
    struct osh_gemca_prepared *gemca;
    unsigned int z_max;
    size_t iz;
    enum osh_status rc;

    if (!beam || !geo || !mat || !scoring || !sim_out) {
        return OSH_EINVAL;
    }
    *sim_out = NULL;
    if (!geo->prepared) {
        osh_error("%s", "simulation: geometry workspace has not been prepared");
        return OSH_ESTATE;
    }
    if (!beam->prepared) {
        osh_error("%s", "simulation: beam workspace has not been prepared");
        return OSH_ESTATE;
    }
    if (!beam->has_primary) {
        osh_error("%s", "simulation: beam workspace has no primary particle defined");
        return OSH_EINVAL;
    }

    sim = (struct osh_simulation *) calloc(1, sizeof(*sim));
    if (!sim) {
        return OSH_ENOMEM;
    }
    sim->beam = beam;
    sim->scoring = scoring;

    /* ---- 1. Zone → material index resolution ----------------------------- */

    gemca = geo->prepared;
    for (iz = 0u; iz < gemca->nzones; ++iz) {
        struct zone *z = gemca->zones[iz];
        struct osh_material const *m;

        if (!z->material_name) {
            osh_error("simulation: zone '%s' has no material assigned", z->name ? z->name : "(unnamed)");
            rc = OSH_EPARSE;
            goto fail;
        }
        m = osh_material_by_name(mat, z->material_name);
        if (!m) {
            osh_error(
                "simulation: zone '%s': unknown material '%s'", z->name ? z->name : "(unnamed)", z->material_name);
            rc = OSH_EPARSE;
            goto fail;
        }
        z->material_idx = m->index;
    }

    /* ---- 2. Geometry runtime -------------------------------------------- */

    rc = osh_gemca_compile(gemca, &sim->geom_rt);
    if (rc != OSH_OK) {
        osh_error("%s", "simulation: failed to compile geometry runtime");
        goto fail;
    }

    /* ---- 3. Transport tables --------------------------------------------- */

    z_max = (beam->primary.z > 0u) ? (unsigned int) beam->primary.z : 1u;
    rc = osh_material_compile(mat, z_max, &sim->transport_tables);
    if (rc != OSH_OK) {
        osh_error("%s", "simulation: failed to compile transport tables");
        goto fail;
    }

    /* ---- 4. Scoring runtime --------------------------------------------- */

    rc = osh_scoring_compile(scoring, &sim->scoring_runtime);
    if (rc != OSH_OK) {
        osh_error("%s", "simulation: failed to compile scoring runtime");
        goto fail;
    }

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

    /* ---- 6. Beam runtime ------------------------------------------------- */

    rc = osh_beam_compile(beam, &sim->beam_rt);
    if (rc != OSH_OK) {
        osh_error("%s", "simulation: failed to initialise beam runtime");
        goto fail;
    }

    *sim_out = sim;
    return OSH_OK;

fail:
    osh_simulation_free(sim);
    return rc;
}

enum osh_status osh_simulation_run(struct osh_simulation *sim, char const *out_dir) {
    struct osh_scoring_save_request save_req;
    enum osh_status rc;

    if (!sim || !out_dir) {
        return OSH_EINVAL;
    }

    rc = osh_transport_run_minimal(
        &sim->transport_ctx, sim->beam_rt, &sim->geom_rt, &sim->transport_tables, &sim->scoring_runtime);
    if (rc != OSH_OK) {
        osh_error("%s", "simulation: transport failed");
        return rc;
    }

    rc = osh_scoring_postprocess(&sim->scoring_runtime);
    if (rc != OSH_OK) {
        osh_error("%s", "simulation: scoring postprocess failed");
        return rc;
    }

    memset(&save_req, 0, sizeof(save_req));
    save_req.out_dir = out_dir;
    save_req.ws = sim->scoring;
    save_req.rt = &sim->scoring_runtime;
    save_req.nstat = sim->beam->nstat;
    save_req.has_nstat = 1;

    rc = osh_scoring_save(&save_req);
    if (rc != OSH_OK) {
        osh_error("%s", "simulation: scoring save failed");
        return rc;
    }

    return OSH_OK;
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
