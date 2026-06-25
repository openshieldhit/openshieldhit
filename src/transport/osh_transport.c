#include "transport/osh_transport.h"

#include "common/osh_diag.h"
#include "transport/osh_neutron_pool.h"
#include "transport/osh_transport_ion.h"
#include "transport/osh_transport_neutron.h"
#include "transport/osh_transport_photon.h"
#include "transport/osh_transport_scheduler.h"

/* ---- Forward declarations ------------------------------------------------ */

static enum osh_status dispatch_transport_family(enum osh_transport_family family,
                                                 struct osh_transport_context *transport_ctx,
                                                 struct osh_beam_runtime *beam_rt,
                                                 struct osh_gemca_runtime const *geom_rt,
                                                 struct osh_material_runtime const *material_rt,
                                                 struct osh_scoring_runtime *score_rt);

/*
 * Orchestrator for the minimal transport run.
 *
 * Drives a scheduler loop over transport families (ions, neutrons, …).
 * dispatch_transport_family() is the natural parallelism boundary: it fully
 * drains one pool before returning, touching no scheduler state.  All
 * has_work bookkeeping lives here so that future multi-pool or GPU-offload
 * refactors only touch this function.
 *
 * has_work discipline:
 *   • Clear the family's own has_work BEFORE the dispatch call (avoids an
 *     accidental self-reschedule if the kernel ever pushes to its own pool).
 *   • After the call, set has_work for peer families based on what changed.
 *   • scheduler_next() does NOT clear has_work on return.
 *
 * Current scheduling policy:
 *   1. Ion pass        — processes all nstat primaries from beam_rt.
 *   2. Neutron pass    — drains neutron_pool if non-empty after the ion pass.
 *   3. Ion feedback    — reschedule ions only if nuclear_neutron_ion_feedback
 *                        is set AND the ion pool received (n,x) secondaries.
 *                        Not yet reachable: neutron transport does not push to
 *                        the ion pool yet (see params.nuclear_neutron_ion_feedback).
 */
enum osh_status osh_transport_run_minimal(struct osh_transport_context *transport_ctx,
                                          struct osh_beam_runtime *beam_rt,
                                          struct osh_gemca_runtime const *geom_rt,
                                          struct osh_material_runtime const *material_rt,
                                          struct osh_scoring_runtime *score_rt) {
    struct osh_transport_scheduler scheduler;
    enum osh_transport_family family;
    enum osh_status rc;
    size_t neutron_capacity;
    int neutron_enabled; /* cache: neutron pool is present for this run */

    if (!transport_ctx) {
        return OSH_EINVAL;
    }

    osh_transport_scheduler_reset(&scheduler);

    rc = osh_transport_scheduler_enable(&scheduler, OSH_TRANSPORT_FAMILY_ION);
    if (rc != OSH_OK) {
        return rc;
    }
    osh_transport_scheduler_set_has_work(&scheduler, OSH_TRANSPORT_FAMILY_ION, 1);

    /* Neutron family: enable when a pool is available; pool starts empty. */
    neutron_enabled = (transport_ctx->neutron_pool != NULL);
    if (neutron_enabled) {
        neutron_capacity = (transport_ctx->params.pool_capacity != 0u) ? transport_ctx->params.pool_capacity
                                                                       : (size_t) OSH_TRANSPORT_POOL_CAPACITY;
        if (transport_ctx->neutron_pool->capacity == 0u) {
            rc = osh_neutron_pool_init(transport_ctx->neutron_pool, neutron_capacity);
            if (rc != OSH_OK) {
                return rc;
            }
        } else {
            osh_neutron_pool_reset(transport_ctx->neutron_pool);
        }
        rc = osh_transport_scheduler_enable(&scheduler, OSH_TRANSPORT_FAMILY_NEUTRON);
        if (rc != OSH_OK) {
            return rc;
        }
        /* has_work stays 0 until the ion pass produces neutrons. */
    }

    while (osh_transport_scheduler_next(&scheduler, &family)) {
        /* Clear before dispatch so the family does not immediately re-schedule
         * itself; it will only run again if a feedback path sets has_work=1. */
        osh_transport_scheduler_set_has_work(&scheduler, family, 0);

        rc = dispatch_transport_family(family, transport_ctx, beam_rt, geom_rt, material_rt, score_rt);
        if (rc != OSH_OK) {
            return rc;
        }

        if (family == OSH_TRANSPORT_FAMILY_ION && neutron_enabled) {
            /* Neutrons produced during the ion pass are now in the pool. */
            osh_transport_scheduler_set_has_work(
                &scheduler, OSH_TRANSPORT_FAMILY_NEUTRON, (transport_ctx->neutron_pool->n > 0u) ? 1 : 0);
        }

        if (family == OSH_TRANSPORT_FAMILY_NEUTRON && transport_ctx->params.nuclear_neutron_ion_feedback) {
            /* TODO: replace 0 with (ion_pool->n > 0) once osh_transport_neutron.c
             * is wired to push (n,p)/(n,α)/compound ion secondaries to the ion pool. */
            osh_transport_scheduler_set_has_work(&scheduler, OSH_TRANSPORT_FAMILY_ION, 0);
        }
    }

    return OSH_OK;
}

/**
 * @brief Dispatch one transport family to its implementation module.
 *
 * @details
 * This switch is the future join point between the scheduler and the concrete
 * transport kernels.  Keeping the mapping here avoids spreading knowledge of
 * per-family source files through the rest of transport/.
 *
 * @param[in]     family     Scheduled family to transport.
 * @param[in]     beam_rt    Hot beam runtime for primary generation.
 * @param[in]     geom_rt    Compiled geometry runtime.
 * @param[in]     material_rt  Hot material runtime tables.
 * @param[in,out] score_rt     Scoring runtime.
 *
 * @returns OSH_OK on success, OSH_ENOTSUP for a family without an
 *          implementation, or another OSH_E* from the family kernel.
 */
static enum osh_status dispatch_transport_family(enum osh_transport_family family,
                                                 struct osh_transport_context *transport_ctx,
                                                 struct osh_beam_runtime *beam_rt,
                                                 struct osh_gemca_runtime const *geom_rt,
                                                 struct osh_material_runtime const *material_rt,
                                                 struct osh_scoring_runtime *score_rt) {
    switch (family) {
    case OSH_TRANSPORT_FAMILY_ION:
        return osh_transport_ion_run_minimal(transport_ctx, beam_rt, geom_rt, material_rt, score_rt);
    case OSH_TRANSPORT_FAMILY_NEUTRON:
        return osh_transport_neutron_run(transport_ctx, beam_rt, geom_rt, material_rt, score_rt);
    case OSH_TRANSPORT_FAMILY_PHOTON:
        return osh_transport_photon_run_minimal(transport_ctx, beam_rt, geom_rt, material_rt, score_rt);
    case OSH_TRANSPORT_FAMILY_ELECTRON:
    case OSH_TRANSPORT_FAMILY_COUNT:
        break;
    }

    OSH_DIAG_ERRORF(transport_ctx ? transport_ctx->diag : NULL,
                    "transport: %s transport is not implemented",
                    osh_transport_family_name(family));
    return OSH_ENOTSUP;
}
