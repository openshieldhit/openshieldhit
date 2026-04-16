#include "transport/osh_transport.h"

#include "common/osh_logger.h"
#include "transport/osh_transport_ion.h"
#include "transport/osh_transport_neutron.h"
#include "transport/osh_transport_photon.h"
#include "transport/osh_transport_scheduler.h"

/* ---- Forward declarations ------------------------------------------------ */

static enum osh_status dispatch_transport_family(enum osh_transport_family family,
                                                 struct beam_workspace const *beam,
                                                 struct gemca_runtime const *geom_rt,
                                                 struct osh_material_runtime const *tables,
                                                 struct osh_scoring_runtime *scoring);

/*
 * Particle-type dispatcher for osh_transport_run_minimal().
 *
 * The current public API still exposes a single "run minimal transport"
 * entry point, but the intended long-term owner of orchestration is a
 * transport-family scheduler with one queue or pool per particle class.
 *
 * The transport modules are split by family:
 *
 *   ions      -> osh_transport_ion.c      (implemented)
 *   neutrons  -> osh_transport_neutron.c  (stub)
 *   photons   -> osh_transport_photon.c   (stub)
 *   electrons -> future module
 *
 * Nuclear reactions and other in-transport physics will eventually feed
 * secondaries into scheduler-owned family queues.  For now only the ion
 * family is enabled and seeded with primary work, so the function preserves
 * the current ion-only behaviour exactly.
 */
enum osh_status osh_transport_run_minimal(struct beam_workspace const *beam,
                                          struct gemca_runtime const *geom_rt,
                                          struct osh_material_runtime const *tables,
                                          struct osh_scoring_runtime *scoring) {
    struct osh_transport_scheduler scheduler;
    enum osh_transport_family family;
    enum osh_status rc;

    osh_transport_scheduler_reset(&scheduler);

    rc = osh_transport_scheduler_enable(&scheduler, OSH_TRANSPORT_FAMILY_ION);
    if (rc != OSH_OK) {
        return rc;
    }
    osh_transport_scheduler_set_has_work(&scheduler, OSH_TRANSPORT_FAMILY_ION, 1);

    if (!osh_transport_scheduler_next(&scheduler, &family)) {
        return OSH_OK;
    }

    return dispatch_transport_family(family, beam, geom_rt, tables, scoring);
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
 * @param[in]     beam       Beam workspace for the run.
 * @param[in]     geom_rt    Compiled geometry runtime.
 * @param[in]     tables     Hot material runtime tables.
 * @param[in,out] scoring    Scoring runtime.
 *
 * @returns OSH_OK on success, OSH_ENOTSUP for a family without an
 *          implementation, or another OSH_E* from the family kernel.
 */
static enum osh_status dispatch_transport_family(enum osh_transport_family family,
                                                 struct beam_workspace const *beam,
                                                 struct gemca_runtime const *geom_rt,
                                                 struct osh_material_runtime const *tables,
                                                 struct osh_scoring_runtime *scoring) {
    switch (family) {
    case OSH_TRANSPORT_FAMILY_ION:
        return osh_transport_ion_run_minimal(beam, geom_rt, tables, scoring);
    case OSH_TRANSPORT_FAMILY_NEUTRON:
        return osh_transport_neutron_run_minimal(beam, geom_rt, tables, scoring);
    case OSH_TRANSPORT_FAMILY_PHOTON:
        return osh_transport_photon_run_minimal(beam, geom_rt, tables, scoring);
    case OSH_TRANSPORT_FAMILY_ELECTRON:
    case OSH_TRANSPORT_FAMILY_COUNT:
        break;
    }

    osh_error("transport: %s transport is not implemented", osh_transport_family_name(family));
    return OSH_ENOTSUP;
}
