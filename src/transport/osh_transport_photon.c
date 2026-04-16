#include "transport/osh_transport_photon.h"

#include "common/osh_logger.h"

enum osh_status osh_transport_photon_run_minimal(struct osh_transport_context *transport_ctx,
                                                 struct osh_beam_runtime *beam_rt,
                                                 struct osh_gemca_runtime const *geom_rt,
                                                 struct osh_material_runtime const *tables,
                                                 struct osh_scoring_runtime *scoring) {
    (void) transport_ctx;
    (void) beam_rt;
    (void) geom_rt;
    (void) tables;
    (void) scoring;
    osh_error("%s", "transport: photon transport is not implemented");
    return OSH_ENOTSUP;
}
