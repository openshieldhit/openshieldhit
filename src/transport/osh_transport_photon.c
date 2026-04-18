#include "transport/osh_transport_photon.h"

#include "common/osh_logger.h"
#include "transport/osh_transport.h"

enum osh_status osh_transport_photon_run_minimal(struct osh_transport_context *transport_ctx,
                                                 struct osh_beam_runtime *beam_rt,
                                                 struct osh_gemca_runtime const *geom_rt,
                                                 struct osh_material_runtime const *tables,
                                                 struct osh_scoring_runtime *scoring) {
    (void) beam_rt;
    (void) geom_rt;
    (void) tables;
    (void) scoring;
    OSH_DIAG_ERRORF(transport_ctx ? transport_ctx->diag : NULL, "%s", "transport: photon transport is not implemented");
    return OSH_ENOTSUP;
}
