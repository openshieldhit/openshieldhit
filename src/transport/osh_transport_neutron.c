#include "transport/osh_transport_neutron.h"

#include "common/osh_diag.h"
#include "transport/osh_transport.h"

enum osh_status osh_transport_neutron_run_minimal(struct osh_transport_context *transport_ctx,
                                                  struct osh_beam_runtime *beam_rt,
                                                  struct osh_gemca_runtime const *geom_rt,
                                                  struct osh_material_runtime const *material_rt,
                                                  struct osh_scoring_runtime *score_rt) {
    (void) beam_rt;
    (void) geom_rt;
    (void) material_rt;
    (void) score_rt;
    OSH_DIAG_ERRORF(
        transport_ctx ? transport_ctx->diag : NULL, "%s", "transport: neutron transport is not implemented");
    return OSH_ENOTSUP;
}
