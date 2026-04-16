#ifndef OSH_TRANSPORT_PHOTON_H
#define OSH_TRANSPORT_PHOTON_H

#include "common/osh_rc.h"

#ifdef __cplusplus
extern "C" {
#endif

struct osh_transport_context;
struct osh_beam_runtime;
struct gemca_runtime;
struct osh_material_runtime;
struct osh_scoring_runtime;

/**
 * @brief Photon transport loop (stub — not yet implemented).
 * @return OSH_ENOTSUP always.
 */
enum osh_status osh_transport_photon_run_minimal(struct osh_transport_context *transport_ctx,
                                                 struct osh_beam_runtime *beam_rt,
                                                 struct gemca_runtime const *geom_rt,
                                                 struct osh_material_runtime const *tables,
                                                 struct osh_scoring_runtime *scoring);

#ifdef __cplusplus
}
#endif

#endif /* OSH_TRANSPORT_PHOTON_H */
