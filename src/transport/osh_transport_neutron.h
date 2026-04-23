#ifndef OSH_TRANSPORT_NEUTRON_H
#define OSH_TRANSPORT_NEUTRON_H

#include "openshieldhit/status.h"

#ifdef __cplusplus
extern "C" {
#endif

struct osh_transport_context;
struct osh_beam_runtime;
struct osh_gemca_runtime;
struct osh_material_runtime;
struct osh_scoring_runtime;

/**
 * @brief Neutron transport loop (stub — not yet implemented).
 * @return OSH_ENOTSUP always.
 */
enum osh_status osh_transport_neutron_run_minimal(struct osh_transport_context *transport_ctx,
                                                  struct osh_beam_runtime *beam_rt,
                                                  struct osh_gemca_runtime const *geom_rt,
                                                  struct osh_material_runtime const *material_rt,
                                                  struct osh_scoring_runtime *score_rt);

#ifdef __cplusplus
}
#endif

#endif /* OSH_TRANSPORT_NEUTRON_H */
