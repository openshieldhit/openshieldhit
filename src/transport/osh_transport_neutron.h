#ifndef OSH_TRANSPORT_NEUTRON_H
#define OSH_TRANSPORT_NEUTRON_H

#include "common/osh_rc.h"

#ifdef __cplusplus
extern "C" {
#endif

struct beam_workspace;
struct gemca_runtime;
struct osh_material_runtime;
struct osh_scoring_runtime;

/**
 * @brief Neutron transport loop (stub — not yet implemented).
 * @return OSH_ENOTSUP always.
 */
enum osh_status osh_transport_neutron_run_minimal(struct beam_workspace const *beam,
                                                  struct gemca_runtime const *geom_rt,
                                                  struct osh_material_runtime const *tables,
                                                  struct osh_scoring_runtime *scoring);

#ifdef __cplusplus
}
#endif

#endif /* OSH_TRANSPORT_NEUTRON_H */
