#ifndef OSH_TRANSPORT_PHOTON_H
#define OSH_TRANSPORT_PHOTON_H

#include "common/osh_rc.h"

#ifdef __cplusplus
extern "C" {
#endif

struct beam_workspace;
struct gemca_runtime;
struct material_workspace;
struct osh_material_runtime;
struct osh_scoring_runtime;

/**
 * @brief Photon transport loop (stub — not yet implemented).
 * @return OSH_ENOTSUP always.
 */
enum osh_status osh_transport_photon_run_minimal(struct beam_workspace const *beam,
                                                 struct gemca_runtime const *geom_rt,
                                                 struct material_workspace const *materials,
                                                 struct osh_material_runtime const *tables,
                                                 struct osh_scoring_runtime *scoring);

#ifdef __cplusplus
}
#endif

#endif /* OSH_TRANSPORT_PHOTON_H */
