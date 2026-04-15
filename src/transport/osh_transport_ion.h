#ifndef OSH_TRANSPORT_ION_H
#define OSH_TRANSPORT_ION_H

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
 * @brief Wavefront CSDA transport loop for ion primaries.
 *
 * @details
 * Transports all beam->nstat ion primaries to termination using a pool-based
 * (BFS) wavefront loop.  Physics: CSDA energy loss, Highland/Molière MCS
 * (random-hinge method), Bohr Gaussian energy straggling.  No nuclear
 * interactions or secondaries.
 *
 * The caller must invoke osh_gemca_runtime_setup() before and
 * osh_gemca_runtime_free() after this function.
 *
 * @sa osh_transport_run_minimal() — public dispatcher that calls this.
 */
enum osh_status osh_transport_ion_run_minimal(struct beam_workspace const *beam,
                                              struct gemca_runtime const *geom_rt,
                                              struct material_workspace const *materials,
                                              struct osh_material_runtime const *tables,
                                              struct osh_scoring_runtime *scoring);

#ifdef __cplusplus
}
#endif

#endif /* OSH_TRANSPORT_ION_H */
