#ifndef OSH_TRANSPORT_ION_H
#define OSH_TRANSPORT_ION_H

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
 * @brief Wavefront CSDA transport loop for ion primaries.
 *
 * @details
 * Transports all params->nstat ion primaries to termination using a pool-based
 * (BFS) wavefront loop.  Physics: CSDA energy loss, Highland/Molière MCS
 * (random-hinge method), Bohr Gaussian energy straggling.  No nuclear
 * interactions or secondaries.
 *
 * The random-hinge treatment follows the fast proton-transport approach
 * described by Fippel and Soukup.
 *
 * @par References
 * Fippel M, Soukup M. A Monte Carlo dose calculation algorithm for proton
 * therapy. Med Phys. 2004;31(8):2263-2273. doi:10.1118/1.1769631.
 *
 * The caller must invoke osh_gemca_runtime_setup() before and
 * osh_gemca_runtime_free() after this function.
 *
 * @sa osh_transport_run_minimal() — public dispatcher that calls this.
 */
enum osh_status osh_transport_ion_run_minimal(struct osh_transport_context *transport_ctx,
                                              struct osh_beam_runtime *beam_rt,
                                              struct osh_gemca_runtime const *geom_rt,
                                              struct osh_material_runtime const *tables,
                                              struct osh_scoring_runtime *scoring);

#ifdef __cplusplus
}
#endif

#endif /* OSH_TRANSPORT_ION_H */
