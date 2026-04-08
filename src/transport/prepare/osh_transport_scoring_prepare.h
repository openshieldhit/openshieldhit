#ifndef OSH_TRANSPORT_SCORING_PREPARE_H
#define OSH_TRANSPORT_SCORING_PREPARE_H

#include "common/osh_rc.h"
#include "scoring/osh_scoring.h"
#include "transport/prepare/osh_transport_scoring_runtime.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Compile parsed scorer configuration into transport-facing runtime objects.
 *
 * @details
 * This resolves raw `detect.dat` names to dense indices, computes geometry bin
 * counts, groups pages by shared geometry, and allocates page-local scoring
 * buffers. Geometry traversal itself is not built here yet; this is the first
 * compiled bridge from readable parser structs to hot-path scorer runtime.
 */
enum osh_status osh_transport_scoring_prepare(struct osh_scoring_workspace const *ws,
                                              struct osh_transport_scoring_runtime *rt);

/**
 * @brief Release a compiled scorer runtime and all owned memory.
 */
void osh_transport_scoring_runtime_free(struct osh_transport_scoring_runtime *rt);

#ifdef __cplusplus
}
#endif

#endif /* OSH_TRANSPORT_SCORING_PREPARE_H */
