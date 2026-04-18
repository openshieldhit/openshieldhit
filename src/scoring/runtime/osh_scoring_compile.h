#ifndef OSH_SCORING_COMPILE_H
#define OSH_SCORING_COMPILE_H

#include "openshieldhit/logger.h"
#include "openshieldhit/status.h"
#include "scoring/osh_scoring.h"
#include "scoring/runtime/osh_scoring_runtime.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Compile parsed scorer configuration into scorer runtime objects.
 *
 * @details
 * This resolves raw `detect.dat` names to dense indices, computes geometry bin
 * counts, groups pages by shared geometry, and allocates page-local scoring
 * buffers. Geometry traversal itself is not built here yet; this is the first
 * compiled bridge from readable parser structs to simulation-ready scorer
 * runtime objects.
 */
enum osh_status osh_scoring_compile(struct osh_scoring_workspace const *ws,
                                    struct osh_diag_sink const *diag,
                                    struct osh_scoring_runtime *rt);

/**
 * @brief Release a compiled scorer runtime and all owned memory.
 */
void osh_scoring_runtime_free(struct osh_scoring_runtime *rt);

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_COMPILE_H */
