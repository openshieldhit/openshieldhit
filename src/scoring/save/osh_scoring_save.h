#ifndef OSH_SCORING_SAVE_H
#define OSH_SCORING_SAVE_H

#include "openshieldhit/status.h"
#include "scoring/osh_scoring.h"
#include "scoring/runtime/osh_scoring_runtime.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Save all compiled outputs described by the workspace/runtime pair.
 *
 * @details
 * Dispatches to format-specific writers. BDO outputs embed @p nstat in the
 * file; ASCII outputs normalise each value by dividing by @p nstat.
 *
 * @param[in] ws     Scoring workspace with output metadata and file paths.
 * @param[in] rt     Compiled scoring runtime with postprocessed accumulators.
 * @param[in] nstat  Actual number of primary particles simulated; must be > 0.
 *
 * @returns OSH_OK on success, OSH_ENOTSUP for an unsupported format, or
 *          another OSH_E* on I/O error.
 */
enum osh_status osh_scoring_save(struct osh_scoring_workspace const *ws,
                                 struct osh_scoring_runtime const *rt,
                                 unsigned long long nstat);

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_SAVE_H */
