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

/**
 * @brief Save a selected subset of outputs (G2), or all of them.
 *
 * @details
 * Same per-output writers as @ref osh_scoring_save, but restricted to the output
 * indices in @p want.  @p want == NULL writes every output (the @ref
 * osh_scoring_save behaviour); otherwise each entry must be a valid output index.
 * This is the granularity a mid-run snapshot uses to dump only the pages backing
 * the outputs of interest.
 *
 * @param[in] ws      Scoring workspace with output metadata and file paths.
 * @param[in] rt      Compiled scoring runtime with postprocessed accumulators.
 * @param[in] nstat   Actual number of primary particles simulated; must be > 0.
 * @param[in] want    Output indices to save, or NULL for all.
 * @param[in] n_want  Number of entries in @p want (ignored when @p want is NULL).
 *
 * @returns OSH_OK on success, OSH_EINVAL on a bad argument or out-of-range index,
 *          OSH_ESTATE on a workspace/runtime output-count mismatch, OSH_ENOTSUP
 *          for an unsupported format, or another OSH_E* on I/O error.
 */
enum osh_status osh_scoring_save_outputs(struct osh_scoring_workspace const *ws,
                                         struct osh_scoring_runtime const *rt,
                                         unsigned long long nstat,
                                         size_t const *want,
                                         size_t n_want);

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_SAVE_H */
