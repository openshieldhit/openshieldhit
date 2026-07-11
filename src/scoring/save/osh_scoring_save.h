#ifndef OSH_SCORING_SAVE_H
#define OSH_SCORING_SAVE_H

#include "openshieldhit/status.h"
#include "scoring/osh_scoring.h"
#include "scoring/runtime/osh_scoring_runtime.h"

#ifdef __cplusplus
extern "C" {
#endif

struct osh_diag_sink;

/**
 * @brief Save all compiled outputs described by the workspace/runtime pair.
 *
 * @details
 * Dispatches to format-specific writers. BDO outputs embed @p nstat in the
 * file; ASCII outputs normalise each value by dividing by @p nstat.
 *
 * Saving is **best-effort** (issue #308): a target that fails to write does not
 * stop the remaining targets — every output is attempted, and the *first* error
 * status observed is returned. A non-OK return therefore means "at least one
 * target failed" (which one is not distinguished here). This lets a multi-format
 * block still emit the targets that can be written when one path is bad, and lets
 * the process exit status reflect that a write failed.
 *
 * @param[in] ws     Scoring workspace with output metadata and file paths.
 * @param[in] rt     Compiled scoring runtime with postprocessed accumulators.
 * @param[in] nstat  Actual number of primary particles simulated; must be > 0.
 *
 * @returns OSH_OK when every target was written; otherwise the first error
 *          observed (OSH_ENOTSUP for an unsupported format, or another OSH_E* on
 *          I/O error), after still attempting the remaining targets.
 */
enum osh_status osh_scoring_save(struct osh_scoring_workspace const *ws,
                                 struct osh_scoring_runtime const *rt,
                                 unsigned long long nstat);

/**
 * @brief Save all compiled outputs, reporting each per-target failure to @p diag.
 *
 * @details
 * Identical to @ref osh_scoring_save (same best-effort semantics and return
 * value), but when a target fails to write, a warning naming the file, the
 * format, and the status is emitted to @p diag before continuing with the
 * remaining targets. This is the entry point the run driver uses so a partial
 * write tells the user *which* output failed, not just that "saving failed".
 *
 * @param[in] ws     Scoring workspace with output metadata and file paths.
 * @param[in] rt     Compiled scoring runtime with postprocessed accumulators.
 * @param[in] nstat  Actual number of primary particles simulated; must be > 0.
 * @param[in] diag   Diagnostic sink for per-target warnings (may be NULL).
 *
 * @returns Same as @ref osh_scoring_save.
 */
enum osh_status osh_scoring_save_diag(struct osh_scoring_workspace const *ws,
                                      struct osh_scoring_runtime const *rt,
                                      unsigned long long nstat,
                                      struct osh_diag_sink const *diag);

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
 * Writing is **best-effort** (issue #308): a failing target does not stop the
 * others; every requested target is attempted and the first error observed is
 * returned. Argument validation (NULL @p ws / @p rt, @p nstat == 0, an
 * out-of-range @p want index, or an output-count mismatch) is still fail-fast and
 * short-circuits before any target is written.
 *
 * @param[in] ws      Scoring workspace with output metadata and file paths.
 * @param[in] rt      Compiled scoring runtime with postprocessed accumulators.
 * @param[in] nstat   Actual number of primary particles simulated; must be > 0.
 * @param[in] want    Output indices to save, or NULL for all.
 * @param[in] n_want  Number of entries in @p want (ignored when @p want is NULL).
 *
 * @returns OSH_OK when every target was written; OSH_EINVAL on a bad argument or
 *          out-of-range index, OSH_ESTATE on a workspace/runtime output-count
 *          mismatch (both checked before writing), otherwise the first per-target
 *          error observed (OSH_ENOTSUP for an unsupported format, or another
 *          OSH_E* on I/O error) after attempting the remaining targets.
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
