#ifndef OSH_SCORING_MULTIFORMAT_H
#define OSH_SCORING_MULTIFORMAT_H

#include "openshieldhit/status.h"

#ifdef __cplusplus
extern "C" {
#endif

struct osh_scoring_workspace;
struct osh_scoring_runtime;
struct osh_diag_sink;

/**
 * @brief Expand a multi-format `Output` block into one runtime output per format.
 *
 * @details
 * Phases 1-6 of @ref osh_scoring_compile compile each cold `Output` block into
 * exactly one runtime output (runtime index == cold index) carrying the block's
 * first requested format and its filename as written.  This phase adds one extra
 * runtime output per additional format the block requests, each sharing the
 * block's `page_indices` — a cheap `size_t` copy; nothing in `rt->pages[]` is
 * duplicated, so scoring memory is independent of the format count (issue #308).
 *
 * With more than one format the block's `Filename` becomes a stem and each target
 * gets a canonical extension (or an explicit per-target name); a single format
 * keeps the filename as written, so an existing input file is unaffected.  RTDOSE
 * targets are reduced to a single page (see the single-format vs mixed-block rules
 * documented in `docs/user/detect.dat.md`).
 *
 * On any error the runtime is left in a state that @ref osh_scoring_runtime_free
 * can release (the caller owns that cleanup); this function never frees @p rt.
 *
 * @param[in]     ws    Cold scoring workspace (source of the format lists).
 * @param[in]     diag  Diagnostic sink for rejection messages (may be NULL).
 * @param[in,out] rt    Runtime built by phases 1-6; grown in place.
 *
 * @returns OSH_OK on success, or OSH_ENOMEM / OSH_EINVAL / OSH_ENOTSUP with a
 *          diagnostic on a rejected block.
 */
enum osh_status osh_scoring_expand_multiformat_outputs(struct osh_scoring_workspace const *ws,
                                                       struct osh_diag_sink const *diag,
                                                       struct osh_scoring_runtime *rt);

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_MULTIFORMAT_H */
