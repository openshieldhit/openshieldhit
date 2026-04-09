#ifndef OSH_SCORING_POSTPROCESS_H
#define OSH_SCORING_POSTPROCESS_H

#include "common/osh_rc.h"
#include "scoring/runtime/osh_scoring_runtime.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Apply postprocessing to accumulated scoring page data.
 *
 * @details
 * This is a separate cold-path phase between raw accumulation and save. The
 * current implementation is intentionally minimal: for the currently
 * implemented page kinds, postprocessing is a no-op.
 */
enum osh_status osh_scoring_postprocess(struct osh_scoring_runtime *rt);

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_POSTPROCESS_H */
