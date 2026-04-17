#ifndef OSH_SCORING_SAVE_H
#define OSH_SCORING_SAVE_H

#include "openshieldhit/status.h"
#include "scoring/osh_scoring.h"
#include "scoring/runtime/osh_scoring_runtime.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Cold-path save request for finalized scoring output.
 *
 * @details
 * Save stays outside `runtime/` because it needs workspace metadata, file
 * format policy, and disk I/O. The request object keeps that cold-path state
 * explicit.
 *
 * TODO: keep this request trivially hand-offable to a background CPU writer
 * thread once asynchronous saving is introduced. The writer should consume
 * finalized runtime buffers and must not depend on transport internals.
 */
struct osh_scoring_save_request {
    char const *out_dir;
    struct osh_scoring_workspace const *ws;
    struct osh_scoring_runtime const *rt;
    unsigned long long nstat;
    char has_nstat;
};

/**
 * @brief Save all compiled outputs described by the request.
 *
 * @details
 * Dispatches to format-specific writers in `src/scoring/save/`.
 */
enum osh_status osh_scoring_save(struct osh_scoring_save_request const *req);

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_SAVE_H */
