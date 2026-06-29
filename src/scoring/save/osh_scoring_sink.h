#ifndef OSH_SCORING_SINK_H
#define OSH_SCORING_SINK_H

#include <stddef.h>

#include "openshieldhit/status.h"
#include "scoring/osh_scoring.h"
#include "scoring/runtime/osh_scoring_runtime.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Abstract destination for a (postprocessed) scoring snapshot — G1.
 *
 * @details
 * The snapshot path writes to this seam, not a hard-coded file path, so the same
 * primitive serves the native CLI (a file sink), POSIX signals, and the planned
 * WASM target (an in-memory buffer for postMessage).  @c save receives the
 * presentation runtime, the completed primary count to normalise by, and an
 * optional output selector @p want (a list of output indices; NULL = all
 * outputs) — G2.
 */
struct osh_scoring_sink {
    enum osh_status (*save)(void *ctx,
                            struct osh_scoring_runtime const *rt,
                            unsigned long long completed_nstat,
                            size_t const *want,
                            size_t n_want);
    void *ctx;
};

/**
 * @brief Caller-owned context for the native file sink.
 *
 * Holds the workspace whose output metadata / file paths the writers consume.
 * Must outlive the @ref osh_scoring_sink initialised from it.
 */
struct osh_scoring_file_sink {
    struct osh_scoring_workspace const *ws;
};

/**
 * @brief Initialise a file sink that writes via @ref osh_scoring_save_outputs.
 *
 * @param[out] fs   Caller-owned context (stores @p ws; must outlive @p out).
 * @param[in]  ws   Scoring workspace.
 * @param[out] out  Sink to populate.
 * @returns OSH_OK, or OSH_EINVAL on a NULL argument.
 */
enum osh_status osh_scoring_file_sink_init(struct osh_scoring_file_sink *fs,
                                           struct osh_scoring_workspace const *ws,
                                           struct osh_scoring_sink *out);

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_SINK_H */
