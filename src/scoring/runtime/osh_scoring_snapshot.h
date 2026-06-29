#ifndef OSH_SCORING_SNAPSHOT_H
#define OSH_SCORING_SNAPSHOT_H

#include <stddef.h>

#include "openshieldhit/status.h"
#include "scoring/runtime/osh_scoring_shadow.h"
#include "scoring/save/osh_scoring_sink.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Non-destructive snapshot + save — primitive P2 (issue #170 / #191).
 *
 * @details
 * Refreshes @p shadow (an out-of-place postprocess of its bound live runtime into
 * private scratch — see @ref osh_scoring_shadow) and hands the resulting
 * presentation view to @p sink, normalising by @p completed_nstat.  The live
 * accumulators are never mutated, so the run can keep accumulating after a dump.
 *
 * @param[in]     sink            Destination (G1); must have a non-NULL @c save.
 * @param[in,out] shadow          Caller-owned shadow bound to the live runtime;
 *                                its scratch is allocated once and reused.
 * @param[in]     completed_nstat Primary count to normalise by; must be > 0 for
 *                                the file sink.
 * @param[in]     want            Output selector (G2): indices to save, NULL = all.
 * @param[in]     n_want          Number of entries in @p want.
 * @returns OSH_OK, OSH_EINVAL on a NULL argument, or the sink's / refresh's error.
 */
enum osh_status osh_scoring_snapshot_save(struct osh_scoring_sink const *sink,
                                          struct osh_scoring_shadow *shadow,
                                          unsigned long long completed_nstat,
                                          size_t const *want,
                                          size_t n_want);

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_SNAPSHOT_H */
