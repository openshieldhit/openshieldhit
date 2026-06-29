#ifndef OSH_SCORING_POSTPROCESS_H
#define OSH_SCORING_POSTPROCESS_H

#include "openshieldhit/status.h"
#include "scoring/runtime/osh_scoring_defs.h"
#include "scoring/runtime/osh_scoring_runtime.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Apply postprocessing to accumulated scoring page data, in place.
 *
 * @details
 * A separate cold-path phase between raw accumulation and save.  For DOSEGY it
 * rescales MeV/g -> Gy; for the LET/Qeff two-pass averages it finalises
 * data/data2; all other kinds are a no-op.  Implemented as the @p dst == @p src
 * case of @ref osh_scoring_postprocess_into, so its (destructive) behaviour is
 * identical to before.
 */
enum osh_status osh_scoring_postprocess(struct osh_scoring_runtime *rt);

/**
 * @brief Out-of-place postprocess: read @p src, write results into @p dst.
 *
 * @details
 * Computes the presentation form of every page of @p src into the matching page
 * of @p dst, **never mutating @p src**.  Only @c acc.data is ever written (DOSEGY
 * rescale, LET/Qeff average); for the LET/Qeff kinds the @c has_data2 / @c divide
 * flags are cleared on the @p dst page only.  Simple scorers are a no-op: when a
 * @p dst page aliases its @p src page's @c data (the shadow's non-transformed
 * pages, or the in-place @p dst == @p src wrapper) nothing is copied; when it owns
 * a distinct buffer the raw values are mirrored so the saved view is complete.
 *
 * @p dst and @p src must have the same @c npages, and matching pages the same
 * @c acc.len.  This is the non-destructive primitive a mid-run snapshot uses
 * (see @ref osh_scoring_shadow): clone only what postprocess writes, alias the
 * rest, and leave the live accumulators byte-identical.
 */
enum osh_status osh_scoring_postprocess_into(struct osh_scoring_runtime *dst, struct osh_scoring_runtime const *src);

/**
 * @brief True for score kinds whose postprocess writes @c acc.data.
 *
 * @details
 * DOSEGY (rescale) and the LET/Qeff averages (data/data2) are the only kinds
 * that mutate @c data; every other kind is a no-op.  A non-destructive snapshot
 * therefore needs a private @c data buffer only for these pages — the single
 * source of truth shared by the shadow and the memory estimate.
 */
int osh_scoring_postprocess_writes_data(enum osh_scoring_score_kind kind);

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_POSTPROCESS_H */
