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
 * A separate cold-path phase between raw accumulation and save.  For DOSE,
 * FLUENCE, and NKERMA it applies the per-bin volume normalisation; for DOSEGY it
 * applies the same volume normalisation and then rescales MeV/g -> Gy; for the
 * LET/Qeff two-pass averages it finalises data/data2.  Score kinds without a
 * postprocess handler are already in final form unless they are differential
 * pages, where postprocess divides additive bin contents by the differential-axis
 * bin width(s) to report differential quantities.  Implemented as the
 * @p dst == @p src case of @ref osh_scoring_postprocess_into, so its
 * (destructive) behaviour is identical to before.
 *
 * @note Single-shot: the transform is destructive and not idempotent, so a
 * second in-place call would double-divide by volume / re-ratio data÷data2.
 * The runtime records that it has run and any further in-place call returns
 * @c OSH_ESTATE.  The guard lives in @ref osh_scoring_postprocess_into (the
 * @p dst == @p src case), so calling that directly with @p dst == @p src is
 * protected identically to this wrapper.  For repeatable (dump / checkpoint /
 * pre-merge) finalisation use @ref osh_scoring_postprocess_into with a
 * **distinct** @p dst — it never mutates @p src and does not trip this guard.
 *
 * @note BDO merge tools must follow the per-page @c OSHBDO_PAG_NORMALIZE tag
 * written from @c page->postproc (assigned at scoring compile time).  That tag,
 * not the unit string, tells the merge whether independent files are summed,
 * nstat-normalised, nstat-weighted averages, or appended.  Differential NORM
 * pages are already bin-width-normalised here, so they still merge by nstat
 * weighting, never by summing reported differential values directly.
 */
enum osh_status osh_scoring_postprocess(struct osh_scoring_runtime *rt);

/**
 * @brief Out-of-place postprocess: read @p src, write results into @p dst.
 *
 * @details
 * Computes the presentation form of every page of @p src into the matching page
 * of @p dst, **never mutating @p src**.  Only @c acc.data is ever written
 * (volume normalisation, DOSEGY rescale, differential bin-width normalisation,
 * LET/Qeff average); for the LET/Qeff kinds the @c has_data2 / @c divide flags
 * are cleared on the @p dst page only.  Score kinds without a postprocess handler
 * are mirrored as raw data, or divided by differential bin width(s) when they are
 * differential additive pages.
 *
 * @p dst and @p src must have the same @c npages, and matching pages the same
 * @c acc.len.  With a **distinct** @p dst this is the non-destructive primitive a
 * mid-run snapshot uses (see @ref osh_scoring_shadow): clone only what postprocess
 * writes, alias the rest, and leave the live accumulators byte-identical.
 *
 * @note The in-place case (@p dst == @p src) is destructive and single-shot: it
 * returns @c OSH_ESTATE if the runtime was already postprocessed, and sets that
 * flag on success.  This is the guard @ref osh_scoring_postprocess relies on, and
 * it protects a direct in-place call here too.  An out-of-place call never touches
 * the flag and stays repeatable.
 */
enum osh_status osh_scoring_postprocess_into(struct osh_scoring_runtime *dst, struct osh_scoring_runtime const *src);

/**
 * @brief Finalise per-bin Monte-Carlo standard error into the M2 arrays, in place.
 *
 * @details
 * Converts each variance-tracking page's batch-means sufficient statistics into a
 * per-bin **relative standard error** of its reported quantity, stored back into
 * @c acc.data_var (spent afterwards).  Storing the *relative* error makes it
 * normalisation-invariant: the save layer emits the absolute error column as
 * @c |value| * data_var, so whatever per-primary / per-volume / unit scaling the
 * writer applies to the value applies to its error automatically.
 *
 * With @c B = @c nbatch batches, @c W = @c weight (Σ history count) and Welford
 * @c M2 in @c data_var, the standard error of the per-primary mean is
 * @c sqrt(M2 / ((B-1)·W)); dividing by the mean @c data/W gives the relative form
 * @c sqrt(M2·W/(B-1)) / |data|.  For the two-pass AVER quantities (DLET/TLET/Qeff)
 * the reported value is the ratio @c data/data2, whose relative error is combined
 * from the numerator (@c data_var) and denominator (@c data2_var) in quadrature
 * (@c rel² = rel_num² + rel_den²).  This ignores the num/den covariance, which for
 * the strongly-correlated LET numerator and denominator makes the estimate
 * *conservative* (an upper bound), never an under-estimate.
 *
 * @b Ordering: this reads the raw sums (@c data, @c data2) and their M2 arrays, so
 * it must run **before** @ref osh_scoring_postprocess rescales @c data or collapses
 * the LET ratio.  @c B < 2 (zero degrees of freedom), a zero mean, or a zero
 * weight yield a 0 error for that bin.  Pages without M2 arrays (variance off, or
 * the dump shadow) are skipped, so this is a no-op unless a page enabled variance
 * tracking via a "Variance On" Settings block.
 *
 * @param[in,out] rt  Scoring runtime whose variance pages are finalised in place.
 * @returns OSH_OK on success, OSH_EINVAL if @p rt is NULL.
 */
enum osh_status osh_scoring_finalize_errors(struct osh_scoring_runtime *rt);

/**
 * @brief True for score kinds whose postprocess writes @c acc.data.
 *
 * @details
 * DOSE, FLUENCE, and NKERMA need volume normalisation; DOSEGY also rescales
 * MeV/g -> Gy; LET/Qeff averages finalise data/data2.  A non-destructive
 * snapshot therefore needs a private @c data buffer only for these pages — the
 * single source of truth shared by the shadow and the memory estimate.
 */
int osh_scoring_postprocess_writes_data(enum osh_scoring_score_kind kind);

/**
 * @brief True when a specific page writes @c acc.data during postprocess.
 *
 * @details
 * Extends @ref osh_scoring_postprocess_writes_data with page-level differential
 * metadata, so additive pages such as Energy get scratch in non-destructive
 * snapshots when they need differential bin-width normalisation.
 */
int osh_scoring_postprocess_page_writes_data(struct osh_scoring_page_runtime const *page);

/* Per-estimator postprocess handlers, registered in the estimator table
 * (osh_scoring_estimator.c).  Each finalises one page's accumulator:
 *   postprocess_volume — DOSE/FLUENCE/NKERMA: multiply each bin by geo->bin_vol_inv.
 *   postprocess_dosegy — as volume, then MeV/g -> Gy.
 *   postprocess_ratio  — DLET/TLET/DQEFF/TQEFF: finalise the data/data2 average. */
enum osh_status postprocess_volume(struct osh_scoring_page_runtime *dst,
                                   struct osh_scoring_page_runtime const *src,
                                   struct osh_scoring_geometry_runtime const *geo);
enum osh_status postprocess_dosegy(struct osh_scoring_page_runtime *dst,
                                   struct osh_scoring_page_runtime const *src,
                                   struct osh_scoring_geometry_runtime const *geo);
enum osh_status postprocess_ratio(struct osh_scoring_page_runtime *dst,
                                  struct osh_scoring_page_runtime const *src,
                                  struct osh_scoring_geometry_runtime const *geo);

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_POSTPROCESS_H */
