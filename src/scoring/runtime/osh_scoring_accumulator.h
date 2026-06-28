#ifndef OSH_SCORING_ACCUMULATOR_H
#define OSH_SCORING_ACCUMULATOR_H

#include <stddef.h>

#include "openshieldhit/status.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Owning storage for one scored page's flat accumulator arrays.
 *
 * @details
 * Separated from the page *descriptor* (@ref osh_scoring_page_runtime, which
 * holds geometry indices, strides and differential-axis configuration) so that
 * the mutable, per-history-written storage has an identity of its own.  This is
 * the unit a future parallel worker owns privately: each worker scores its
 * assigned histories into a private accumulator set and then folds them into the
 * master set with @ref osh_scoring_accumulator_merge.  In the current
 * single-worker build every page simply owns one accumulator inline.
 *
 * Array layout matches the page descriptor: row-major spatial bins, with the
 * differential axes appended as outer strides (see osh_scoring_page_runtime).
 *
 *   - @p data      primary accumulator (always allocated): the running sum of
 *                  per-history deposits Σx in each bin (additive across batches).
 *   - @p data2     secondary weight accumulator for two-pass averages
 *                  (dose/track-averaged LET and Qeff); NULL for simple scorers.
 *   - @p data_var  Welford M2 (sum of squared deviations from the mean) for the
 *                  @p data estimator; NULL until the variance feature allocates
 *                  it.  See "Uncertainty representation" below.
 *   - @p data2_var Welford M2 for the @p data2 estimator; NULL until then.
 *   - @p len       element count of every allocated array.
 *   - @p weight    statistical weight (history count W) of the batch these
 *                  accumulators represent — the Schubert–Gertz merge weight.
 *   - @p nbatch    number of independent batches folded in so far (B); the
 *                  degrees of freedom for the error estimate are B − 1.
 *
 * @par Uncertainty representation (batch-means, Welford/Schubert–Gertz)
 * Monte-Carlo uncertainty here is a **batch-means** estimate: a *batch* is one
 * independent unit of work — a parallel worker's history range, a periodic
 * partial-dump interval, an internal sub-split of a serial run, or a separate
 * run — and each batch is treated as a single weighted observation of the
 * per-primary mean @c data/@c weight.  True per-history variance would need an
 * @c nbins × live-histories scratch buffer (prohibitive for the wavefront
 * design), so the spread *between batches* is the estimator instead.  This is
 * why @p weight / @p nbatch are per-accumulator scalars, not per-bin arrays:
 * every bin is exposed to the same set of histories (most depositing zero), so
 * the observation count is identical across bins.
 *
 * The variance state is combined with the numerically-stable parallel formula of
 * Schubert & Gertz (2018): merging batch A (weight @c wA, mean @c mA, M2 @c M2A)
 * with batch B yields
 *
 *   w  = wA + wB
 *   δ  = mB − mA                       (means derived as data/weight)
 *   M2 = M2A + M2B + δ² · wA·wB / w
 *
 * The δ²·wA·wB/w cross-term is exactly what a plain element-wise @c += cannot
 * express, and it handles **unequal-size batches by construction** — the normal
 * case under heterogeneous CPU cores and arbitrary dump boundaries.  @p data and
 * @p data2 (raw sums) and @p weight / @p nbatch stay additive; only the M2 arrays
 * use the cross-term.  Folding a brand-new batch (one observation, M2 = 0) is the
 * single-observation special case, equivalent to West's (1979) weighted update.
 *
 * A single batch (B = 1, e.g. a plain serial run with no sub-splitting) has zero
 * degrees of freedom and therefore no error estimate — at least two batches
 * (threads, dumps, sub-splits, or independent runs) are required.  At finalize
 * time the standard error of the per-primary mean in a bin is
 * @c sqrt(M2 / ((nbatch − 1) · weight)); this conversion lives in the (not-yet-
 * wired) variance feature, not here — this struct stores only the sufficient
 * statistics and the merge contract that combines them.
 */
struct osh_scoring_accumulator {
    double *data;
    double *data2;
    double *data_var;
    double *data2_var;
    size_t len;
    double weight;             /* Σ history weight (W) of the batches folded in; 0 when variance inactive. */
    unsigned long long nbatch; /* Number of batches folded in (B); error-estimate dof is B − 1. */
};

/**
 * @brief Single accumulation seam: add @p value into @p arr[@p idx].
 *
 * @details
 * Every scoring deposit funnels through this one inline so that the write
 * policy lives in exactly one place.  Today it is a plain `+=`, which the
 * compiler inlines to the same machine code as a direct array store (zero
 * overhead in the single-worker build).  A future parallel backend can change
 * the body here — to a relaxed atomic fetch-add, a per-worker private store, or
 * a lock-guarded update — without touching any of the call sites in the scoring
 * hot path.  @p arr is the raw array (e.g. @c page->acc.data or
 * @c page->acc.data2); callers compute @p idx including any differential-axis
 * offsets.
 */
static inline void osh_score_deposit(double *arr, size_t idx, double value) {
    arr[idx] += value;
}

/**
 * @brief Allocate zero-initialised accumulator arrays.
 *
 * @details
 * Allocates @p data always, and @p data2 when @p want_data2 is non-zero.  The
 * variance arrays are left NULL (the variance feature will allocate them when it
 * lands).  A @p len of 0 is rounded up to a single element so the pointers are
 * never NULL on success, matching the previous inline behaviour in
 * osh_scoring_compile().
 *
 * @param[out] acc        Accumulator to populate (overwritten; must not be NULL).
 * @param[in]  len        Element count for each array.
 * @param[in]  want_data2 Non-zero to also allocate the secondary @p data2 array.
 * @returns OSH_OK on success, OSH_EINVAL if @p acc is NULL, OSH_ENOMEM on
 *          allocation failure (no arrays leak; @p acc is left freed).
 */
enum osh_status osh_scoring_accumulator_alloc(struct osh_scoring_accumulator *acc, size_t len, int want_data2);

/**
 * @brief Reset every allocated array to zero without freeing.
 *
 * Intended for reusing a worker's accumulator set across runs.  NULL arrays are
 * skipped; safe to call on a zeroed struct.
 */
void osh_scoring_accumulator_zero(struct osh_scoring_accumulator *acc);

/**
 * @brief Multiply every element of the primary @p data array by @p factor in place.
 *
 * @details
 * The post-process scaling step (e.g. MeV/g -> Gy for absorbed dose) expressed as
 * an accumulator operation, so callers do not index the arrays themselves.  Only
 * @p data is touched; @p data2 and the variance arrays are left unchanged.  No-op
 * when @p acc or @p acc->data is NULL.
 */
void osh_scoring_accumulator_rescale(struct osh_scoring_accumulator *acc, double factor);

/**
 * @brief Finalise a two-pass weighted average in place: @c data[i] /= @c data2[i].
 *
 * @details
 * For track-/dose-averaged quantities (LET, Qeff), @p data holds the running
 * weighted sum and @p data2 the running weight.  Each bin becomes
 * @c data2[i] > @p eps ? data[i]/data2[i] : 0, so empty bins read as zero instead
 * of dividing by ~0.  @p data2 is left untouched (the caller decides whether to
 * keep or drop it).  Both @p data and @p data2 must be allocated.
 *
 * @param[in,out] acc  Accumulator to finalise (must be non-NULL with data + data2).
 * @param[in]     eps  Lower threshold below which a weight counts as zero.
 * @returns OSH_OK on success, OSH_EINVAL if @p acc, @p acc->data or @p acc->data2
 *          is NULL.
 */
enum osh_status osh_scoring_accumulator_finalize_average(struct osh_scoring_accumulator *acc, double eps);

/**
 * @brief Free all arrays and zero the struct.
 *
 * Safe to call with @p acc NULL and on a partially-allocated struct.
 */
void osh_scoring_accumulator_free(struct osh_scoring_accumulator *acc);

/**
 * @brief Reduce @p src into @p dst: combine two batches' statistics.
 *
 * @details
 * The reduction at the heart of any parallel scoring scheme: each worker scores
 * its histories into a private accumulator (one batch), then the master set is
 * formed by merging every worker's set into one.  Because Monte Carlo histories
 * are independent and the deposits commute, merging in any order yields the same
 * totals (up to floating-point summation order).
 *
 * The reduction is **representation-aware**, not a blanket element-wise add:
 *   - @c data, @c data2 (raw sums) and @c weight, @c nbatch are **additive** —
 *     @p dst += @p src.
 *   - @c data_var, @c data2_var hold Welford M2 and are combined with the
 *     Schubert–Gertz cross-term @c δ²·wA·wB/w (see @ref osh_scoring_accumulator),
 *     using the pre-merge sums and weights to derive the batch means.  A plain
 *     @c += over these arrays would silently drop the cross-term and corrupt the
 *     variance, so it is deliberately *not* used.
 *
 * The variance arrays are combined first (they need @p dst's pre-merge @c data /
 * @c weight), then the additive fields are folded.  Empty batches (@c weight 0)
 * are handled as the identity.  When the variance arrays are NULL — the current
 * default, since the variance feature is unwired — the merge reduces to the
 * familiar additive behaviour over @c data / @c data2.
 *
 * All optional arrays must have matching presence between @p dst and @p src: an
 * array allocated on one side but NULL on the other is rejected (OSH_EINVAL)
 * rather than silently dropped.  Accumulators compiled from the same page
 * descriptor always agree.
 *
 * @note MPI: the additive fields can ride @c MPI_SUM, but the M2 arrays cannot —
 *       a rank-level reduction needs a custom @c MPI_Op_create that applies this
 *       same combine to the whole accumulator (data + weight + M2 together).
 *
 * @param[in,out] dst  Destination accumulator (must be non-NULL).
 * @param[in]     src  Source accumulator to fold in (must be non-NULL).
 * @returns OSH_OK on success, OSH_EINVAL if either pointer is NULL, the two
 *          accumulators have different @c len, or their optional arrays
 *          (@c data2 / @c data_var / @c data2_var) disagree on presence.
 */
enum osh_status osh_scoring_accumulator_merge(struct osh_scoring_accumulator *dst,
                                              struct osh_scoring_accumulator const *src);

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_ACCUMULATOR_H */
