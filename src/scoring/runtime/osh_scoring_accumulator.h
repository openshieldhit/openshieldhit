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
 *   - @p data      primary accumulator (always allocated).
 *   - @p data2     secondary weight accumulator for two-pass averages
 *                  (dose/track-averaged LET and Qeff); NULL for simple scorers.
 *   - @p data_var  running sum-of-squares for @p data variance; NULL until the
 *                  variance (x²/σ) feature allocates it.
 *   - @p data2_var running sum-of-squares for @p data2 variance; NULL until then.
 *   - @p len       element count of every allocated array.
 */
struct osh_scoring_accumulator {
    double *data;
    double *data2;
    double *data_var;
    double *data2_var;
    size_t len;
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
 * @brief Free all arrays and zero the struct.
 *
 * Safe to call with @p acc NULL and on a partially-allocated struct.
 */
void osh_scoring_accumulator_free(struct osh_scoring_accumulator *acc);

/**
 * @brief Element-wise reduce: @p dst += @p src over every array.
 *
 * @details
 * The reduction at the heart of any parallel scoring scheme: each worker scores
 * its histories into a private accumulator, then the master set is formed by
 * merging every worker's set into one.  Because Monte Carlo histories are
 * independent and the deposits commute, merging in any order yields the same
 * totals (up to floating-point summation order).
 *
 * Adds @p src into @p dst for @c data, @c data2, @c data_var and @c data2_var.
 * An array is summed only when *both* sides have it allocated; an array present
 * on one side but NULL on the other is skipped (mismatched two-pass/variance
 * configuration is the caller's responsibility — accumulators compiled from the
 * same page descriptor always agree).
 *
 * @param[in,out] dst  Destination accumulator (must be non-NULL).
 * @param[in]     src  Source accumulator to add (must be non-NULL).
 * @returns OSH_OK on success, OSH_EINVAL if either pointer is NULL or the two
 *          accumulators have different @c len.
 */
enum osh_status osh_scoring_accumulator_merge(struct osh_scoring_accumulator *dst,
                                              struct osh_scoring_accumulator const *src);

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_ACCUMULATOR_H */
