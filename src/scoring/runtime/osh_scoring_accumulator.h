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

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_ACCUMULATOR_H */
