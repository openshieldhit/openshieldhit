#ifndef OSH_RNG_GAUSS_TRUNC_H
#define OSH_RNG_GAUSS_TRUNC_H

/**
 * @file osh_rng_gauss_trunc.h
 * @brief Truncated normal distribution by inverse-CDF sampling.
 *
 * @details
 * Draws from N(mu, sigma^2) restricted to [lo, hi] using **exactly one**
 * uniform deviate per sample, via
 *
 *     x = mu + sigma * Phi^-1( Phi(alpha) + u * (Phi(beta) - Phi(alpha)) )
 *
 * with the standardised cuts alpha = (lo - mu)/sigma, beta = (hi - mu)/sigma.
 * The two interval constants Phi(alpha) and the span Phi(beta) - Phi(alpha)
 * depend only on the truncation, so they are computed once by
 * osh_gauss_trunc_prepare() and reused for every draw.
 *
 * Nothing loops: no draw is ever discarded, so the cost is constant and the
 * result is exact wherever the window sits.  Discarding out-of-window draws
 * instead would cost 1/Z samples, where Z is the acceptance mass
 * (@ref osh_gauss_trunc::span) -- already ~740 per accepted sample for a window
 * 3 sigma off the mean -- and any bounded retry count would run out and have to
 * fall back to something that is no longer the requested distribution.  The
 * fixed deviate count also keeps a history's stream position independent of the
 * truncation.
 *
 * Accuracy is set by Acklam's rational approximation to the probit,
 * ~1.15e-9 relative worst case.  No Newton/Halley polish step is applied: it
 * would pull erfc() and exp() into the per-particle path to refine a value
 * already far below any physically meaningful energy resolution.  The reference
 * for Acklam's coefficients is given at osh_norm_ppf() in the .c file.
 */

#include "openshieldhit/status.h"

struct osh_rng;

/**
 * @struct osh_gauss_trunc
 *
 * @brief Interval constants for one (mu, sigma, lo, hi) truncated normal.
 *
 * @details
 * Fill with osh_gauss_trunc_prepare() once per distinct parameter set, then
 * pass to osh_rng_gauss_trunc() / osh_gauss_trunc_from_uniform() per draw.
 * Plain value type: copy it freely, no ownership, no allocation.
 *
 * @p p_lo and @p q_lo are maintained as two independently computed values
 * rather than one and its complement, so each keeps full *relative* accuracy
 * on its own side of the distribution.  Carrying both is what lets the probit
 * stay accurate as the argument approaches 1, and removes the "mirror the
 * interval into the lower half" trick that would otherwise be needed.
 */
struct osh_gauss_trunc {
    double mu;      /**< Untruncated mean. */
    double sigma;   /**< Untruncated standard deviation, > 0. */
    double lo;      /**< Lower truncation bound; -HUGE_VAL for none. */
    double hi;      /**< Upper truncation bound; +HUGE_VAL for none. */
    double p_lo;    /**< Phi(alpha), the probability mass below @p lo. */
    double q_lo;    /**< 1 - Phi(alpha), computed directly from erfc, not as 1 - p_lo. */
    double span;    /**< Phi(beta) - Phi(alpha): the acceptance mass Z in (0, 1]. */
    int degenerate; /**< Non-zero when @p span underflowed and the uniform fallback is used. */
};

/**
 * @brief Standard normal CDF Phi(x).
 *
 * @details
 * Evaluated as erfc(-x/sqrt(2))/2, which keeps full relative accuracy for
 * x << 0 where Phi(x) is tiny.
 *
 * @param[in] x Standardised deviate.
 *
 * @returns Phi(x) in [0, 1].
 */
double osh_norm_cdf(double x);

/**
 * @brief Standard normal complementary CDF 1 - Phi(x).
 *
 * @details
 * Evaluated as erfc(x/sqrt(2))/2 — never as 1 - osh_norm_cdf(x), which would
 * lose every significant digit for x >> 0.
 *
 * @param[in] x Standardised deviate.
 *
 * @returns 1 - Phi(x) in [0, 1].
 */
double osh_norm_cdf_upper(double x);

/**
 * @brief Standard normal quantile (probit) Phi^-1(p), given p and 1 - p.
 *
 * @details
 * Takes the complement @p q as a separate argument instead of deriving it,
 * because the caller can usually compute both with better relative accuracy
 * than a subtraction allows (see @ref osh_gauss_trunc).  The two tails are
 * handled by one code path: the upper tail is the lower-tail formula applied
 * to min(p, q) with the sign flipped.
 *
 * @param[in] p Probability in [0, 1].
 * @param[in] q Complement 1 - p, computed independently where possible.
 *
 * @returns z such that Phi(z) = p, accurate to ~1.15e-9 relative.  Saturates
 *          near +/-37.5 instead of returning an infinity when min(p, q)
 *          underflows to zero.
 */
double osh_norm_ppf(double p, double q);

/**
 * @brief Precompute the interval constants for a truncated normal.
 *
 * @details
 * Costs up to four erfc() evaluations, so call it once per distinct
 * (mu, sigma, lo, hi) — typically once per beam spot at setup — and never
 * per particle.  Either bound may be infinite for a one-sided truncation.
 *
 * The span is formed by subtracting the two *smaller* of the four CDF values,
 * choosing the side by the sign of alpha + beta.  Subtracting the two large
 * values instead would cancel away the significant digits whenever both cuts
 * sit deep in the same tail.
 *
 * If the span underflows to zero (a window so narrow, or so far out in a
 * tail, that its probability mass is not representable), @p tg is marked
 * @ref osh_gauss_trunc::degenerate and sampling falls back to a uniform draw
 * on [lo, hi] — the correct limit for a vanishingly narrow window.  This is
 * still reported as OSH_OK; inspect the flag if the caller cares.
 *
 * @param[out] tg    State to fill; untouched on error.
 * @param[in]  mu    Untruncated mean.
 * @param[in]  sigma Untruncated standard deviation; must be > 0.
 * @param[in]  lo    Lower bound, or -HUGE_VAL for none.
 * @param[in]  hi    Upper bound, or +HUGE_VAL for none; must be >= @p lo.
 *
 * @returns OSH_OK on success, OSH_EINVAL if @p tg is NULL, @p sigma <= 0,
 *          either bound is NaN, or @p hi < @p lo.
 */
enum osh_status osh_gauss_trunc_prepare(struct osh_gauss_trunc *tg, double mu, double sigma, double lo, double hi);

/**
 * @brief Map one uniform deviate to a truncated normal draw.
 *
 * @details
 * Split out from osh_rng_gauss_trunc() so callers that already hold a uniform
 * (quasi-random or stratified sequences, replay from a recorded stream, unit
 * tests probing specific quantiles) can use the same transform.
 *
 * @param[in] tg Prepared state from osh_gauss_trunc_prepare().
 * @param[in] u  Uniform deviate in [0, 1]; both endpoints are safe and map to
 *               @p tg->lo and @p tg->hi respectively.
 *
 * @returns A value in [tg->lo, tg->hi].
 */
double osh_gauss_trunc_from_uniform(struct osh_gauss_trunc const *tg, double u);

/**
 * @brief Draw one truncated normal variate, consuming exactly one uniform.
 *
 * @param[in]     tg  Prepared state from osh_gauss_trunc_prepare().
 * @param[in,out] rng Random-number generator; advanced by one uniform draw.
 *
 * @returns A value in [tg->lo, tg->hi].
 */
double osh_rng_gauss_trunc(struct osh_gauss_trunc const *tg, struct osh_rng *rng);

#endif /* OSH_RNG_GAUSS_TRUNC_H */
