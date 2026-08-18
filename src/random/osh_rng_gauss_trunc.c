#include "random/osh_rng_gauss_trunc.h"

#include <float.h>
#include <math.h>

#include "openshieldhit/const.h"
#include "random/osh_rng.h"

/* ---- Standard normal CDF and its complement -------------------------------- */

double osh_norm_cdf(double x) {
    return 0.5 * erfc(-x * OSH_M_SQRT1_2);
}

double osh_norm_cdf_upper(double x) {
    return 0.5 * erfc(x * OSH_M_SQRT1_2);
}

/* ---- Probit -------------------------------------------------------------- */

/*
 * Acklam's rational approximation to Phi^-1, ~1.15e-9 relative worst case.
 * Two regions: a rational in (p - 1/2)^2 near the centre, and a rational in
 * sqrt(-2 log p) in the tail, switching at p = 0.02425.
 *
 * The coefficients below are Acklam's; his own page has been offline for years,
 * so the reference for them (and for the derivation) is
 * L. M. Barros, "Acklam's Algorithm for the Inverse Normal CDF", 2017:
 * https://stackedboxes.org/2017/05/01/acklams-normal-quantile-function/
 *
 * The region test is a real branch rather than a blend of both rationals.
 * That is the right trade for scalar code: the beam sampler is called once
 * per primary from a loop that seeds a fresh RNG stream per history
 * (see fill_from_spots() in beam/runtime/), so there is no batch to
 * vectorise and nothing to gain from making the body branch-free.
 */
static double const ppf_a[6] = {-3.969683028665376e+01,
                                2.209460984245205e+02,
                                -2.759285104469687e+02,
                                1.383577518672690e+02,
                                -3.066479806614716e+01,
                                2.506628277459239e+00};

static double const ppf_b[5] = {-5.447609879822406e+01,
                                1.615858368580409e+02,
                                -1.556989798598866e+02,
                                6.680131188771972e+01,
                                -1.328068155288572e+01};

static double const ppf_c[6] = {-7.784894002430293e-03,
                                -3.223964580411365e-01,
                                -2.400758277161838e+00,
                                -2.549732539343734e+00,
                                4.374664141464968e+00,
                                2.938163982698783e+00};

static double const ppf_d[4] = {
    7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00};

/* Region boundary; below this, min(p, 1-p) uses the tail rational. */
#define OSH_PPF_TAIL 0.02425

double osh_norm_ppf(double p, double q) {
    double pm;
    double sgn;
    double y;
    double r;
    double s;
    double num;
    double den;

    /* Fold the upper tail onto the lower one: work with the smaller of the
     * two probabilities, which is the one that still has relative accuracy,
     * and undo the reflection with a sign. */
    if (p < q) {
        pm = p;
        sgn = 1.0;
    } else {
        pm = q;
        sgn = -1.0;
    }

    if (pm >= OSH_PPF_TAIL) {
        /* Central region.  0.5 * (p - q) equals p - 0.5 whenever p + q == 1,
         * and is better conditioned when the caller's p and q were computed
         * independently. */
        y = 0.5 * (p - q);
        r = y * y;
        num = (((((ppf_a[0] * r + ppf_a[1]) * r + ppf_a[2]) * r + ppf_a[3]) * r + ppf_a[4]) * r + ppf_a[5]) * y;
        den = ((((ppf_b[0] * r + ppf_b[1]) * r + ppf_b[2]) * r + ppf_b[3]) * r + ppf_b[4]) * r + 1.0;
        return num / den;
    }

    /* Tail region.  Clamp the argument away from zero: an underflowed pm
     * would give log(0) = -inf and then inf/inf = NaN in the rational.  The
     * clamp saturates the result near +/-37.5 sigma, which is where the
     * double-precision normal CDF runs out of range anyway. */
    if (pm < DBL_MIN) {
        pm = DBL_MIN;
    }
    s = sqrt(-2.0 * log(pm));
    num = ((((ppf_c[0] * s + ppf_c[1]) * s + ppf_c[2]) * s + ppf_c[3]) * s + ppf_c[4]) * s + ppf_c[5];
    den = (((ppf_d[0] * s + ppf_d[1]) * s + ppf_d[2]) * s + ppf_d[3]) * s + 1.0;
    return sgn * (num / den);
}

/* ---- Truncated normal ---------------------------------------------------- */

enum osh_status osh_gauss_trunc_prepare(struct osh_gauss_trunc *tg, double mu, double sigma, double lo, double hi) {
    double alpha;
    double beta;

    if (!tg || !(sigma > 0.0) || isnan(mu) || isnan(lo) || isnan(hi) || hi < lo) {
        return OSH_EINVAL;
    }

    alpha = (lo - mu) / sigma;
    beta = (hi - mu) / sigma;

    tg->mu = mu;
    tg->sigma = sigma;
    tg->lo = lo;
    tg->hi = hi;
    tg->p_lo = osh_norm_cdf(alpha);
    tg->q_lo = osh_norm_cdf_upper(alpha);

    /* Subtract the two small values, whichever tail those happen to be in,
     * so the difference keeps its significant digits. */
    if (alpha + beta <= 0.0) {
        tg->span = osh_norm_cdf(beta) - tg->p_lo;
    } else {
        tg->span = tg->q_lo - osh_norm_cdf_upper(beta);
    }

    tg->degenerate = !(tg->span > 0.0);
    return OSH_OK;
}

double osh_gauss_trunc_from_uniform(struct osh_gauss_trunc const *tg, double u) {
    double a;
    double b;
    double x;

    if (tg->degenerate) {
        /* The window carries no representable probability mass, so the
         * distribution over it is indistinguishable from uniform.  A window
         * this extreme can still have one infinite bound (e.g. (-inf, -50 sigma],
         * whose mass underflows); collapse to the finite edge in that case,
         * since that is where all the mass sits to within double precision.
         * Both bounds infinite gives span == 1 and never reaches here. */
        a = tg->lo;
        b = tg->hi;
        if (a == -HUGE_VAL) {
            a = b;
        }
        if (b == HUGE_VAL) {
            b = a;
        }
        return a + u * (b - a);
    }

    /* p and q are each accumulated on their own side rather than one being
     * derived from the other, so both keep relative accuracy for the probit. */
    x = tg->mu + tg->sigma * osh_norm_ppf(tg->p_lo + u * tg->span, tg->q_lo - u * tg->span);

    /* Guard only: the probit saturates near +/-37.5 sigma and the arithmetic
     * above can round a hair past a bound.  Unlike the clamp in a bounded
     * rejection sampler this creates no point mass at the cut -- it fires on
     * rounding, not on a routine failure to accept a draw. */
    if (x < tg->lo) {
        x = tg->lo;
    } else if (x > tg->hi) {
        x = tg->hi;
    }
    return x;
}

double osh_rng_gauss_trunc(struct osh_gauss_trunc const *tg, struct osh_rng *rng) {
    return osh_gauss_trunc_from_uniform(tg, osh_rng_double(rng));
}
