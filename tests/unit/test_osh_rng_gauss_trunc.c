#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "openshieldhit/status.h"
#include "random/osh_rng.h"
#include "random/osh_rng_gauss_trunc.h"

#define ASSERT_TRUE(cond)                                                                                              \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            fprintf(stderr, "ASSERT FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__);                                 \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

#define NSAMP 200000

/* Standard normal density. */
static double npdf(double z) {
    return exp(-0.5 * z * z) / 2.5066282746310002;
}

/*
 * Closed-form mean and standard deviation of the truncated normal, used as the
 * reference the sampler is measured against.  Terms of the form z*phi(z) are
 * zero at an infinite cut.
 */
static void trunc_moments(double mu, double sg, double lo, double hi, double *mean_out, double *sd_out) {
    double a;
    double b;
    double z;
    double pa;
    double pb;
    double apa;
    double bpb;
    double lam;

    a = (lo - mu) / sg;
    b = (hi - mu) / sg;
    z = osh_norm_cdf(b) - osh_norm_cdf(a);
    pa = isinf(a) ? 0.0 : npdf(a);
    pb = isinf(b) ? 0.0 : npdf(b);
    apa = isinf(a) ? 0.0 : a * pa;
    bpb = isinf(b) ? 0.0 : b * pb;
    lam = (pa - pb) / z;
    *mean_out = mu + sg * lam;
    *sd_out = sg * sqrt(1.0 + (apa - bpb) / z - lam * lam);
}

static int cmp_double(void const *lhs, void const *rhs) {
    double x = *(double const *) lhs;
    double y = *(double const *) rhs;

    return (x > y) - (x < y);
}

/*
 * Probit accuracy, checked in the direction where the reference is exact:
 * pick z, form p and q from erfc, and require osh_norm_ppf() to recover z to
 * Acklam's stated ~1.15e-9 relative.  Covers both tails out to the point where
 * the double-precision normal CDF underflows.
 */
static void test_probit_roundtrip_accuracy(void) {
    double const zs[] = {-37.0,
                         -30.0,
                         -12.0,
                         -8.0,
                         -5.0,
                         -3.0,
                         -1.96,
                         -1.0,
                         -0.25,
                         0.0,
                         0.25,
                         1.0,
                         1.96,
                         3.0,
                         5.0,
                         8.0,
                         12.0,
                         30.0,
                         37.0};
    double worst;
    double z;
    double zr;
    double err;
    size_t i;

    worst = 0.0;
    for (i = 0; i < sizeof zs / sizeof zs[0]; i++) {
        z = zs[i];
        zr = osh_norm_ppf(osh_norm_cdf(z), osh_norm_cdf_upper(z));
        err = fabs(zr - z) / (fabs(z) > 1.0 ? fabs(z) : 1.0);
        if (err > worst) {
            worst = err;
        }
    }
    ASSERT_TRUE(worst < 2.0e-9);

    /* Median is exact, and the tail clamp keeps an underflowed argument finite
     * instead of producing inf/inf = NaN. */
    ASSERT_TRUE(osh_norm_ppf(0.5, 0.5) == 0.0);
    ASSERT_TRUE(osh_norm_ppf(0.0, 1.0) < -37.0 && osh_norm_ppf(0.0, 1.0) > -38.0);
    ASSERT_TRUE(osh_norm_ppf(1.0, 0.0) > 37.0 && osh_norm_ppf(1.0, 0.0) < 38.0);
}

/*
 * Sampled mean and standard deviation against the closed forms, plus a
 * hard out-of-bounds count, over windows from mildly truncating to 8 sigma off
 * the mean.  The 8 sigma case is the one a bounded-retry rejection sampler
 * cannot do at all: its acceptance probability is 6.2e-16, so every draw would
 * exhaust the retries and clamp to the window edge.
 */
static void test_moments_match_closed_form(void) {
    struct {
        double mu;
        double sg;
        double lo;
        double hi;
    } const cases[] = {
        {60.0, 5.0, -HUGE_VAL, 61.0}, /* upper cut just above the mean */
        {60.0, 5.0, 59.0, 61.0},      /* narrow two-sided band around the mean */
        {60.0, 5.0, -HUGE_VAL, 45.0}, /* 3 sigma below: 0.13% acceptance */
        {60.0, 5.0, -HUGE_VAL, 20.0}, /* 8 sigma below: 6.2e-16 acceptance */
        {60.0, 5.0, 61.0, 62.0},      /* band entirely in the upper tail */
    };
    struct osh_gauss_trunc tg;
    struct osh_rng rng;
    double s1;
    double s2;
    double x;
    double mean;
    double sd;
    double exp_mean;
    double exp_sd;
    double se;
    size_t oob;
    size_t k;
    int i;
    enum osh_status rc;

    for (k = 0; k < sizeof cases / sizeof cases[0]; k++) {
        rc = osh_gauss_trunc_prepare(&tg, cases[k].mu, cases[k].sg, cases[k].lo, cases[k].hi);
        ASSERT_TRUE(rc == OSH_OK);
        ASSERT_TRUE(!tg.degenerate);

        osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, 20260818u, (uint64_t) k);
        s1 = 0.0;
        s2 = 0.0;
        oob = 0u;
        for (i = 0; i < NSAMP; i++) {
            x = osh_rng_gauss_trunc(&tg, &rng);
            s1 += x;
            s2 += x * x;
            if (x < cases[k].lo || x > cases[k].hi) {
                oob++;
            }
        }
        ASSERT_TRUE(oob == 0u);

        mean = s1 / (double) NSAMP;
        sd = sqrt(s2 / (double) NSAMP - mean * mean);
        trunc_moments(cases[k].mu, cases[k].sg, cases[k].lo, cases[k].hi, &exp_mean, &exp_sd);

        /* Five standard errors on the mean; the run is deterministic (fixed
         * seed), so this is a fixed margin, not a flakiness budget. */
        se = exp_sd / sqrt((double) NSAMP);
        ASSERT_TRUE(fabs(mean - exp_mean) < 5.0 * se);
        ASSERT_TRUE(fabs(sd - exp_sd) < 0.01 * exp_sd);
    }
}

/* Kolmogorov-Smirnov against the exact truncated CDF. */
static void test_ks_against_truncated_cdf(void) {
    struct osh_gauss_trunc tg;
    struct osh_rng rng;
    double *x;
    double cdf;
    double d;
    double ks;
    int i;
    enum osh_status rc;

    rc = osh_gauss_trunc_prepare(&tg, 60.0, 5.0, -HUGE_VAL, 61.0);
    ASSERT_TRUE(rc == OSH_OK);

    x = (double *) malloc((size_t) NSAMP * sizeof(double));
    ASSERT_TRUE(x != NULL);

    osh_rng_init(&rng, OSH_RNG_TYPE_XOSHIRO256SS, 4242u, 7u);
    for (i = 0; i < NSAMP; i++) {
        x[i] = osh_rng_gauss_trunc(&tg, &rng);
    }
    qsort(x, (size_t) NSAMP, sizeof(double), cmp_double);

    ks = 0.0;
    for (i = 0; i < NSAMP; i++) {
        cdf = (osh_norm_cdf((x[i] - 60.0) / 5.0) - tg.p_lo) / tg.span;
        d = fabs(cdf - (double) i / (double) NSAMP);
        if (d > ks) {
            ks = d;
        }
        d = fabs(cdf - (double) (i + 1) / (double) NSAMP);
        if (d > ks) {
            ks = d;
        }
    }
    free(x);

    /* 95% critical value for the two-sided one-sample KS statistic. */
    ASSERT_TRUE(ks * sqrt((double) NSAMP) < 1.36);
}

/* The quantile transform is monotone and hits both window edges exactly. */
static void test_uniform_endpoints_and_monotonicity(void) {
    struct osh_gauss_trunc tg;
    double prev;
    double x;
    double u;
    int i;
    enum osh_status rc;

    rc = osh_gauss_trunc_prepare(&tg, 60.0, 5.0, 58.0, 62.0);
    ASSERT_TRUE(rc == OSH_OK);

    ASSERT_TRUE(fabs(osh_gauss_trunc_from_uniform(&tg, 0.0) - 58.0) < 1e-9);
    ASSERT_TRUE(fabs(osh_gauss_trunc_from_uniform(&tg, 1.0) - 62.0) < 1e-9);
    /* u = 0.5 lands on the median of the truncated law, which is t0 for a
     * window centred on t0. */
    ASSERT_TRUE(fabs(osh_gauss_trunc_from_uniform(&tg, 0.5) - 60.0) < 1e-8);

    prev = -HUGE_VAL;
    for (i = 0; i <= 1000; i++) {
        u = (double) i / 1000.0;
        x = osh_gauss_trunc_from_uniform(&tg, u);
        ASSERT_TRUE(x >= prev);
        ASSERT_TRUE(x >= 58.0 && x <= 62.0);
        prev = x;
    }
}

/*
 * Exactly one deviate per draw, always.  This is what makes a run reproducible
 * across truncation settings: a rejection loop consumed a data-dependent number
 * of deviates, so changing the window shifted every subsequent draw in the
 * history -- position, direction, everything downstream.
 */
static void test_fixed_stream_consumption(void) {
    struct osh_gauss_trunc narrow;
    struct osh_gauss_trunc wide;
    struct osh_rng rng_a;
    struct osh_rng rng_b;
    struct osh_rng rng_ref;
    int i;
    enum osh_status rc;

    /* A window 3 sigma off the mean: acceptance 1.3e-3, i.e. ~740 draws per
     * accepted sample under rejection. */
    rc = osh_gauss_trunc_prepare(&narrow, 60.0, 5.0, 44.0, 45.0);
    ASSERT_TRUE(rc == OSH_OK);
    rc = osh_gauss_trunc_prepare(&wide, 60.0, 5.0, 50.0, 70.0);
    ASSERT_TRUE(rc == OSH_OK);

    osh_rng_init(&rng_a, OSH_RNG_TYPE_PCG32, 99u, 1u);
    osh_rng_init(&rng_b, OSH_RNG_TYPE_PCG32, 99u, 1u);
    osh_rng_init(&rng_ref, OSH_RNG_TYPE_PCG32, 99u, 1u);

    for (i = 0; i < 64; i++) {
        (void) osh_rng_gauss_trunc(&narrow, &rng_a);
        (void) osh_rng_gauss_trunc(&wide, &rng_b);
        (void) osh_rng_double(&rng_ref);
    }

    /* All three streams sit at the same position, so whatever is drawn next is
     * identical regardless of which window was sampled. */
    ASSERT_TRUE(osh_rng_double(&rng_a) == osh_rng_double(&rng_ref));
    osh_rng_init(&rng_ref, OSH_RNG_TYPE_PCG32, 99u, 1u);
    for (i = 0; i < 64; i++) {
        (void) osh_rng_double(&rng_ref);
    }
    ASSERT_TRUE(osh_rng_double(&rng_b) == osh_rng_double(&rng_ref));
}

/*
 * Windows whose probability mass underflows to zero: the sampler must stay
 * inside the window and stay finite, including when one bound is infinite.
 */
static void test_degenerate_windows(void) {
    struct osh_gauss_trunc tg;
    double x;
    int i;
    enum osh_status rc;

    /* Two-sided, 40 sigma up: Phi complement underflows, so the span does too. */
    rc = osh_gauss_trunc_prepare(&tg, 0.0, 1.0, 40.0, 41.0);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(tg.degenerate);
    for (i = 0; i <= 10; i++) {
        x = osh_gauss_trunc_from_uniform(&tg, (double) i / 10.0);
        ASSERT_TRUE(x >= 40.0 && x <= 41.0);
    }
    /* The fallback is uniform on the window, so it does spread. */
    ASSERT_TRUE(osh_gauss_trunc_from_uniform(&tg, 1.0) - osh_gauss_trunc_from_uniform(&tg, 0.0) > 0.9);

    /* One-sided and unrepresentable: everything collapses onto the finite edge
     * rather than producing NaN from (-inf) + u * (hi - -inf). */
    rc = osh_gauss_trunc_prepare(&tg, 60.0, 1.0, -HUGE_VAL, 15.0);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(tg.degenerate);
    for (i = 0; i <= 10; i++) {
        x = osh_gauss_trunc_from_uniform(&tg, (double) i / 10.0);
        ASSERT_TRUE(x == 15.0);
    }
}

/* A degenerate window is a zero-width one, and it is not an error. */
static void test_zero_width_window(void) {
    struct osh_gauss_trunc tg;
    enum osh_status rc;

    rc = osh_gauss_trunc_prepare(&tg, 60.0, 5.0, 61.0, 61.0);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(tg.degenerate);
    ASSERT_TRUE(osh_gauss_trunc_from_uniform(&tg, 0.0) == 61.0);
    ASSERT_TRUE(osh_gauss_trunc_from_uniform(&tg, 1.0) == 61.0);
}

/* An untruncated "window" must reproduce the plain normal quantiles. */
static void test_unbounded_window_is_plain_normal(void) {
    struct osh_gauss_trunc tg;
    enum osh_status rc;

    rc = osh_gauss_trunc_prepare(&tg, 60.0, 5.0, -HUGE_VAL, HUGE_VAL);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(!tg.degenerate);
    ASSERT_TRUE(fabs(tg.span - 1.0) < 1e-15);
    ASSERT_TRUE(fabs(osh_gauss_trunc_from_uniform(&tg, 0.5) - 60.0) < 1e-9);
    /* 0.9750021048517795 is Phi(1.96). */
    ASSERT_TRUE(fabs(osh_gauss_trunc_from_uniform(&tg, 0.9750021048517795) - (60.0 + 5.0 * 1.96)) < 1e-7);
}

static void test_prepare_rejects_bad_parameters(void) {
    struct osh_gauss_trunc tg;

    ASSERT_TRUE(osh_gauss_trunc_prepare(NULL, 60.0, 5.0, 58.0, 62.0) == OSH_EINVAL);
    ASSERT_TRUE(osh_gauss_trunc_prepare(&tg, 60.0, 0.0, 58.0, 62.0) == OSH_EINVAL);
    ASSERT_TRUE(osh_gauss_trunc_prepare(&tg, 60.0, -1.0, 58.0, 62.0) == OSH_EINVAL);
    ASSERT_TRUE(osh_gauss_trunc_prepare(&tg, 60.0, 5.0, 62.0, 58.0) == OSH_EINVAL);
    ASSERT_TRUE(osh_gauss_trunc_prepare(&tg, 60.0, 5.0, 58.0, (double) NAN) == OSH_EINVAL);
}

int main(void) {
    test_probit_roundtrip_accuracy();
    test_moments_match_closed_form();
    test_ks_against_truncated_cdf();
    test_uniform_endpoints_and_monotonicity();
    test_fixed_stream_consumption();
    test_degenerate_windows();
    test_zero_width_window();
    test_unbounded_window_is_plain_normal();
    test_prepare_rejects_bad_parameters();

    return 0;
}
