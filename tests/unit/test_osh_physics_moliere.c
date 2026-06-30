/* Unit tests for the Highland/Bethe-Molière multiple-scattering physics.
 *
 * The Molière pieces are validated against analytic constraints derived from
 * the published theory (Bethe, Phys. Rev. 89, 1256 (1953)) — f⁽⁰⁾ closed form,
 * the normalisation moments of f⁽¹⁾/f⁽²⁾, and the B equation — rather than
 * against any copied tabulated numbers. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "physics/atomic/osh_physics_scat.h"
#include "physics/atomic/osh_physics_scat_highland.h"
#include "physics/atomic/osh_physics_scat_moliere.h"
#include "random/osh_rng.h"
#include "test_assert.h"

/* ---- Mode 1 (Highland) path-scale additivity ----------------------------- */

/* Splitting a step into N substeps and summing the per-substep Gaussian
 * variances must reproduce the single full-step variance, because the log
 * correction is evaluated at the (shared) macroscopic path scale. */
static void test_theta0_path_scale_additivity(void) {
    double t; /* MeV */
    double m; /* MeV */
    double z;
    double x0; /* g/cm^2 (water) */
    double l;  /* total areal thickness g/cm^2 */
    double full;
    double sum;
    double th;
    int n;
    int i;

    t = 150.0;
    m = 938.272;
    z = 1.0;
    x0 = 36.08;
    l = 5.0;
    n = 50;
    full = osh_physics_highland_theta0(t, m, z, l, l, x0);
    sum = 0.0;
    for (i = 0; i < n; ++i) {
        th = osh_physics_highland_theta0(t, m, z, l / (double) n, l, x0);
        sum += th * th;
    }
    /* sqrt(sum of substep variances) == full-step theta0 */
    ASSERT_TRUE(full > 0.0);
    ASSERT_TRUE(fabs(sqrt(sum) - full) < 1e-9 * full);
}

/* The legacy single-step behaviour (path_scale <= 0 falls back to thickness). */
static void test_theta0_fallback_to_thickness(void) {
    double t;
    double m;
    double a;
    double b;

    t = 150.0;
    m = 938.272;
    a = osh_physics_highland_theta0(t, m, 1.0, 0.5, 0.5, 36.08);
    b = osh_physics_highland_theta0(t, m, 1.0, 0.5, -1.0, 36.08);
    ASSERT_TRUE(fabs(a - b) < 1e-12);
}

/* ---- Molière B equation -------------------------------------------------- */

static void test_solve_b_equation(void) {
    double omega[] = {3.5, 1.0e2, 1.0e4, 1.0e6};
    double b;
    int ok;
    size_t i;

    for (i = 0; i < sizeof(omega) / sizeof(omega[0]); ++i) {
        ok = osh_physics_moliere_solve_b(omega[i], &b);
        ASSERT_TRUE(ok == 1);
        ASSERT_TRUE(b > 1.0);
        /* B - ln B == ln(omega) */
        ASSERT_TRUE(fabs((b - log(b)) - log(omega[i])) < 1e-6);
    }
}

static void test_solve_b_threshold(void) {
    double b;

    b = -1.0;
    /* ln(omega) <= 1.0816 -> no valid Molière distribution */
    ASSERT_TRUE(osh_physics_moliere_solve_b(2.0, &b) == 0);
    ASSERT_TRUE(osh_physics_moliere_solve_b(1.0, &b) == 0);
    ASSERT_TRUE(osh_physics_moliere_solve_b(0.0, &b) == 0);
}

/* ---- Reduced-angle functions f⁽ⁿ⁾ --------------------------------------- */

/* f⁽⁰⁾(ϑ) = 2·exp(-ϑ²) in closed form. */
static void test_reduced_f0_analytic(void) {
    double th;
    double got;
    double want;

    for (th = 0.0; th <= 4.0; th += 0.5) {
        got = osh_physics_moliere_reduced_f(0, th);
        want = 2.0 * exp(-th * th);
        ASSERT_TRUE(fabs(got - want) < 5e-3);
    }
}

/* Normalisation is carried entirely by f⁽⁰⁾: ∫ ϑ·f⁽⁰⁾ dϑ = 1 and
 * ∫ ϑ·f⁽ⁿ⁾ dϑ = 0 for n = 1, 2 (so the total is normalised for any B). */
static void test_reduced_f_normalization(void) {
    double dth;
    double th;
    double m0;
    double m1;
    double m2;

    dth = 0.01;
    m0 = 0.0;
    m1 = 0.0;
    m2 = 0.0;
    for (th = 0.5 * dth; th < 12.0; th += dth) {
        m0 += th * osh_physics_moliere_reduced_f(0, th) * dth;
        m1 += th * osh_physics_moliere_reduced_f(1, th) * dth;
        m2 += th * osh_physics_moliere_reduced_f(2, th) * dth;
    }
    ASSERT_TRUE(fabs(m0 - 1.0) < 1e-2);
    ASSERT_TRUE(fabs(m1) < 1e-2);
    ASSERT_TRUE(fabs(m2) < 1e-2);
}

/* ---- Sampling moments ---------------------------------------------------- */

static void test_sample_reduced_moments(void) {
    struct osh_rng rng;
    double b;
    double sum_sq;
    double th;
    double mean_sq;
    double tail_frac;
    int n;
    int i;
    long tail;

    n = 400000;
    sum_sq = 0.0;
    tail = 0;
    ASSERT_TRUE(osh_physics_moliere_solve_b(1.0e6, &b) == 1);
    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, 12345u, 67u);

    for (i = 0; i < n; ++i) {
        th = osh_physics_moliere_sample_reduced(b, &rng);
        ASSERT_TRUE(th >= 0.0);
        sum_sq += th * th;
        if (th > 2.5) {
            ++tail;
        }
    }
    /* Core variance ≈ 1 (Gaussian part) with a small positive tail correction. */
    mean_sq = sum_sq / (double) n;
    tail_frac = (double) tail / (double) n;
    ASSERT_TRUE(mean_sq > 0.8 && mean_sq < 3.0);
    /* Heavier tail than a pure Gaussian (P(ϑ>2.5) = e^-6.25 ≈ 0.0019). */
    ASSERT_TRUE(tail_frac > 0.0019);
    ASSERT_TRUE(tail_frac < 0.2);
}

/* The precomputed inverse-CDF must be monotone non-decreasing in u and bounded
 * to the tabulated reduced-angle range for every B. */
static void test_inv_cdf_monotone(void) {
    double bvals[] = {2.0, 5.0, 12.0, 30.0};
    double prev;
    double th;
    double u;
    size_t bi;
    int k;

    for (bi = 0; bi < sizeof(bvals) / sizeof(bvals[0]); ++bi) {
        prev = osh_physics_moliere_inv_cdf(bvals[bi], 0.0);
        ASSERT_TRUE(prev >= 0.0 && prev < 1e-6); /* u=0 -> theta=0 */
        for (k = 1; k <= 1000; ++k) {
            u = (double) k / 1000.0;
            th = osh_physics_moliere_inv_cdf(bvals[bi], u);
            ASSERT_TRUE(th >= prev - 1e-9);  /* monotone non-decreasing */
            ASSERT_TRUE(th >= 0.0 && th <= 12.0); /* within [0, THMAX] */
            prev = th;
        }
    }
}

/* ---- Dispatch: osh_physics_mcs_scatter ----------------------------------- */

static int is_unit(double const w[3]) {
    double n;

    n = sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
    return fabs(n - 1.0) < 1e-9;
}

static void test_mcs_scatter_moliere_deflects(void) {
    struct osh_rng rng;
    double v[3];
    double w[3];
    /* water-like medium constants at 200 MeV proton */
    double chic2; /* 0.157 * Sum w_i Z_i(Z_i+1)/A_i */
    double screen_z;
    double x0;
    double sum_cos;
    int n;
    int i;
    int deflected;
    int rc;

    v[0] = 0.0;
    v[1] = 0.0;
    v[2] = 1.0;
    chic2 = 0.6625;
    screen_z = 7.18;
    x0 = 36.08;
    n = 20000;
    deflected = 0;
    sum_cos = 0.0;
    osh_physics_moliere_init();
    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, 999u, 1u);

    for (i = 0; i < n; ++i) {
        rc = osh_physics_mcs_scatter(
            OSH_MCS_MOLIERE, v, w, 200.0, 938.272, 1.0, 1.0 /*d*/, 26.0 /*path*/, x0, chic2, screen_z, &rng);
        ASSERT_TRUE(is_unit(w));
        if (rc == 1) {
            ++deflected;
        }
        sum_cos += w[2];
    }
    /* Most steps should deflect, and the mean forward cosine stays high but < 1. */
    ASSERT_TRUE(deflected > n / 2);
    ASSERT_TRUE(sum_cos / (double) n > 0.9);
    ASSERT_TRUE(sum_cos / (double) n < 1.0);
}

static void test_mcs_scatter_gaussian_unit(void) {
    struct osh_rng rng;
    double v[3];
    double w[3];
    int rc;

    v[0] = 0.0;
    v[1] = 0.0;
    v[2] = 1.0;
    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, 7u, 3u);
    rc = osh_physics_mcs_scatter(OSH_MCS_GAUSSIAN, v, w, 200.0, 938.272, 1.0, 1.0, 26.0, 36.08, 0.6625, 7.18, &rng);
    ASSERT_TRUE(rc == 1);
    ASSERT_TRUE(is_unit(w));
}

static void test_mcs_scatter_off_copies(void) {
    struct osh_rng rng;
    double v[3];
    double w[3];
    int rc;

    v[0] = 0.0;
    v[1] = 0.0;
    v[2] = 1.0;
    w[0] = 9.0;
    w[1] = 9.0;
    w[2] = 9.0;
    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, 1u, 1u);
    rc = osh_physics_mcs_scatter(OSH_MCS_OFF, v, w, 200.0, 938.272, 1.0, 1.0, 26.0, 36.08, 0.6625, 7.18, &rng);
    ASSERT_TRUE(rc == 0);
    ASSERT_TRUE(w[0] == v[0] && w[1] == v[1] && w[2] == v[2]);

    rc = osh_physics_mcs_scatter(OSH_MCS_WENTZEL, v, w, 200.0, 938.272, 1.0, 1.0, 26.0, 36.08, 0.6625, 7.18, &rng);
    ASSERT_TRUE(rc == 0);
    ASSERT_TRUE(w[0] == v[0] && w[1] == v[1] && w[2] == v[2]);
}

static int run_named_test(char const *name) {
    if (strcmp(name, "theta0_path_scale_additivity") == 0) {
        test_theta0_path_scale_additivity();
        return 0;
    }
    if (strcmp(name, "theta0_fallback_to_thickness") == 0) {
        test_theta0_fallback_to_thickness();
        return 0;
    }
    if (strcmp(name, "solve_b_equation") == 0) {
        test_solve_b_equation();
        return 0;
    }
    if (strcmp(name, "solve_b_threshold") == 0) {
        test_solve_b_threshold();
        return 0;
    }
    if (strcmp(name, "reduced_f0_analytic") == 0) {
        test_reduced_f0_analytic();
        return 0;
    }
    if (strcmp(name, "reduced_f_normalization") == 0) {
        test_reduced_f_normalization();
        return 0;
    }
    if (strcmp(name, "sample_reduced_moments") == 0) {
        test_sample_reduced_moments();
        return 0;
    }
    if (strcmp(name, "inv_cdf_monotone") == 0) {
        test_inv_cdf_monotone();
        return 0;
    }
    if (strcmp(name, "mcs_scatter_moliere_deflects") == 0) {
        test_mcs_scatter_moliere_deflects();
        return 0;
    }
    if (strcmp(name, "mcs_scatter_gaussian_unit") == 0) {
        test_mcs_scatter_gaussian_unit();
        return 0;
    }
    if (strcmp(name, "mcs_scatter_off_copies") == 0) {
        test_mcs_scatter_off_copies();
        return 0;
    }
    return 1;
}

int main(int argc, char *argv[]) {
    if (argc == 2) {
        return run_named_test(argv[1]);
    }
    test_theta0_path_scale_additivity();
    test_theta0_fallback_to_thickness();
    test_solve_b_equation();
    test_solve_b_threshold();
    test_reduced_f0_analytic();
    test_reduced_f_normalization();
    test_sample_reduced_moments();
    test_inv_cdf_monotone();
    test_mcs_scatter_moliere_deflects();
    test_mcs_scatter_gaussian_unit();
    test_mcs_scatter_off_copies();
    printf("All osh_physics_scat tests passed.\n");
    return 0;
}
