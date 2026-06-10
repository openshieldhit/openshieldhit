#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "openshieldhit/const.h"
#include "physics/nuclear/osh_nuclear_pp.h"
#include "random/osh_rng.h"

#define ASSERT_TRUE(cond)                                                                                              \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            fprintf(stderr, "ASSERT FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__);                                 \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

#define ASSERT_NEAR(a, b, tol) ASSERT_TRUE(fabs((a) - (b)) <= (tol))

/* ---- osh_nuclear_pp_sigma_el --------------------------------------------- */

static void test_sigma_el_known_values(void) {
    double sigma;

    /* Below table: must be zero */
    sigma = osh_nuclear_pp_sigma_el(5.0);
    ASSERT_TRUE(sigma == 0.0);

    /* At 100 MeV: σ_tot = 32.8 mb → 32.8e-27 cm² */
    sigma = osh_nuclear_pp_sigma_el(100.0);
    ASSERT_NEAR(sigma, 32.8 * OSH_MB_TO_CM2, OSH_MB_TO_CM2);

    /* At 5000 MeV (last bin): σ_tot = 33.603 mb → clamp returns last value */
    sigma = osh_nuclear_pp_sigma_el(5000.0);
    ASSERT_NEAR(sigma, 33.603 * OSH_MB_TO_CM2, OSH_MB_TO_CM2);

    /* Well above table: should return last-bin value (te clamped to 1) */
    double sigma_clamped = osh_nuclear_pp_sigma_el(9999.0);
    ASSERT_NEAR(sigma_clamped, sigma, OSH_MB_TO_CM2);
}

static void test_sigma_el_positive(void) {
    int i;
    for (i = 1; i <= 100; ++i) {
        double e = 10.0 + (5000.0 - 10.0) * i / 100.0;
        ASSERT_TRUE(osh_nuclear_pp_sigma_el(e) > 0.0);
    }
}

/* ---- osh_nuclear_pp_sample_cos_theta_cm ---------------------------------- */

static void test_sample_cos_theta_cm_range(void) {
    struct osh_rng rng;
    double cos_theta;
    int i;

    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, 42u, 0u);
    for (i = 0; i < 1000; ++i) {
        cos_theta = osh_nuclear_pp_sample_cos_theta_cm(200.0, &rng);
        ASSERT_TRUE(cos_theta >= -1.0 && cos_theta <= 1.0);
    }
}

static void test_sample_cos_theta_cm_multiple_energies(void) {
    struct osh_rng rng;
    double cos_theta;
    double energies[] = {10.0, 50.0, 100.0, 300.0, 1000.0, 5000.0};
    int ne = 6;
    int i, j;

    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, 99u, 0u);
    for (i = 0; i < ne; ++i) {
        for (j = 0; j < 200; ++j) {
            cos_theta = osh_nuclear_pp_sample_cos_theta_cm(energies[i], &rng);
            ASSERT_TRUE(cos_theta >= -1.0 && cos_theta <= 1.0);
        }
    }
}

/* ---- osh_nuclear_pp_lambda_gcm2 ------------------------------------------ */

static void test_lambda_boundary_conditions(void) {
    double sigma = osh_nuclear_pp_sigma_el(100.0);

    /* Zero hydrogen fraction → large mean free path */
    ASSERT_TRUE(osh_nuclear_pp_lambda_gcm2(0.0, sigma) >= 1.0e29);

    /* Zero sigma → large mean free path */
    ASSERT_TRUE(osh_nuclear_pp_lambda_gcm2(1.0, 0.0) >= 1.0e29);

    /* Pure hydrogen at 100 MeV: λ ≈ 1.00794 / (N_A × σ) */
    double lambda = osh_nuclear_pp_lambda_gcm2(1.0, sigma);
    ASSERT_TRUE(lambda > 0.0 && lambda < 1.0e10);

    /* Half hydrogen → double mean free path */
    double lambda_half = osh_nuclear_pp_lambda_gcm2(0.5, sigma);
    ASSERT_NEAR(lambda_half / lambda, 2.0, 0.01);
}

int main(void) {
    test_sigma_el_known_values();
    test_sigma_el_positive();
    test_sample_cos_theta_cm_range();
    test_sample_cos_theta_cm_multiple_energies();
    test_lambda_boundary_conditions();
    return 0;
}
