/* Validate the runtime Vavilov λ(u; κ, β²) evaluator against the committed
 * golden inverse-CDF fixture (tests/fixtures/vavilov/vavilov_ppf.csv), which was
 * generated from the exact Vavilov distribution at (κ, β²) nodes that are
 * independent of the fit's cache — so this also exercises interpolation. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "physics/atomic/osh_physics_strag_vavilov.h"
#include "test_assert.h"

#ifndef OSH_TEST_FIXTURES_DIR
#define OSH_TEST_FIXTURES_DIR "tests/fixtures"
#endif

/* Max |λ_fit − λ_ref| tolerance.  After the reference-generation fix the fit is
 * accurate to <1e-3 λ where transport operates (κ≳0.1) and to ~0.07 for u≤0.99
 * across the whole grid; the worst ~0.20 is confined to the extreme high-loss
 * tail (u∈(0.99,0.995]) at small κ — a rare×rare, sub-keV region at the edge of
 * the sampler's range that the polynomial+Chebyshev form cannot fit tighter.
 * There OSH is still an order of magnitude closer to the exact Vavilov than the
 * reference libvav (verified: OSH ≤0.2 λ vs libvav 5–19 λ, the latter being
 * outside its own stated validity u≤0.990 for κ<0.04).  0.25 covers the worst
 * with headroom while still guarding against gross coefficient regressions. */
#define VAV_TOL 0.25

static void test_vavilov_matches_fixture(void) {
    char path[1024];
    FILE *fp;
    char line[512];
    double kappa;
    double beta2;
    double u;
    double lam_ref;
    double worst = 0.0;
    long checked = 0;

    snprintf(path, sizeof(path), "%s/vavilov/vavilov_ppf.csv", OSH_TEST_FIXTURES_DIR);
    fp = fopen(path, "r");
    ASSERT_TRUE(fp != NULL);

    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#' || line[0] == '\n') {
            continue;
        }
        if (sscanf(line, "%lf %lf %lf %lf", &kappa, &beta2, &u, &lam_ref) != 4) {
            continue;
        }
        double lam = osh_physics_strag_vavilov_lambda(kappa, beta2, u);
        double dev = fabs(lam - lam_ref);
        if (dev > worst) {
            worst = dev;
        }
        ASSERT_TRUE(dev < VAV_TOL);
        ++checked;
    }
    fclose(fp);
    ASSERT_TRUE(checked > 100);
    printf("vavilov: %ld fixture rows, worst |dlambda| = %.4f\n", checked, worst);
}

/* Sampler must be monotone increasing in u (a valid inverse CDF). */
static void test_vavilov_monotone_in_u(void) {
    double kappas[] = {0.02, 0.1, 0.5, 2.0, 8.0};
    double betas[] = {0.1, 0.5, 0.9};
    size_t ki;
    size_t bi;

    for (ki = 0; ki < sizeof(kappas) / sizeof(kappas[0]); ++ki) {
        for (bi = 0; bi < sizeof(betas) / sizeof(betas[0]); ++bi) {
            double prev = -1e30;
            double u;
            for (u = 0.005; u <= 0.99; u += 0.005) {
                double lam = osh_physics_strag_vavilov_lambda(kappas[ki], betas[bi], u);
                ASSERT_TRUE(lam > prev);
                prev = lam;
            }
        }
    }
}

static void test_vavilov_log_kappa_matches_kappa(void) {
    double kappas[] = {0.02, 0.1, 0.5, 2.0, 8.0};
    double betas[] = {0.1, 0.5, 0.9};
    size_t ki;
    size_t bi;
    double by_kappa;
    double by_log;

    for (ki = 0; ki < sizeof(kappas) / sizeof(kappas[0]); ++ki) {
        for (bi = 0; bi < sizeof(betas) / sizeof(betas[0]); ++bi) {
            by_kappa = osh_physics_strag_vavilov_lambda(kappas[ki], betas[bi], 0.42);
            by_log = osh_physics_strag_vavilov_lambda_log(log(kappas[ki]), betas[bi], 0.42);
            ASSERT_TRUE(fabs(by_kappa - by_log) < 1.0e-12);
        }
    }
}

static int run_named_test(char const *name) {
    if (strcmp(name, "vavilov_matches_fixture") == 0) {
        test_vavilov_matches_fixture();
        return 0;
    }
    if (strcmp(name, "vavilov_monotone_in_u") == 0) {
        test_vavilov_monotone_in_u();
        return 0;
    }
    if (strcmp(name, "vavilov_log_kappa_matches_kappa") == 0) {
        test_vavilov_log_kappa_matches_kappa();
        return 0;
    }
    return 1;
}

int main(int argc, char *argv[]) {
    if (argc == 2) {
        return run_named_test(argv[1]);
    }
    test_vavilov_matches_fixture();
    test_vavilov_monotone_in_u();
    test_vavilov_log_kappa_matches_kappa();
    printf("All osh_physics_strag_vavilov tests passed.\n");
    return 0;
}
