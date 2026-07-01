/* Validate the runtime universal-Landau λ(u) evaluator against the committed
 * golden fixture (tests/fixtures/vavilov/landau_ppf.csv), our clean-room
 * DENLAN reference. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "physics/atomic/osh_physics_strag_landau.h"
#include "test_assert.h"

#ifndef OSH_TEST_FIXTURES_DIR
#define OSH_TEST_FIXTURES_DIR "tests/fixtures"
#endif

/* Fit worst is ~0.01 λ; 0.05 leaves headroom. */
#define LAN_TOL 0.05

static void test_landau_matches_fixture(void) {
    char path[1024];
    FILE *fp;
    char line[256];
    double u;
    double lam_ref;
    double worst = 0.0;
    long checked = 0;

    snprintf(path, sizeof(path), "%s/vavilov/landau_ppf.csv", OSH_TEST_FIXTURES_DIR);
    fp = fopen(path, "r");
    ASSERT_TRUE(fp != NULL);

    while (fgets(line, sizeof(line), fp)) {
        double lam;
        double dev;
        if (line[0] == '#' || line[0] == '\n') {
            continue;
        }
        if (sscanf(line, "%lf %lf", &u, &lam_ref) != 2) {
            continue;
        }
        lam = osh_physics_strag_landau_lambda(u);
        dev = fabs(lam - lam_ref);
        if (dev > worst) {
            worst = dev;
        }
        ASSERT_TRUE(dev < LAN_TOL);
        ++checked;
    }
    fclose(fp);
    ASSERT_TRUE(checked > 30);
    printf("landau: %ld fixture rows, worst |dlambda| = %.4f\n", checked, worst);
}

static void test_landau_monotone_in_u(void) {
    double prev = -1e30;
    double u;
    for (u = 0.002; u <= 0.99; u += 0.002) {
        double lam = osh_physics_strag_landau_lambda(u);
        ASSERT_TRUE(lam > prev);
        prev = lam;
    }
}

static int run_named_test(char const *name) {
    if (strcmp(name, "landau_matches_fixture") == 0) {
        test_landau_matches_fixture();
        return 0;
    }
    if (strcmp(name, "landau_monotone_in_u") == 0) {
        test_landau_monotone_in_u();
        return 0;
    }
    return 1;
}

int main(int argc, char *argv[]) {
    if (argc == 2) {
        return run_named_test(argv[1]);
    }
    test_landau_matches_fixture();
    test_landau_monotone_in_u();
    printf("All osh_physics_strag_landau tests passed.\n");
    return 0;
}
