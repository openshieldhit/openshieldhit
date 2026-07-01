/* Unit tests for the energy-straggling kinematics helpers (ξ, E_max, κ, λ̄) and
 * the Bohr Gaussian width.  Reference values are hand-computed from the analytic
 * definitions (proton, T = 200 MeV, water-like Z/A), not copied from any
 * tabulated third-party numbers. */

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "physics/atomic/osh_physics_strag.h"
#include "test_assert.h"

/* Relative closeness with an absolute floor of 1, for near-zero comparands. */
static int approx(double a, double b, double rtol) {
    double s = fabs(b);
    return fabs(a - b) <= rtol * (s > 1.0 ? s : 1.0);
}

/* E_max = 2 mₑ β²γ² / (1 + 2γ mₑ/M + (mₑ/M)²); proton 200 MeV → 0.4815 MeV. */
static void test_emax_proton_200mev(void) {
    double emax = osh_physics_strag_emax(200.0, 938.272);
    ASSERT_TRUE(approx(emax, 0.4815, 3.0e-3));
    ASSERT_TRUE(osh_physics_strag_emax(-1.0, 938.272) == 0.0);
    ASSERT_TRUE(osh_physics_strag_emax(200.0, 0.0) == 0.0);
}

/* ξ = (K/2)(Z/A)(z²/β²)·d; z=1, Z/A=0.5551, d=0.1 g/cm², β²=0.320525. */
static void test_xi_water_step(void) {
    double xi = osh_physics_strag_xi(1.0, 0.5551, 0.1, 0.320525);
    double xi2 = osh_physics_strag_xi(1.0, 0.5551, 0.2, 0.320525);
    ASSERT_TRUE(approx(xi, 0.026590, 3.0e-3));
    ASSERT_TRUE(approx(xi2, 2.0 * xi, 1.0e-9)); /* linear in thickness */
    ASSERT_TRUE(osh_physics_strag_xi(1.0, 0.5551, 0.1, 0.0) == 0.0);
    ASSERT_TRUE(osh_physics_strag_xi(0.0, 0.5551, 0.1, 0.320525) == 0.0);
}

/* κ = ξ / E_max → 0.0266 / 0.4815 ≈ 0.05522 (Vavilov regime). */
static void test_kappa(void) {
    double k = osh_physics_strag_kappa(0.026590, 0.481494);
    ASSERT_TRUE(approx(k, 0.055225, 3.0e-3));
    ASSERT_TRUE(osh_physics_strag_kappa(0.02, 0.0) == 0.0);
    ASSERT_TRUE(osh_physics_strag_kappa(0.0, 0.4) == 0.0);
}

/* λ̄ = −β² − ln κ + (γ_Euler − 1); β²=0.320525, κ=0.055225 → 2.1529. */
static void test_lambda_bar(void) {
    double lb = osh_physics_strag_lambda_bar(0.055225, 0.320525);
    ASSERT_TRUE(approx(lb, 2.1529, 1.0e-3));
    ASSERT_TRUE(osh_physics_strag_lambda_bar(0.0, 0.32) == 0.0);
}

/* Bohr σ ∝ sqrt(d); guards on non-physical inputs. */
static void test_sigma_guards_and_scaling(void) {
    double s1 = osh_physics_strag_sigma(1.0, 0.5551, 0.1);
    double s4 = osh_physics_strag_sigma(1.0, 0.5551, 0.4);
    ASSERT_TRUE(s1 > 0.0);
    ASSERT_TRUE(approx(s4, 2.0 * s1, 1.0e-9));
    ASSERT_TRUE(osh_physics_strag_sigma(0.0, 0.5551, 0.1) == 0.0);
    ASSERT_TRUE(osh_physics_strag_sigma(1.0, 0.5551, 0.0) == 0.0);
}

static int run_named_test(char const *name) {
    if (strcmp(name, "emax_proton_200mev") == 0) {
        test_emax_proton_200mev();
        return 0;
    }
    if (strcmp(name, "xi_water_step") == 0) {
        test_xi_water_step();
        return 0;
    }
    if (strcmp(name, "kappa") == 0) {
        test_kappa();
        return 0;
    }
    if (strcmp(name, "lambda_bar") == 0) {
        test_lambda_bar();
        return 0;
    }
    if (strcmp(name, "sigma_guards_and_scaling") == 0) {
        test_sigma_guards_and_scaling();
        return 0;
    }
    return 1;
}

int main(int argc, char *argv[]) {
    if (argc == 2) {
        return run_named_test(argv[1]);
    }
    test_emax_proton_200mev();
    test_xi_water_step();
    test_kappa();
    test_lambda_bar();
    test_sigma_guards_and_scaling();
    printf("All osh_physics_strag tests passed.\n");
    return 0;
}
