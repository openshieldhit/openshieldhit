#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "physics/osh_kinematics.h"
#include "random/osh_rng.h"

#define ASSERT_TRUE(cond)                                                                                              \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            fprintf(stderr, "ASSERT FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__);                                 \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

#define ASSERT_NEAR(a, b, tol) ASSERT_TRUE(fabs((a) - (b)) <= (tol))

/* ---- osh_kinematics_rotate_dir_cos --------------------------------------- */

static void test_rotate_identity(void) {
    double v[3] = {0.0, 0.0, 1.0};
    double w[3];
    double len;

    /* sin_theta == 0 → w must equal v */
    osh_kinematics_rotate_dir_cos(v, w, 1.0, 0.0, 1.0, 0.0);
    ASSERT_NEAR(w[0], v[0], 1e-12);
    ASSERT_NEAR(w[1], v[1], 1e-12);
    ASSERT_NEAR(w[2], v[2], 1e-12);

    /* Result must be unit length after any rotation */
    osh_kinematics_rotate_dir_cos(v, w, 0.5, sqrt(1.0 - 0.25), 1.0, 0.0);
    len = sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
    ASSERT_NEAR(len, 1.0, 1e-12);
}

static void test_rotate_unit_length(void) {
    double v[3];
    double w[3];
    double len;
    double cos_t, sin_t;
    int i;
    struct osh_rng rng;
    double cos_phi, sin_phi;

    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, 42u, 0u);

    for (i = 0; i < 1000; ++i) {
        /* Random incident direction */
        double u1 = osh_rng_double(&rng);
        double u2 = osh_rng_double(&rng);
        cos_t = 2.0 * u1 - 1.0;
        sin_t = sqrt(1.0 - cos_t * cos_t);

        double phi = 2.0 * 3.141592653589793 * u2;
        v[0] = sin_t * cos(phi);
        v[1] = sin_t * sin(phi);
        v[2] = cos_t;

        cos_t = 2.0 * osh_rng_double(&rng) - 1.0;
        sin_t = sqrt(1.0 - cos_t * cos_t);
        osh_kinematics_azimuth(&rng, &cos_phi, &sin_phi);

        osh_kinematics_rotate_dir_cos(v, w, cos_t, sin_t, cos_phi, sin_phi);
        len = sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
        ASSERT_NEAR(len, 1.0, 1e-12);
    }
}

static void test_rotate_aliasing(void) {
    double v[3] = {0.5773502691896258, 0.5773502691896258, 0.5773502691896258};
    double len;

    /* w may alias v — must not corrupt the result */
    osh_kinematics_rotate_dir_cos(v, v, 0.6, 0.8, 1.0, 0.0);
    len = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    ASSERT_NEAR(len, 1.0, 1e-12);
}

/* ---- osh_kinematics_azimuth ---------------------------------------------- */

static void test_azimuth_unit_circle(void) {
    struct osh_rng rng;
    double cos_phi, sin_phi;
    double r2;
    int i;

    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, 123u, 0u);
    for (i = 0; i < 10000; ++i) {
        osh_kinematics_azimuth(&rng, &cos_phi, &sin_phi);
        r2 = cos_phi * cos_phi + sin_phi * sin_phi;
        ASSERT_NEAR(r2, 1.0, 1e-12);
        ASSERT_TRUE(cos_phi >= -1.0 && cos_phi <= 1.0);
        ASSERT_TRUE(sin_phi >= -1.0 && sin_phi <= 1.0);
    }
}

/* ---- osh_kinematics_elastic_equal_mass_lab -------------------------------- */

#define PROTON_MASS_MEV 938.2720882

static void test_elastic_energy_conservation(void) {
    double e1, e2, cos1, cos2;
    double energies[] = {10.0, 100.0, 500.0, 1000.0, 2000.0};
    double cos_cms[] = {1.0, 0.0, -1.0, 0.5, -0.5};
    int i, j;

    for (i = 0; i < 5; ++i) {
        for (j = 0; j < 5; ++j) {
            osh_kinematics_elastic_equal_mass_lab(
                energies[i], PROTON_MASS_MEV, cos_cms[j], &cos1, &e1, &cos2, &e2);
            /* Energy conservation: e1 + e2 == T_lab within 1e-6 relative */
            ASSERT_NEAR(e1 + e2, energies[i], energies[i] * 1e-6);
        }
    }
}

static void test_elastic_nonrelativistic_limit(void) {
    /* At low T, θ1 + θ2 ≈ π/2: cos_theta1 + cos_theta2 ≈ 0
     * (exact identity for non-relativistic equal-mass elastic) */
    double e1, e2, cos1, cos2;
    double cos_cm = 0.5; /* arbitrary CM angle */

    osh_kinematics_elastic_equal_mass_lab(
        10.0, PROTON_MASS_MEV, cos_cm, &cos1, &e1, &cos2, &e2);

    /* Non-relativistic: cos1^2 + cos2^2 ≈ 1 (Pythagorean identity)
     * This holds because θ1 + θ2 = π/2 implies cos^2(θ1) + cos^2(θ2) = 1. */
    ASSERT_NEAR(cos1 * cos1 + cos2 * cos2, 1.0, 0.01);
}

static void test_elastic_recoil_always_forward(void) {
    double e1, e2, cos1, cos2;
    double cos_cm;
    int i;

    /* Recoil proton is always in the forward hemisphere for all CM angles.
     * Analytically exact: p2z_lab = γ*p_cm*(1-cos_theta_cm) >= 0.
     * Allow 1e-10 tolerance for floating-point rounding near cos_cm = 1. */
    for (i = 0; i <= 100; ++i) {
        cos_cm = -1.0 + 2.0 * i / 100.0;
        osh_kinematics_elastic_equal_mass_lab(
            500.0, PROTON_MASS_MEV, cos_cm, &cos1, &e1, &cos2, &e2);
        ASSERT_TRUE(cos2 >= -1e-10);
    }
}

static void test_elastic_energies_non_negative(void) {
    double e1, e2, cos1, cos2;
    double cos_cm;
    int i;

    for (i = 0; i <= 100; ++i) {
        cos_cm = -1.0 + 2.0 * i / 100.0;
        osh_kinematics_elastic_equal_mass_lab(
            200.0, PROTON_MASS_MEV, cos_cm, &cos1, &e1, &cos2, &e2);
        ASSERT_TRUE(e1 >= 0.0);
        ASSERT_TRUE(e2 >= 0.0);
    }
}

int main(void) {
    test_rotate_identity();
    test_rotate_unit_length();
    test_rotate_aliasing();
    test_azimuth_unit_circle();
    test_elastic_energy_conservation();
    test_elastic_nonrelativistic_limit();
    test_elastic_recoil_always_forward();
    test_elastic_energies_non_negative();
    return 0;
}
