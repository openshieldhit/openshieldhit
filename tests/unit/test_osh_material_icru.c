/*
 * test_osh_material_icru.c
 *
 * Unit tests for osh_material_icru_lookup().  Values are cross-checked against
 * the NIST ESTAR/PSTAR material database (https://physics.nist.gov/Star/).
 *
 * Run individual cases via:
 *   test_osh_material_icru <test_name>
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "material/osh_material_icru.h"
#include "openshieldhit/status.h"

/* Tolerance for density/mean-excitation comparisons.
 * The embedded tables are stored as float (~6 significant digits). */
#define FTOL 1e-5

#define ASSERT_TRUE(cond)                                                                                              \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            fprintf(stderr, "ASSERT FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__);                                 \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

#define ASSERT_NEAR(a, b, tol) ASSERT_TRUE(fabs((double) (a) - (double) (b)) < (tol))

/* ---- Elements ------------------------------------------------------------ */

static void test_icru_hydrogen_gas(void) {
    struct osh_material_icru_entry e;

    ASSERT_TRUE(osh_material_icru_lookup(1, &e) == OSH_OK);
    ASSERT_TRUE(e.icru_id == 1);
    ASSERT_TRUE(e.nelements == 1u);
    ASSERT_TRUE(e.elements[0].z == 1u);
    ASSERT_TRUE(e.elements[0].a == 0u);
    ASSERT_NEAR(e.elements[0].mass_fraction, 1.0, FTOL);
    ASSERT_NEAR(e.rho, 8.3748e-5, 1e-8); /* NIST: 8.37480e-5 g/cm³ */
    ASSERT_NEAR(e.mean_excitation_energy, 19.2, FTOL);
    ASSERT_TRUE(e.state == OSH_MATERIAL_STATE_GAS);
}

static void test_icru_carbon_condensed(void) {
    struct osh_material_icru_entry e;

    ASSERT_TRUE(osh_material_icru_lookup(6, &e) == OSH_OK);
    ASSERT_TRUE(e.icru_id == 6);
    ASSERT_TRUE(e.nelements == 1u);
    ASSERT_TRUE(e.elements[0].z == 6u);
    ASSERT_TRUE(e.elements[0].a == 0u);
    ASSERT_NEAR(e.elements[0].mass_fraction, 1.0, FTOL);
    ASSERT_NEAR(e.rho, 2.0, FTOL); /* NIST: 2.0 g/cm³ */
    ASSERT_NEAR(e.mean_excitation_energy, 81.0, FTOL);
    ASSERT_TRUE(e.state == OSH_MATERIAL_STATE_CONDENSED);
}

static void test_icru_aluminum_condensed(void) {
    struct osh_material_icru_entry e;

    ASSERT_TRUE(osh_material_icru_lookup(13, &e) == OSH_OK);
    ASSERT_TRUE(e.icru_id == 13);
    ASSERT_TRUE(e.nelements == 1u);
    ASSERT_TRUE(e.elements[0].z == 13u);
    ASSERT_NEAR(e.rho, 2.6989, FTOL); /* NIST: 2.6989 g/cm³ */
    ASSERT_NEAR(e.mean_excitation_energy, 166.0, FTOL);
    ASSERT_TRUE(e.state == OSH_MATERIAL_STATE_CONDENSED);
}

/* ---- Compounds ----------------------------------------------------------- */

static void test_icru_air_dry(void) {
    struct osh_material_icru_entry e;

    /* NIST ESTAR id 104: Air, Dry (near sea level) */
    ASSERT_TRUE(osh_material_icru_lookup(104, &e) == OSH_OK);
    ASSERT_TRUE(e.icru_id == 104);
    ASSERT_TRUE(e.state == OSH_MATERIAL_STATE_GAS);
    ASSERT_NEAR(e.rho, 1.20479e-3, 1e-7); /* NIST: 1.20479e-3 g/cm³ */
    ASSERT_NEAR(e.mean_excitation_energy, 85.7, FTOL);
    /* Composition: C, N, O, Ar (4 elements) */
    ASSERT_TRUE(e.nelements == 4u);
    /* Nitrogen dominates by mass */
    {
        size_t i;
        double mfrac_n = -1.0;
        double mfrac_o = -1.0;
        double mfrac_sum = 0.0;
        for (i = 0; i < e.nelements; i++) {
            if (e.elements[i].z == 7u)
                mfrac_n = e.elements[i].mass_fraction;
            if (e.elements[i].z == 8u)
                mfrac_o = e.elements[i].mass_fraction;
            mfrac_sum += e.elements[i].mass_fraction;
        }
        ASSERT_NEAR(mfrac_n, 0.755267, FTOL); /* NIST: ~75.5% N */
        ASSERT_NEAR(mfrac_o, 0.231781, FTOL); /* NIST: ~23.2% O */
        ASSERT_NEAR(mfrac_sum, 1.0, 1e-4);    /* fractions sum to 1 */
    }
}

static void test_icru_pmma(void) {
    struct osh_material_icru_entry e;

    /* NIST ESTAR id 223: PMMA (poly methyl methacrylate / Lucite / Plexiglas)
     * Composition: H, C, O  (C5H8O2)n */
    ASSERT_TRUE(osh_material_icru_lookup(223, &e) == OSH_OK);
    ASSERT_TRUE(e.icru_id == 223);
    ASSERT_TRUE(e.state == OSH_MATERIAL_STATE_CONDENSED);
    ASSERT_NEAR(e.rho, 1.19, FTOL); /* NIST: 1.19 g/cm³ */
    ASSERT_NEAR(e.mean_excitation_energy, 74.0, FTOL);
    ASSERT_TRUE(e.nelements == 3u);
    {
        size_t i;
        double mfrac_h = -1.0;
        double mfrac_c = -1.0;
        double mfrac_o = -1.0;
        for (i = 0; i < e.nelements; i++) {
            if (e.elements[i].z == 1u)
                mfrac_h = e.elements[i].mass_fraction;
            if (e.elements[i].z == 6u)
                mfrac_c = e.elements[i].mass_fraction;
            if (e.elements[i].z == 8u)
                mfrac_o = e.elements[i].mass_fraction;
        }
        ASSERT_NEAR(mfrac_h, 0.080538, FTOL);
        ASSERT_NEAR(mfrac_c, 0.599848, FTOL);
        ASSERT_NEAR(mfrac_o, 0.319614, FTOL);
    }
}

static void test_icru_water(void) {
    struct osh_material_icru_entry e;

    /* NIST ESTAR id 276: Water, Liquid */
    ASSERT_TRUE(osh_material_icru_lookup(276, &e) == OSH_OK);
    ASSERT_TRUE(e.icru_id == 276);
    ASSERT_TRUE(e.state == OSH_MATERIAL_STATE_CONDENSED);
    ASSERT_NEAR(e.rho, 1.0, FTOL);
    ASSERT_NEAR(e.mean_excitation_energy, 75.0, FTOL);
    ASSERT_TRUE(e.nelements == 2u);
    {
        size_t i;
        double mfrac_h = -1.0;
        double mfrac_o = -1.0;
        for (i = 0; i < e.nelements; i++) {
            if (e.elements[i].z == 1u)
                mfrac_h = e.elements[i].mass_fraction;
            if (e.elements[i].z == 8u)
                mfrac_o = e.elements[i].mass_fraction;
        }
        ASSERT_NEAR(mfrac_h, 0.111894, FTOL);
        ASSERT_NEAR(mfrac_o, 0.888106, FTOL);
    }
}

/* ---- Edge cases ---------------------------------------------------------- */

static void test_icru_invalid_id_returns_error(void) {
    struct osh_material_icru_entry e;

    ASSERT_TRUE(osh_material_icru_lookup(0, &e) == OSH_ENOTSUP);
    ASSERT_TRUE(osh_material_icru_lookup(-1, &e) == OSH_ENOTSUP);
    ASSERT_TRUE(osh_material_icru_lookup(9999, &e) == OSH_ENOTSUP);
}

static void test_icru_null_entry_returns_error(void) {
    ASSERT_TRUE(osh_material_icru_lookup(276, NULL) == OSH_EINVAL);
}

static void test_icru_graphite_special_id(void) {
    struct osh_material_icru_entry e;

    /* Graphite uses the non-sequential external id 906 */
    ASSERT_TRUE(osh_material_icru_lookup(906, &e) == OSH_OK);
    ASSERT_TRUE(e.icru_id == 906);
    ASSERT_TRUE(e.nelements == 1u);
    ASSERT_TRUE(e.elements[0].z == 6u); /* pure carbon */
    ASSERT_NEAR(e.mean_excitation_energy, 78.0, FTOL);
    ASSERT_TRUE(e.state == OSH_MATERIAL_STATE_CONDENSED);
}

/* ---- Runner -------------------------------------------------------------- */

int main(int argc, char *argv[]) {
    /* If a test name is passed, run only that one (for CTest per-function mode). */
#define RUN(name)                                                                                                      \
    if (argc < 2 || strcmp(argv[1], #name) == 0) {                                                                     \
        name();                                                                                                        \
    }

    RUN(test_icru_hydrogen_gas)
    RUN(test_icru_carbon_condensed)
    RUN(test_icru_aluminum_condensed)
    RUN(test_icru_air_dry)
    RUN(test_icru_pmma)
    RUN(test_icru_water)
    RUN(test_icru_invalid_id_returns_error)
    RUN(test_icru_null_entry_returns_error)
    RUN(test_icru_graphite_special_id)

#undef RUN
    return 0;
}
