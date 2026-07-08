/*
 * Tests for the evaluated-table reaction cross section (issue #277) and the
 * energy-dependent sigma_el/sigma_R ratio.
 *
 * The table values asserted here are the condensed ENDF/B-VIII.0 proton
 * sublibrary (LA150h) nonelastic curves committed in
 * src/physics/nuclear/osh_nuclear_sigma_reac_data.h; grid points are exact,
 * off-grid energies are lin-lin interpolated, and energies above the
 * evaluation's end (150 MeV) clamp flat (Renberg 1972: p+O 296 mb at
 * 230 MeV vs the 295 mb clamp).
 */
#include <math.h>
#include <stdio.h>

#include "physics/nuclear/osh_nuclear_elastic.h"
#include "physics/nuclear/osh_nuclear_sigma_reac.h"
#include "physics/nuclear/osh_nuclear_tripathi.h"
#include "test_assert.h"

#define MB 1.0e-27

static void test_table_hit_o16(void) {
    double sigma;

    /* Grid point: LA150 p+O-16 nonelastic at 100 MeV is 297 mb. */
    sigma = osh_nuclear_sigma_reac(1u, 1u, 8.0, 16.0, 100.0);
    ASSERT_TRUE(fabs(sigma / MB - 297.0) < 1.0);

    /* Above the evaluation's end: flat clamp at the 150 MeV value (295 mb). */
    sigma = osh_nuclear_sigma_reac(1u, 1u, 8.0, 16.0, 200.0);
    ASSERT_TRUE(fabs(sigma / MB - 295.0) < 1.0);

    /* Below the evaluation's threshold (7 MeV) the channel is closed. */
    sigma = osh_nuclear_sigma_reac(1u, 1u, 8.0, 16.0, 3.0);
    ASSERT_TRUE(sigma == 0.0);
}

static void test_table_hit_c12(void) {
    double sigma;

    /* Grid point: LA150 p+C-12 nonelastic at 100 MeV is 227 mb. */
    sigma = osh_nuclear_sigma_reac(1u, 1u, 6.0, 12.0, 100.0);
    ASSERT_TRUE(fabs(sigma / MB - 227.0) < 1.0);
}

static void test_interpolation_is_bounded(void) {
    double lo;
    double hi;
    double mid;

    lo = osh_nuclear_sigma_reac(1u, 1u, 8.0, 16.0, 40.0);
    hi = osh_nuclear_sigma_reac(1u, 1u, 8.0, 16.0, 41.0);
    mid = osh_nuclear_sigma_reac(1u, 1u, 8.0, 16.0, 40.5);
    ASSERT_TRUE((mid >= fmin(lo, hi)) && (mid <= fmax(lo, hi)));
}

static void test_fallback_is_tripathi(void) {
    double via_entry;
    double via_tripathi;

    /* Untabulated target (Al-27): the entry point must be exactly Tripathi. */
    via_entry = osh_nuclear_sigma_reac(1u, 1u, 13.0, 27.0, 100.0);
    via_tripathi = osh_nuclear_tripathi_sigma(1u, 1u, 13.0, 27.0, 100.0);
    ASSERT_TRUE(via_entry == via_tripathi);

    /* Non-proton projectile (alpha) on a tabulated target: Tripathi too. */
    via_entry = osh_nuclear_sigma_reac(2u, 4u, 8.0, 16.0, 100.0);
    via_tripathi = osh_nuclear_tripathi_sigma(2u, 4u, 8.0, 16.0, 100.0);
    ASSERT_TRUE(via_entry == via_tripathi);

    /* Compound-average pseudo-nucleus (non-integer A) never hits a table. */
    via_entry = osh_nuclear_sigma_reac(1u, 1u, 7.42, 14.89, 100.0);
    via_tripathi = osh_nuclear_tripathi_sigma(1u, 1u, 7.42, 14.89, 100.0);
    ASSERT_TRUE(via_entry == via_tripathi);
}

static void test_elastic_ratio_garron_anchor(void) {
    double sigma_el;

    /* Garron 1962: integrated nuclear elastic p+C-12 at 155 MeV = 75 +- 7 mb.
     * The default anchors must land within a generous window around it. */
    sigma_el = osh_nuclear_elastic_sigma(1u, 1u, 6.0, 12.0, 155.0);
    ASSERT_TRUE(sigma_el / MB > 55.0);
    ASSERT_TRUE(sigma_el / MB < 100.0);
}

static void test_elastic_ratio_monotone(void) {
    double sig_r_30;
    double sig_el_30;
    double r30;
    double r100;
    double r200;

    /* At and below the low anchor the ratio is the black-disk scale (1.0). */
    sig_r_30 = osh_nuclear_sigma_reac(1u, 1u, 8.0, 16.0, 30.0);
    sig_el_30 = osh_nuclear_elastic_sigma(1u, 1u, 8.0, 16.0, 30.0);
    ASSERT_TRUE(fabs(sig_el_30 / sig_r_30 - 1.0) < 1.0e-12);

    /* The ratio decreases with energy and settles at the high anchor. */
    r30 = sig_el_30 / sig_r_30;
    r100 = osh_nuclear_elastic_sigma(1u, 1u, 8.0, 16.0, 100.0) / osh_nuclear_sigma_reac(1u, 1u, 8.0, 16.0, 100.0);
    r200 = osh_nuclear_elastic_sigma(1u, 1u, 8.0, 16.0, 200.0) / osh_nuclear_sigma_reac(1u, 1u, 8.0, 16.0, 200.0);
    ASSERT_TRUE(r100 < r30);
    ASSERT_TRUE(r200 < r100);
    ASSERT_TRUE(fabs(r200 - 0.32) < 1.0e-12);
}

int main(void) {
    test_table_hit_o16();
    test_table_hit_c12();
    test_interpolation_is_bounded();
    test_fallback_is_tripathi();
    test_elastic_ratio_garron_anchor();
    test_elastic_ratio_monotone();
    printf("test_osh_nuclear_sigma_reac: all tests passed\n");
    return 0;
}
