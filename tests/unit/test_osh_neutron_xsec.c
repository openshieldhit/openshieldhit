#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "openshieldhit/status.h"
#include "physics/neutron/osh_neutron_xsec.h"

#define ASSERT_TRUE(cond)                                                                                              \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            fprintf(stderr, "ASSERT FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__);                                 \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

#define ASSERT_NEAR(a, b, tol) ASSERT_TRUE(fabs((a) - (b)) <= (tol))

/* H-1: σ_tot ≥ σ_el at all energies; σ_el > 0; σ_nonel = σ_tot − σ_el. */
static void test_h1_channels_consistent(void) {
    struct osh_neutron_xsec xsec;
    struct osh_neutron_xsec_result r;
    double energies[] = {1e-3, 1e-2, 0.1, 1.0, 5.0, 14.0};
    int i;

    ASSERT_TRUE(osh_neutron_xsec_compile(NULL, &xsec) == OSH_OK);

    for (i = 0; i < 6; ++i) {
        osh_neutron_xsec_lookup(&xsec, 1, 1, energies[i], &r);
        ASSERT_TRUE(r.el > 0.0);
        ASSERT_TRUE(r.tot >= r.el);
        ASSERT_NEAR(r.nonel, r.tot - r.el, 1e-6 * r.tot);
    }
    osh_neutron_xsec_free(&xsec);
}

/*
 * B-10(n,α) follows 1/v below ~100 keV: σ(E₁)/σ(E₂) = sqrt(E₂/E₁).
 * Check two adjacent decades on our thermal grid (1 meV and 10 meV).
 */
static void test_b10_na_one_over_v(void) {
    struct osh_neutron_xsec xsec;
    struct osh_neutron_xsec_result r1, r2;
    double e1 = 1.000e-9; /* 1 meV */
    double e2 = 1.000e-8; /* 10 meV */
    double ratio_got, ratio_expected;

    ASSERT_TRUE(osh_neutron_xsec_compile(NULL, &xsec) == OSH_OK);
    osh_neutron_xsec_lookup(&xsec, 5, 10, e1, &r1);
    osh_neutron_xsec_lookup(&xsec, 5, 10, e2, &r2);

    ASSERT_TRUE(r1.na > 0.0 && r2.na > 0.0);
    ratio_got = r1.na / r2.na;
    ratio_expected = sqrt(e2 / e1);                                /* = sqrt(10) ≈ 3.162 */
    ASSERT_NEAR(ratio_got, ratio_expected, 0.03 * ratio_expected); /* 3% */

    osh_neutron_xsec_free(&xsec);
}

/*
 * A nuclide not in the JEFF-4.0 table (Ta-181, Z=73) must fall back to the
 * optical/Tripathi model and return σ_tot > 0.
 */
static void test_optical_fallback_gives_nonzero_tot(void) {
    struct osh_neutron_xsec xsec;
    struct osh_neutron_xsec_result r;

    ASSERT_TRUE(osh_neutron_xsec_compile(NULL, &xsec) == OSH_OK);
    osh_neutron_xsec_lookup(&xsec, 73, 181, 10.0, &r);
    ASSERT_TRUE(r.tot > 0.0);
    osh_neutron_xsec_free(&xsec);
}

static void test_natural_element_resolves_to_representative_isotope(void) {
    struct osh_neutron_xsec xsec;
    struct osh_neutron_xsec_result natural_o;
    struct osh_neutron_xsec_result o16;

    ASSERT_TRUE(osh_neutron_xsec_compile(NULL, &xsec) == OSH_OK);
    ASSERT_TRUE(osh_neutron_xsec_resolve_a(8u, 0u) == 16u);

    osh_neutron_xsec_lookup(&xsec, 8, 0, 1.0, &natural_o);
    osh_neutron_xsec_lookup(&xsec, 8, 16, 1.0, &o16);
    ASSERT_NEAR(natural_o.tot, o16.tot, 1e-12 * o16.tot);
    ASSERT_NEAR(natural_o.el, o16.el, 1e-12 * o16.el);

    osh_neutron_xsec_free(&xsec);
}

int main(void) {
    test_h1_channels_consistent();
    test_b10_na_one_over_v();
    test_optical_fallback_gives_nonzero_tot();
    test_natural_element_resolves_to_representative_isotope();
    return 0;
}
