#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "openshieldhit/status.h"
#include "physics/neutron/osh_neutron_reaction.h"
#include "physics/neutron/osh_neutron_xsec.h"
#include "physics/nuclear/osh_nuclear_handler.h"
#include "random/osh_rng.h"

#define ASSERT_TRUE(cond)                                                                                              \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            fprintf(stderr, "ASSERT FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__);                                 \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

#define ASSERT_NEAR(a, b, tol) ASSERT_TRUE(fabs((a) - (b)) <= (tol))
#define N_SAMPLES 500

/*
 * Build a minimal single-material handler for one element.
 * The fbu field is zero-initialised; it is safe as long as the test energies
 * keep the compound channel probability negligibly small (H-1 at 5 MeV has
 * σ(n,n') ≈ 0; O-16 at 1 MeV is below the (n,n') threshold of ~6 MeV).
 */
static void make_handler(struct osh_nuclear_handler *h,
                         struct osh_nuclear_elem *elem,
                         size_t *offset,
                         size_t *count,
                         unsigned int z,
                         unsigned int a) {
    memset(h, 0, sizeof(*h));
    elem->z = z;
    elem->a = a;
    elem->mass_fraction = 1.0f;
    *offset = 0u;
    *count = 1u;
    h->elem_pool = elem;
    h->elem_offset = offset;
    h->elem_count = count;
}

/*
 * H-1 elastic at 5 MeV: proton secondary always present; its kinetic energy
 * plus the scattered neutron energy must equal the incident energy within the
 * non-relativistic approximation (error < 0.5% at 5 MeV/938 MeV).
 */
static void test_elastic_h1_secondary_energy(void) {
    struct osh_neutron_xsec xsec;
    struct osh_nuclear_handler handler;
    struct osh_nuclear_elem elem;
    size_t offset, count;
    struct osh_rng rng;
    struct osh_neutron_reaction_event ev;
    double dir[3] = {0.0, 0.0, 1.0};
    double e_mev = 5.0;
    int i, n_elastic = 0;

    ASSERT_TRUE(osh_neutron_xsec_compile(NULL, &xsec) == OSH_OK);
    make_handler(&handler, &elem, &offset, &count, 1u, 1u);
    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, 1u, 0u);

    for (i = 0; i < N_SAMPLES; ++i) {
        osh_neutron_reaction_sample(&xsec, &handler, 0u, 1.0, e_mev, dir, &rng, &ev);
        if (ev.kind != OSH_NEUTRON_REACTION_ELASTIC) {
            continue;
        }
        ++n_elastic;
        ASSERT_TRUE(ev.n_secondaries == 1u);
        ASSERT_TRUE(ev.neutron_e_mev >= 0.0 && ev.neutron_e_mev <= e_mev);
        /* proton + neutron energies sum to incident (non-relativistic; 1% tol) */
        ASSERT_NEAR(ev.neutron_e_mev + ev.secondaries[0].energy, e_mev, 0.01 * e_mev);
    }

    ASSERT_TRUE(n_elastic > N_SAMPLES / 2); /* elastic should dominate at 5 MeV */
    osh_neutron_xsec_free(&xsec);
}

/*
 * Natural O (A=0 resolved to O-16) elastic at 1 MeV: scattered neutron energy
 * must fall in the kinematic window [E_min, E_n] where
 * E_min = E_n × ((A−1)/(A+1))².
 * For A = 16: E_min/E_n = (15/17)² ≈ 0.780.
 */
static void test_elastic_o16_energy_bounds(void) {
    struct osh_neutron_xsec xsec;
    struct osh_nuclear_handler handler;
    struct osh_nuclear_elem elem;
    size_t offset, count;
    struct osh_rng rng;
    struct osh_neutron_reaction_event ev;
    double dir[3] = {0.0, 0.0, 1.0};
    double e_mev = 1.0;
    double e_min = e_mev * (15.0 / 17.0) * (15.0 / 17.0); /* ≈ 0.780 MeV */
    int i, n_elastic = 0;

    ASSERT_TRUE(osh_neutron_xsec_compile(NULL, &xsec) == OSH_OK);
    make_handler(&handler, &elem, &offset, &count, 8u, 0u);
    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, 2u, 0u);

    for (i = 0; i < N_SAMPLES; ++i) {
        osh_neutron_reaction_sample(&xsec, &handler, 0u, 1.0, e_mev, dir, &rng, &ev);
        if (ev.kind != OSH_NEUTRON_REACTION_ELASTIC) {
            continue;
        }
        ++n_elastic;
        ASSERT_TRUE(ev.neutron_e_mev >= e_min - 1e-9);
        ASSERT_TRUE(ev.neutron_e_mev <= e_mev + 1e-9);
    }

    ASSERT_TRUE(n_elastic > N_SAMPLES / 2);
    osh_neutron_xsec_free(&xsec);
}

/*
 * After any elastic scatter the neutron direction must remain a unit vector.
 */
static void test_elastic_dir_normalized(void) {
    struct osh_neutron_xsec xsec;
    struct osh_nuclear_handler handler;
    struct osh_nuclear_elem elem;
    size_t offset, count;
    struct osh_rng rng;
    struct osh_neutron_reaction_event ev;
    double dir[3] = {0.0, 0.0, 1.0};
    double norm2;
    int i;

    ASSERT_TRUE(osh_neutron_xsec_compile(NULL, &xsec) == OSH_OK);
    make_handler(&handler, &elem, &offset, &count, 1u, 1u);
    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, 3u, 0u);

    for (i = 0; i < N_SAMPLES; ++i) {
        osh_neutron_reaction_sample(&xsec, &handler, 0u, 1.0, 5.0, dir, &rng, &ev);
        if (ev.kind != OSH_NEUTRON_REACTION_ELASTIC) {
            continue;
        }
        norm2 = ev.neutron_dir[0] * ev.neutron_dir[0] + ev.neutron_dir[1] * ev.neutron_dir[1]
                + ev.neutron_dir[2] * ev.neutron_dir[2];
        ASSERT_NEAR(norm2, 1.0, 1e-10);
    }
    osh_neutron_xsec_free(&xsec);
}

int main(void) {
    test_elastic_h1_secondary_energy();
    test_elastic_o16_energy_bounds();
    test_elastic_dir_normalized();
    return 0;
}
