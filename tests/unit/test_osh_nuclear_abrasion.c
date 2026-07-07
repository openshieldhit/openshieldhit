#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "openshieldhit/const.h"
#include "particle/osh_particle.h"
#include "particle/osh_particle_const.h"
#include "physics/nuclear/osh_nuclear_abrasion.h"
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

/*
 * O-16 at 140 MeV: Tripathi sigma_R is approximately 450 mb.
 * Expected mean participants: sigma_pN * A / sigma_pA = 30mb * 16 / 450mb ≈ 1.07
 */
#define O16_A 16.0
#define O16_Z 8.0
#define SIGMA_PA_CM2 (450.0 * OSH_MB_TO_CM2)

static void test_event_kind(void) {
    struct osh_rng rng;
    struct osh_nuclear_event ev;
    double dir[3];
    int i;

    dir[0] = 0.0;
    dir[1] = 0.0;
    dir[2] = 1.0;

    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, 55u, 0u);
    for (i = 0; i < 100; ++i) {
        osh_nuclear_abrasion_step(140.0, dir, O16_A, O16_Z, SIGMA_PA_CM2, &rng, &ev);
        ASSERT_TRUE(ev.kind == OSH_NUCLEAR_EVENT_ABRASION);
        ASSERT_TRUE(ev.primary_energy == 0.0);
        ASSERT_TRUE(ev.n_fragments == 1u);
        ASSERT_TRUE(ev.fragments[0].a <= (unsigned int) O16_A);
        ASSERT_TRUE(ev.fragments[0].z <= (unsigned int) O16_Z);
    }
}

static void test_energy_conservation(void) {
    /* Exact kinetic-energy budget: T_in = sum(KE_secondaries) + E*.
     * The knocked-out nucleons, the escaping cascade proton, and the
     * prefragment excitation account for the full incident energy. */
    struct osh_rng rng;
    struct osh_nuclear_event ev;
    double dir[3];
    double T_incident;
    double sum_e;
    size_t j;
    int i;

    dir[0] = 0.0;
    dir[1] = 0.0;
    dir[2] = 1.0;
    T_incident = 140.0;

    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, 99u, 0u);
    for (i = 0; i < 1000; ++i) {
        sum_e = 0.0;
        osh_nuclear_abrasion_step(T_incident, dir, O16_A, O16_Z, SIGMA_PA_CM2, &rng, &ev);
        for (j = 0u; j < ev.n_secondaries; ++j) {
            sum_e += ev.secondaries[j].energy;
        }
        ASSERT_TRUE(sum_e <= T_incident + 1.0e-6);
        if (ev.n_fragments > 0u) {
            ASSERT_NEAR(sum_e + ev.fragments[0].excitation_energy, T_incident, 1.0e-6);
        }
    }
}

static void test_secondary_dir_unit(void) {
    struct osh_rng rng;
    struct osh_nuclear_event ev;
    double dir[3];
    double norm2;
    size_t j;
    int i;

    dir[0] = 0.0;
    dir[1] = 0.0;
    dir[2] = 1.0;

    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, 7u, 0u);
    for (i = 0; i < 1000; ++i) {
        osh_nuclear_abrasion_step(140.0, dir, O16_A, O16_Z, SIGMA_PA_CM2, &rng, &ev);
        for (j = 0u; j < ev.n_secondaries; ++j) {
            double *d = ev.secondaries[j].dir;
            norm2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
            ASSERT_NEAR(norm2, 1.0, 1.0e-10);
        }
    }
}

static void test_species_ratio(void) {
    /* Knock-outs are 50/50 n/p for O-16 (N/A = 0.5), but the escaping cascade
     * proton dilutes the neutron fraction: with mean knock-out count m and
     * escape probability close to 1, expect ~ 0.5·m / (m + 1). */
    struct osh_rng rng;
    struct osh_nuclear_event ev;
    double dir[3];
    int n_neutron;
    int n_total;
    double frac;
    size_t j;
    int i;

    dir[0] = 0.0;
    dir[1] = 0.0;
    dir[2] = 1.0;
    n_neutron = 0;
    n_total = 0;

    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, 13u, 0u);
    for (i = 0; i < 50000; ++i) {
        osh_nuclear_abrasion_step(140.0, dir, O16_A, O16_Z, SIGMA_PA_CM2, &rng, &ev);
        for (j = 0u; j < ev.n_secondaries; ++j) {
            ++n_total;
            if (ev.secondaries[j].species->pdg == 2112) { /* OSH_PART_PDG_NEUTRON */
                ++n_neutron;
            }
        }
    }
    ASSERT_TRUE(n_total > 0);
    frac = (double) n_neutron / (double) n_total;
    ASSERT_TRUE(frac > 0.20 && frac < 0.40);
}

static void test_mean_nu(void) {
    /* Knock-out count follows a zero-truncated Poisson with mean
     * λ/(1−e^{−λ}); the escaping cascade proton adds at most one more
     * secondary (slightly less due to absorption below ~1 MeV). */
    struct osh_rng rng;
    struct osh_nuclear_event ev;
    double dir[3];
    double total_nu;
    double lambda;
    double expected_ztp;
    double mean;
    int n_events;
    int i;

    dir[0] = 0.0;
    dir[1] = 0.0;
    dir[2] = 1.0;
    total_nu = 0.0;
    n_events = 20000;

    lambda = (OSH_ABRASION_SIGMA_PN_MB * OSH_MB_TO_CM2 * O16_A) / SIGMA_PA_CM2;
    expected_ztp = lambda / (1.0 - exp(-lambda));

    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, 42u, 0u);
    for (i = 0; i < n_events; ++i) {
        osh_nuclear_abrasion_step(140.0, dir, O16_A, O16_Z, SIGMA_PA_CM2, &rng, &ev);
        total_nu += (double) ev.n_secondaries;
    }
    mean = total_nu / n_events;
    ASSERT_TRUE(mean > expected_ztp * 0.9);
    ASSERT_TRUE(mean < expected_ztp + 1.0 + 1.0e-9);
}

static void test_excitation_energy(void) {
    /* E* accumulates the per-hole charge taken from the cascade proton (plus
     * an absorbed sub-MeV remnant): never negative, bounded by the per-hole
     * cost of the knock-outs unless the cascade proton was absorbed, and
     * always closing the exact energy budget (checked in
     * test_energy_conservation). */
    struct osh_rng rng;
    struct osh_nuclear_event ev;
    double dir[3];
    double e_star;
    int i;

    dir[0] = 0.0;
    dir[1] = 0.0;
    dir[2] = 1.0;

    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, 17u, 0u);
    for (i = 0; i < 2000; ++i) {
        osh_nuclear_abrasion_step(140.0, dir, O16_A, O16_Z, SIGMA_PA_CM2, &rng, &ev);
        if (ev.n_fragments == 0u) {
            continue;
        }
        e_star = ev.fragments[0].excitation_energy;
        ASSERT_TRUE(e_star >= 0.0);
        ASSERT_TRUE(e_star < 140.0);
        /* At least one knock-out happened (zero-truncated ν), so some hole
         * excitation must have been booked. */
        ASSERT_TRUE(ev.n_secondaries >= 1u);
        ASSERT_TRUE(e_star > 0.0);
    }
}

static void test_fragment_momentum_balance(void) {
    /* p_fragment = p_incident − Σ p_emitted (the absorbed cascade proton's
     * momentum is captured by the fragment). */
    struct osh_rng rng;
    struct osh_nuclear_event ev;
    double dir[3];
    double p_in;
    double p_expect[3];
    double p_sec;
    size_t j;
    int i;

    dir[0] = 0.6;
    dir[1] = 0.0;
    dir[2] = 0.8;
    p_in = sqrt(140.0 * (140.0 + 2.0 * OSH_PART_MASS_PROTON));

    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, 29u, 0u);
    for (i = 0; i < 2000; ++i) {
        osh_nuclear_abrasion_step(140.0, dir, O16_A, O16_Z, SIGMA_PA_CM2, &rng, &ev);
        if (ev.n_fragments == 0u) {
            continue;
        }
        p_expect[0] = p_in * dir[0];
        p_expect[1] = p_in * dir[1];
        p_expect[2] = p_in * dir[2];
        for (j = 0u; j < ev.n_secondaries; ++j) {
            p_sec = sqrt(ev.secondaries[j].energy * (ev.secondaries[j].energy + 2.0 * OSH_PART_MASS_PROTON));
            p_expect[0] -= p_sec * ev.secondaries[j].dir[0];
            p_expect[1] -= p_sec * ev.secondaries[j].dir[1];
            p_expect[2] -= p_sec * ev.secondaries[j].dir[2];
        }
        ASSERT_NEAR(ev.fragments[0].p[0], p_expect[0], 1.0e-9);
        ASSERT_NEAR(ev.fragments[0].p[1], p_expect[1], 1.0e-9);
        ASSERT_NEAR(ev.fragments[0].p[2], p_expect[2], 1.0e-9);
    }
}

static void test_exciton_bookkeeping(void) {
    /* Exciton export (issue #263).  Per event: every cascade collision
     * leaves a hole; a retained knock-out is a particle exciton and does NOT
     * reduce A; an absorbed cascade proton is one more particle exciton.
     * Externally observable invariants:
     *   escaped  = A_initial - A_fragment
     *   retained = excitons_h - escaped            (>= 0)
     *   excitons_p - retained ∈ {0, 1}             (1 iff cascade absorbed) */
    struct osh_rng rng;
    struct osh_nuclear_event ev;
    double dir[3];
    double total_p;
    unsigned int escaped;
    unsigned int retained;
    unsigned int slack;
    int i;

    dir[0] = 0.0;
    dir[1] = 0.0;
    dir[2] = 1.0;
    total_p = 0.0;

    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, 61u, 0u);
    for (i = 0; i < 20000; ++i) {
        osh_nuclear_abrasion_step(140.0, dir, O16_A, O16_Z, SIGMA_PA_CM2, &rng, &ev);
        if (ev.n_fragments == 0u) {
            continue;
        }
        ASSERT_TRUE(ev.fragments[0].a <= (unsigned int) O16_A);
        escaped = (unsigned int) O16_A - ev.fragments[0].a;
        ASSERT_TRUE(ev.fragments[0].excitons_h >= escaped);
        retained = ev.fragments[0].excitons_h - escaped;
        ASSERT_TRUE(ev.fragments[0].excitons_p >= retained);
        slack = ev.fragments[0].excitons_p - retained;
        ASSERT_TRUE(slack == 0u || slack == 1u);
        total_p += (double) ev.fragments[0].excitons_p;
    }
    /* With the default threshold a few % of 140 MeV knock-outs are retained:
     * the particle-exciton channel must actually fire. */
    ASSERT_TRUE(total_p > 0.0);
}

static void test_retention_low_energy(void) {
    /* At T = 3 MeV a single collision runs; the cascade proton can never pay
     * the 13.3 MeV hole cost, so it is always absorbed (one particle
     * exciton), and the knock-out (e_sec <= 3 MeV, well under the threshold)
     * is retained with high probability.  Retained events keep A/Z intact
     * and put the entire budget into E*. */
    struct osh_rng rng;
    struct osh_nuclear_event ev;
    double dir[3];
    double t_in;
    double mean_p;
    double total_p;
    int n_events;
    int i;

    dir[0] = 0.0;
    dir[1] = 0.0;
    dir[2] = 1.0;
    t_in = 3.0;
    total_p = 0.0;
    n_events = 2000;

    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, 71u, 0u);
    for (i = 0; i < n_events; ++i) {
        osh_nuclear_abrasion_step(t_in, dir, O16_A, O16_Z, SIGMA_PA_CM2, &rng, &ev);
        ASSERT_TRUE(ev.n_fragments == 1u);
        ASSERT_TRUE(ev.fragments[0].excitons_h == 1u);
        ASSERT_TRUE(ev.fragments[0].excitation_energy <= t_in + 1.0e-9);
        if (ev.n_secondaries == 0u) {
            /* Knock-out retained AND cascade proton absorbed: the full
             * kinetic budget must sit in E* and A/Z must be untouched. */
            ASSERT_NEAR(ev.fragments[0].excitation_energy, t_in, 1.0e-9);
            ASSERT_TRUE(ev.fragments[0].a == (unsigned int) O16_A);
            ASSERT_TRUE(ev.fragments[0].z == (unsigned int) O16_Z);
            ASSERT_TRUE(ev.fragments[0].excitons_p == 2u);
        }
        total_p += (double) ev.fragments[0].excitons_p;
    }
    /* Retention probability at ~1.5 MeV is >0.9, so the mean particle count
     * sits close to 2 (retained knock-out + absorbed cascade proton). */
    mean_p = total_p / (double) n_events;
    ASSERT_TRUE(mean_p > 1.7);
    ASSERT_TRUE(mean_p <= 2.0);
}

int main(void) {
    test_event_kind();
    test_energy_conservation();
    test_secondary_dir_unit();
    test_species_ratio();
    test_mean_nu();
    test_excitation_energy();
    test_fragment_momentum_balance();
    test_exciton_bookkeeping();
    test_retention_low_energy();
    return 0;
}
