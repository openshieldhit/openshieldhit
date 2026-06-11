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
    /* For O-16, N/A = 0.5 → expect ~50% neutron secondaries. */
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
    if (n_total > 0) {
        frac = (double) n_neutron / (double) n_total;
        ASSERT_NEAR(frac, 0.5, 0.05);
    }
}

static void test_mean_nu(void) {
    /* Mean ν should match (sigma_pN * A / sigma_pA) within 10% over many events. */
    struct osh_rng rng;
    struct osh_nuclear_event ev;
    double dir[3];
    double total_nu;
    double expected_nu;
    int n_events;
    int i;

    dir[0] = 0.0;
    dir[1] = 0.0;
    dir[2] = 1.0;
    total_nu = 0.0;
    n_events = 20000;

    expected_nu = (OSH_ABRASION_SIGMA_PN_MB * OSH_MB_TO_CM2 * O16_A) / SIGMA_PA_CM2;

    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, 42u, 0u);
    for (i = 0; i < n_events; ++i) {
        osh_nuclear_abrasion_step(140.0, dir, O16_A, O16_Z, SIGMA_PA_CM2, &rng, &ev);
        total_nu += (double) ev.n_secondaries;
    }
    ASSERT_NEAR(total_nu / n_events, expected_nu, expected_nu * 0.10);
}

static void test_excitation_energy(void) {
    /* E* = per-hole constant × abraded nucleons, clamped to the leftover
     * cascade energy: never negative, never above the per-hole estimate, and
     * exactly 0 when no nucleon was knocked out. */
    struct osh_rng rng;
    struct osh_nuclear_event ev;
    double dir[3];
    double e_star;
    double sum_e;
    size_t j;
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
        ASSERT_TRUE(e_star <= OSH_ABRASION_EXCITATION_PER_HOLE_MEV * (double) ev.n_secondaries + 1.0e-9);
        if (ev.n_secondaries == 0u) {
            ASSERT_TRUE(e_star == 0.0);
        }

        /* Four-vector budget at the abrasion → break-up boundary: the emitted
         * kinetic energy plus the excitation may never exceed the incident
         * kinetic energy (no double counting of the leftover cascade energy). */
        sum_e = 0.0;
        for (j = 0u; j < ev.n_secondaries; ++j) {
            sum_e += ev.secondaries[j].energy;
        }
        ASSERT_TRUE(sum_e + e_star <= 140.0 + 1.0e-6);
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

int main(void) {
    test_event_kind();
    test_energy_conservation();
    test_secondary_dir_unit();
    test_species_ratio();
    test_mean_nu();
    test_excitation_energy();
    test_fragment_momentum_balance();
    return 0;
}
