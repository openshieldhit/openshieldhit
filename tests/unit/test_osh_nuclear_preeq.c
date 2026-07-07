#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "openshieldhit/const.h"
#include "particle/osh_particle.h"
#include "particle/osh_particle_pdg.h"
#include "physics/nuclear/osh_nuclear_handler.h"
#include "physics/nuclear/osh_nuclear_preeq.h"
#include "random/osh_rng.h"

#define ASSERT_TRUE(cond)                                                                                              \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            fprintf(stderr, "ASSERT FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__);                                 \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

#define ASSERT_NEAR(a, b, tol) ASSERT_TRUE(fabs((a) - (b)) <= (tol))

/* Whitelist check: (z, a) of every emittable species. */
static int is_whitelisted(struct particle const *sp) {
    if (sp->pdg == OSH_PART_PDG_NEUTRON || sp->pdg == OSH_PART_PDG_PROTON) {
        return 1;
    }
    if (sp->z == 1u && (sp->a == 2u || sp->a == 3u)) {
        return 1;
    }
    if (sp->z == 2u && (sp->a == 3u || sp->a == 4u)) {
        return 1;
    }
    return 0;
}

/* Build a bare excited O-16-like fragment with a given exciton config. */
static void make_fragment(
    struct osh_nuclear_fragment *f, unsigned int z, unsigned int a, double e_star, unsigned int np, unsigned int nh) {
    memset(f, 0, sizeof(*f));
    f->z = z;
    f->a = a;
    f->excitation_energy = e_star;
    f->excitons_p = np;
    f->excitons_h = nh;
}

/* Compare the physically meaningful fragment state field-by-field.
 * Padding bytes are not stable, so memcmp() is not valid here. */
static void assert_fragment_equal(struct osh_nuclear_fragment const *lhs, struct osh_nuclear_fragment const *rhs) {
    ASSERT_NEAR(lhs->excitation_energy, rhs->excitation_energy, 1.0e-12);
    ASSERT_NEAR(lhs->p[0], rhs->p[0], 1.0e-12);
    ASSERT_NEAR(lhs->p[1], rhs->p[1], 1.0e-12);
    ASSERT_NEAR(lhs->p[2], rhs->p[2], 1.0e-12);
    ASSERT_TRUE(lhs->z == rhs->z);
    ASSERT_TRUE(lhs->a == rhs->a);
    ASSERT_TRUE(lhs->excitons_p == rhs->excitons_p);
    ASSERT_TRUE(lhs->excitons_h == rhs->excitons_h);
}

static void test_thermalized_noop(void) {
    /* (p, h) = (0, 0) is the thermalized contract: compound nuclei and
     * break-up residues must pass through untouched. */
    struct osh_nuclear_preeq model;
    struct osh_rng rng;
    struct osh_nuclear_event ev;
    struct osh_nuclear_fragment f;
    struct osh_nuclear_fragment before;
    int i;

    ASSERT_TRUE(osh_nuclear_preeq_compile(&model) == OSH_OK);
    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, 11u, 0u);

    for (i = 0; i < 100; ++i) {
        memset(&ev, 0, sizeof(ev));
        make_fragment(&f, 8u, 16u, 40.0, 0u, 0u);
        before = f;
        osh_nuclear_preeq_step(&model, &f, &rng, &ev);
        assert_fragment_equal(&before, &f);
        ASSERT_TRUE(ev.n_secondaries == 0u);
    }
}

static void test_below_threshold_noop(void) {
    /* Excitation below every separation energy: nothing can be emitted;
     * the fragment may only thermalize (excitons change, E* and A/Z not). */
    struct osh_nuclear_preeq model;
    struct osh_rng rng;
    struct osh_nuclear_event ev;
    struct osh_nuclear_fragment f;
    int i;

    ASSERT_TRUE(osh_nuclear_preeq_compile(&model) == OSH_OK);
    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, 13u, 0u);

    for (i = 0; i < 200; ++i) {
        memset(&ev, 0, sizeof(ev));
        make_fragment(&f, 8u, 16u, 3.0, 1u, 1u);
        osh_nuclear_preeq_step(&model, &f, &rng, &ev);
        ASSERT_TRUE(ev.n_secondaries == 0u);
        ASSERT_TRUE(f.z == 8u && f.a == 16u);
        ASSERT_NEAR(f.excitation_energy, 3.0, 1.0e-12);
    }
}

static void test_conservation(void) {
    /* Exact per-event bookkeeping with masses:
     *   M(parent) + E*_before == sum(m_j + T_j) + M(residue) + E*_after
     * and fragment momentum == -(sum of emitted momenta) for a parent at
     * rest.  Uses the same isotope masses the model compiled. */
    struct osh_nuclear_preeq model;
    struct osh_rng rng;
    struct osh_nuclear_event ev;
    struct osh_nuclear_fragment f;
    double m_parent;
    double m_res;
    double m_j;
    double lhs;
    double rhs;
    double p_expect[3];
    double p_j;
    size_t j;
    int i;

    ASSERT_TRUE(osh_nuclear_preeq_compile(&model) == OSH_OK);
    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, 17u, 0u);
    ASSERT_TRUE(osh_particle_nuclear_mass_mev_from_za(8u, 16u, &m_parent));

    for (i = 0; i < 5000; ++i) {
        memset(&ev, 0, sizeof(ev));
        make_fragment(&f, 8u, 16u, 45.0, 2u, 2u);
        osh_nuclear_preeq_step(&model, &f, &rng, &ev);

        ASSERT_TRUE(osh_particle_nuclear_mass_mev_from_za(f.z, f.a, &m_res));
        lhs = m_parent + 45.0;
        rhs = m_res + f.excitation_energy;
        p_expect[0] = 0.0;
        p_expect[1] = 0.0;
        p_expect[2] = 0.0;
        for (j = 0u; j < ev.n_secondaries; ++j) {
            m_j = ev.secondaries[j].species->mass;
            rhs += m_j + ev.secondaries[j].energy;
            p_j = sqrt(ev.secondaries[j].energy * (ev.secondaries[j].energy + 2.0 * m_j));
            p_expect[0] -= p_j * ev.secondaries[j].dir[0];
            p_expect[1] -= p_j * ev.secondaries[j].dir[1];
            p_expect[2] -= p_j * ev.secondaries[j].dir[2];
        }
        ASSERT_NEAR(lhs, rhs, 1.0e-6);
        ASSERT_NEAR(f.p[0], p_expect[0], 1.0e-9);
        ASSERT_NEAR(f.p[1], p_expect[1], 1.0e-9);
        ASSERT_NEAR(f.p[2], p_expect[2], 1.0e-9);
        ASSERT_TRUE(f.excitation_energy >= 0.0);
    }
}

static void test_whitelist_and_kind(void) {
    /* Only {n, p, d, t, He-3, He-4} may be emitted; FRAGMENTATION is set
     * exactly when something was emitted. */
    struct osh_nuclear_preeq model;
    struct osh_rng rng;
    struct osh_nuclear_event ev;
    struct osh_nuclear_fragment f;
    size_t j;
    int i;

    ASSERT_TRUE(osh_nuclear_preeq_compile(&model) == OSH_OK);
    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, 19u, 0u);

    for (i = 0; i < 5000; ++i) {
        memset(&ev, 0, sizeof(ev));
        make_fragment(&f, 8u, 16u, 50.0, 3u, 2u);
        osh_nuclear_preeq_step(&model, &f, &rng, &ev);
        for (j = 0u; j < ev.n_secondaries; ++j) {
            ASSERT_TRUE(is_whitelisted(ev.secondaries[j].species));
            ASSERT_TRUE(ev.secondaries[j].energy > 0.0);
        }
        if (ev.n_secondaries > 0u) {
            ASSERT_TRUE(ev.kind == OSH_NUCLEAR_EVENT_FRAGMENTATION);
        }
    }
}

static void test_emission_statistics(void) {
    /* Physics smoke test at the #212-relevant configuration: O-16 with
     * E* = 45 MeV and a (2, 2) exciton config must emit fast ejectiles in a
     * sizeable fraction of events, including clusters, with mean nucleon
     * energy well above the FBU evaporation-like scale, and must always
     * lower the excitation it hands to the equilibrium stage. */
    struct osh_nuclear_preeq model;
    struct osh_rng rng;
    struct osh_nuclear_event ev;
    struct osh_nuclear_fragment f;
    double n_events;
    double n_emitting;
    double n_cluster;
    double e_nucleon_sum;
    double n_nucleon;
    size_t j;
    int i;

    ASSERT_TRUE(osh_nuclear_preeq_compile(&model) == OSH_OK);
    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, 23u, 0u);

    n_events = 20000.0;
    n_emitting = 0.0;
    n_cluster = 0.0;
    e_nucleon_sum = 0.0;
    n_nucleon = 0.0;

    for (i = 0; i < (int) n_events; ++i) {
        memset(&ev, 0, sizeof(ev));
        make_fragment(&f, 8u, 16u, 45.0, 2u, 2u);
        osh_nuclear_preeq_step(&model, &f, &rng, &ev);
        ASSERT_TRUE(f.excitation_energy <= 45.0 + 1.0e-9);
        if (ev.n_secondaries > 0u) {
            n_emitting += 1.0;
        }
        for (j = 0u; j < ev.n_secondaries; ++j) {
            if (ev.secondaries[j].species->a >= 2u) {
                n_cluster += 1.0;
            } else {
                e_nucleon_sum += ev.secondaries[j].energy;
                n_nucleon += 1.0;
            }
        }
    }

    /* Loose statistical bands: the exact fractions are calibration targets
     * (#225 gamma scales), but the channels must exist and be fast. */
    ASSERT_TRUE(n_emitting / n_events > 0.10);
    ASSERT_TRUE(n_cluster > 0.0);
    ASSERT_TRUE(n_nucleon > 0.0);
    ASSERT_TRUE(e_nucleon_sum / n_nucleon > 5.0);
}

int main(void) {
    test_thermalized_noop();
    test_below_threshold_noop();
    test_conservation();
    test_whitelist_and_kind();
    test_emission_statistics();
    return 0;
}
