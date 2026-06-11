#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "openshieldhit/const.h"
#include "particle/osh_particle.h"
#include "particle/osh_particle_const.h"
#include "particle/osh_particle_pdg.h"
#include "physics/nuclear/osh_nuclear_fermi_breakup.h"
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

/* Map a final-product species descriptor back to (z, a).  The registry
 * descriptors for n and p carry z = a = 0 (non-ion branch), so they need an
 * explicit mapping. */
static void species_za(struct particle const *sp, unsigned int *z, unsigned int *a) {
    if (sp->pdg == OSH_PART_PDG_NEUTRON) {
        *z = 0u;
        *a = 1u;
        return;
    }
    if (sp->pdg == OSH_PART_PDG_PROTON) {
        *z = 1u;
        *a = 1u;
        return;
    }
    *z = sp->z;
    *a = sp->a;
}

/* Ground-state nuclear mass [MeV/c²] from the same sources the model uses. */
static double za_mass(unsigned int z, unsigned int a) {
    double m;

    if (z == 0u && a == 1u) {
        return OSH_PART_MASS_NEUTRON;
    }
    ASSERT_TRUE(osh_particle_nuclear_mass_mev_from_za(z, a, &m));
    return m;
}

/* Set up an event as the abrasion stage leaves it: one prefragment, no
 * secondaries. */
static void
init_abrasion_event(struct osh_nuclear_event *ev, unsigned int z, unsigned int a, double e_star, double const p[3]) {
    memset(ev, 0, sizeof(*ev));
    ev->kind = OSH_NUCLEAR_EVENT_ABRASION;
    ev->n_fragments = 1u;
    ev->fragments[0].z = z;
    ev->fragments[0].a = a;
    ev->fragments[0].excitation_energy = e_star;
    ev->fragments[0].p[0] = p[0];
    ev->fragments[0].p[1] = p[1];
    ev->fragments[0].p[2] = p[2];
}

static int is_whitelisted_pdg(int pdg) {
    return pdg == OSH_PART_PDG_NEUTRON || pdg == OSH_PART_PDG_PROTON || pdg == OSH_PART_PDG_DEUTERON
           || pdg == OSH_PART_PDG_TRITON || pdg == OSH_PART_PDG_HE3 || pdg == OSH_PART_PDG_HE4;
}

static void test_compile_channels(void) {
    struct osh_nuclear_fermi_breakup model;
    size_t idx;
    uint16_t off;
    uint16_t cnt;
    uint16_t c;
    int found_2alpha;

    memset(&model, 0, sizeof(model));
    ASSERT_TRUE(osh_nuclear_fermi_breakup_compile(&model) == OSH_OK);
    ASSERT_TRUE(model.nchannels > 0u);

    /* Be-8 must have an open 2-alpha channel at E* = 0 with Q ≈ +0.092 MeV. */
    idx = 4u * (OSH_FERMI_BREAKUP_AMAX + 1u) + 8u;
    off = model.chan_offset[idx];
    cnt = model.chan_count[idx];
    ASSERT_TRUE(cnt > 0u);
    found_2alpha = 0;
    for (c = 0u; c < cnt; ++c) {
        struct osh_fermi_channel const *ch = &model.channel_pool[off + c];
        if (ch->z1 == 2u && ch->a1 == 4u && ch->z2 == 2u && ch->a2 == 4u) {
            found_2alpha = 1;
            ASSERT_TRUE(ch->q_mev > 0.05f && ch->q_mev < 0.15f);
        }
    }
    ASSERT_TRUE(found_2alpha);

    osh_nuclear_fermi_breakup_free(&model);
}

static void test_be8_breaks_to_two_alphas(void) {
    struct osh_nuclear_fermi_breakup model;
    struct osh_rng rng;
    struct osh_nuclear_event ev;
    double p0[3];
    double q_expected;
    double sum_ke;
    double sum_p[3];
    double pj;
    size_t j;
    int i;

    memset(&model, 0, sizeof(model));
    ASSERT_TRUE(osh_nuclear_fermi_breakup_compile(&model) == OSH_OK);
    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, 11u, 0u);

    p0[0] = 0.0;
    p0[1] = 0.0;
    p0[2] = 0.0;
    q_expected = za_mass(4u, 8u) - 2.0 * za_mass(2u, 4u);

    for (i = 0; i < 200; ++i) {
        init_abrasion_event(&ev, 4u, 8u, 0.0, p0);
        osh_nuclear_fermi_breakup_step(&model, &ev.fragments[0], &rng, &ev);

        ASSERT_TRUE(ev.kind == OSH_NUCLEAR_EVENT_FRAGMENTATION);
        ASSERT_TRUE(ev.n_secondaries == 2u);
        ASSERT_TRUE(ev.n_fragments == 0u);
        ASSERT_TRUE(ev.secondaries[0].species->pdg == OSH_PART_PDG_HE4);
        ASSERT_TRUE(ev.secondaries[1].species->pdg == OSH_PART_PDG_HE4);

        sum_ke = 0.0;
        sum_p[0] = 0.0;
        sum_p[1] = 0.0;
        sum_p[2] = 0.0;
        for (j = 0u; j < ev.n_secondaries; ++j) {
            double m = za_mass(2u, 4u);
            sum_ke += ev.secondaries[j].energy;
            pj = sqrt(ev.secondaries[j].energy * (ev.secondaries[j].energy + 2.0 * m));
            sum_p[0] += pj * ev.secondaries[j].dir[0];
            sum_p[1] += pj * ev.secondaries[j].dir[1];
            sum_p[2] += pj * ev.secondaries[j].dir[2];
        }
        ASSERT_NEAR(sum_ke, q_expected, 1.0e-3);
        ASSERT_NEAR(sum_p[0], 0.0, 1.0e-6);
        ASSERT_NEAR(sum_p[1], 0.0, 1.0e-6);
        ASSERT_NEAR(sum_p[2], 0.0, 1.0e-6);
    }

    osh_nuclear_fermi_breakup_free(&model);
}

static void test_o16_below_threshold_untouched(void) {
    struct osh_nuclear_fermi_breakup model;
    struct osh_rng rng;
    struct osh_nuclear_event ev;
    double p0[3];

    memset(&model, 0, sizeof(model));
    ASSERT_TRUE(osh_nuclear_fermi_breakup_compile(&model) == OSH_OK);
    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, 23u, 0u);

    /* Lowest O-16 particle-separation threshold is the alpha channel at
     * ≈ 7.16 MeV; E* = 1 MeV must leave the event untouched. */
    p0[0] = 0.0;
    p0[1] = 0.0;
    p0[2] = 0.0;
    init_abrasion_event(&ev, 8u, 16u, 1.0, p0);
    osh_nuclear_fermi_breakup_step(&model, &ev.fragments[0], &rng, &ev);

    ASSERT_TRUE(ev.kind == OSH_NUCLEAR_EVENT_ABRASION);
    ASSERT_TRUE(ev.n_secondaries == 0u);
    ASSERT_TRUE(ev.n_fragments == 1u);
    ASSERT_TRUE(ev.fragments[0].z == 8u && ev.fragments[0].a == 16u);
    ASSERT_NEAR(ev.fragments[0].excitation_energy, 1.0, 1.0e-12);

    osh_nuclear_fermi_breakup_free(&model);
}

static void test_conservation_c12(void) {
    static double const e_stars[3] = {10.0, 30.0, 100.0};
    struct osh_nuclear_fermi_breakup model;
    struct osh_rng rng;
    struct osh_nuclear_event ev;
    double p0[3];
    double e_total_expected;
    double e_total;
    double sum_p[3];
    double m;
    double pj;
    unsigned int sum_z;
    unsigned int sum_a;
    unsigned int z;
    unsigned int a;
    size_t j;
    size_t k;
    int i;

    memset(&model, 0, sizeof(model));
    ASSERT_TRUE(osh_nuclear_fermi_breakup_compile(&model) == OSH_OK);
    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, 37u, 0u);

    p0[0] = 0.0;
    p0[1] = 0.0;
    p0[2] = 0.0;

    for (k = 0u; k < 3u; ++k) {
        e_total_expected = za_mass(6u, 12u) + e_stars[k];
        for (i = 0; i < 2000; ++i) {
            init_abrasion_event(&ev, 6u, 12u, e_stars[k], p0);
            osh_nuclear_fermi_breakup_step(&model, &ev.fragments[0], &rng, &ev);

            /* All listed E* are above the lowest C-12 threshold (3-alpha via
             * alpha + Be-8 at ≈ 7.27 MeV), so the break-up must fire. */
            ASSERT_TRUE(ev.n_secondaries > 0u || ev.n_fragments > 1u);

            sum_z = 0u;
            sum_a = 0u;
            e_total = 0.0;
            sum_p[0] = 0.0;
            sum_p[1] = 0.0;
            sum_p[2] = 0.0;

            for (j = 0u; j < ev.n_secondaries; ++j) {
                species_za(ev.secondaries[j].species, &z, &a);
                sum_z += z;
                sum_a += a;
                m = za_mass(z, a);
                e_total += ev.secondaries[j].energy + m;
                pj = sqrt(ev.secondaries[j].energy * (ev.secondaries[j].energy + 2.0 * m));
                sum_p[0] += pj * ev.secondaries[j].dir[0];
                sum_p[1] += pj * ev.secondaries[j].dir[1];
                sum_p[2] += pj * ev.secondaries[j].dir[2];
            }
            for (j = 0u; j < ev.n_fragments; ++j) {
                struct osh_nuclear_fragment const *f = &ev.fragments[j];
                sum_z += f->z;
                sum_a += f->a;
                m = za_mass(f->z, f->a) + f->excitation_energy;
                e_total += sqrt(m * m + f->p[0] * f->p[0] + f->p[1] * f->p[1] + f->p[2] * f->p[2]);
                sum_p[0] += f->p[0];
                sum_p[1] += f->p[1];
                sum_p[2] += f->p[2];
            }

            ASSERT_TRUE(sum_z == 6u);
            ASSERT_TRUE(sum_a == 12u);
            ASSERT_NEAR(e_total, e_total_expected, 1.0e-2);
            ASSERT_NEAR(sum_p[0], 0.0, 1.0e-3);
            ASSERT_NEAR(sum_p[1], 0.0, 1.0e-3);
            ASSERT_NEAR(sum_p[2], 0.0, 1.0e-3);
        }
    }

    osh_nuclear_fermi_breakup_free(&model);
}

static void test_whitelist_only(void) {
    struct osh_nuclear_fermi_breakup model;
    struct osh_rng rng;
    struct osh_nuclear_event ev;
    double p0[3];
    size_t j;
    int i;

    memset(&model, 0, sizeof(model));
    ASSERT_TRUE(osh_nuclear_fermi_breakup_compile(&model) == OSH_OK);
    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, 71u, 0u);

    p0[0] = 0.0;
    p0[1] = 0.0;
    p0[2] = 0.0;

    for (i = 0; i < 5000; ++i) {
        init_abrasion_event(&ev, 8u, 16u, (i % 2 == 0) ? 50.0 : 120.0, p0);
        osh_nuclear_fermi_breakup_step(&model, &ev.fragments[0], &rng, &ev);
        for (j = 0u; j < ev.n_secondaries; ++j) {
            ASSERT_TRUE(is_whitelisted_pdg(ev.secondaries[j].species->pdg));
        }
    }

    osh_nuclear_fermi_breakup_free(&model);
}

static void test_c12_three_alpha_channel(void) {
    struct osh_nuclear_fermi_breakup model;
    struct osh_rng rng;
    struct osh_nuclear_event ev;
    double p0[3];
    size_t j;
    int i;
    int n_three_alpha;

    memset(&model, 0, sizeof(model));
    ASSERT_TRUE(osh_nuclear_fermi_breakup_compile(&model) == OSH_OK);
    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, 5u, 0u);

    p0[0] = 0.0;
    p0[1] = 0.0;
    p0[2] = 0.0;
    n_three_alpha = 0;

    for (i = 0; i < 5000; ++i) {
        int n_alpha;

        init_abrasion_event(&ev, 6u, 12u, 15.0, p0);
        osh_nuclear_fermi_breakup_step(&model, &ev.fragments[0], &rng, &ev);
        n_alpha = 0;
        for (j = 0u; j < ev.n_secondaries; ++j) {
            if (ev.secondaries[j].species->pdg == OSH_PART_PDG_HE4) {
                ++n_alpha;
            }
        }
        if (n_alpha == 3 && ev.n_secondaries == 3u) {
            ++n_three_alpha;
        }
    }
    /* The 3-alpha final state must occur with non-negligible frequency just
     * above its threshold. */
    ASSERT_TRUE(n_three_alpha > 50);

    osh_nuclear_fermi_breakup_free(&model);
}

static void test_capacity_truncation(void) {
    struct osh_nuclear_fermi_breakup model;
    struct osh_rng rng;
    struct osh_nuclear_event ev;
    double p0[3];
    size_t prefill;
    unsigned int sum_a;
    unsigned int z;
    unsigned int a;
    size_t j;
    int i;

    memset(&model, 0, sizeof(model));
    ASSERT_TRUE(osh_nuclear_fermi_breakup_compile(&model) == OSH_OK);
    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, 91u, 0u);

    p0[0] = 0.0;
    p0[1] = 0.0;
    p0[2] = 0.0;
    prefill = (size_t) OSH_NUCLEAR_MAX_SECONDARIES - 2u;

    for (i = 0; i < 500; ++i) {
        init_abrasion_event(&ev, 8u, 16u, 300.0, p0);
        /* Simulate an abrasion stage that already nearly filled the array. */
        ev.n_secondaries = prefill;
        osh_nuclear_fermi_breakup_step(&model, &ev.fragments[0], &rng, &ev);

        ASSERT_TRUE(ev.n_secondaries <= (size_t) OSH_NUCLEAR_MAX_SECONDARIES);
        ASSERT_TRUE(ev.n_fragments <= (size_t) OSH_NUCLEAR_MAX_FRAGMENTS);

        /* Whatever was emitted plus the unprocessed residues must never
         * exceed the parent nucleon count. */
        sum_a = 0u;
        for (j = prefill; j < ev.n_secondaries; ++j) {
            species_za(ev.secondaries[j].species, &z, &a);
            sum_a += a;
        }
        for (j = 0u; j < ev.n_fragments; ++j) {
            sum_a += ev.fragments[j].a;
        }
        ASSERT_TRUE(sum_a <= 16u);
    }

    osh_nuclear_fermi_breakup_free(&model);
}

static void test_moving_parent(void) {
    struct osh_nuclear_fermi_breakup model;
    struct osh_rng rng;
    struct osh_nuclear_event ev;
    double p0[3];
    double m_parent;
    double e_total_expected;
    double e_total;
    double sum_p[3];
    double m;
    double pj;
    size_t j;
    int i;

    memset(&model, 0, sizeof(model));
    ASSERT_TRUE(osh_nuclear_fermi_breakup_compile(&model) == OSH_OK);
    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, 101u, 0u);

    /* Be-8 in flight along +z: products must conserve the parent 4-momentum
     * and both alphas must be forward-boosted (beta_parent >> beta_cm). */
    p0[0] = 0.0;
    p0[1] = 0.0;
    p0[2] = 500.0;
    m_parent = za_mass(4u, 8u);
    e_total_expected = sqrt(m_parent * m_parent + p0[2] * p0[2]);

    for (i = 0; i < 500; ++i) {
        init_abrasion_event(&ev, 4u, 8u, 0.0, p0);
        osh_nuclear_fermi_breakup_step(&model, &ev.fragments[0], &rng, &ev);

        ASSERT_TRUE(ev.kind == OSH_NUCLEAR_EVENT_FRAGMENTATION);
        ASSERT_TRUE(ev.n_secondaries == 2u);

        e_total = 0.0;
        sum_p[0] = 0.0;
        sum_p[1] = 0.0;
        sum_p[2] = 0.0;
        for (j = 0u; j < ev.n_secondaries; ++j) {
            m = za_mass(2u, 4u);
            e_total += ev.secondaries[j].energy + m;
            pj = sqrt(ev.secondaries[j].energy * (ev.secondaries[j].energy + 2.0 * m));
            sum_p[0] += pj * ev.secondaries[j].dir[0];
            sum_p[1] += pj * ev.secondaries[j].dir[1];
            sum_p[2] += pj * ev.secondaries[j].dir[2];
            ASSERT_TRUE(pj * ev.secondaries[j].dir[2] > 0.0);
        }
        ASSERT_NEAR(e_total, e_total_expected, 1.0e-6);
        ASSERT_NEAR(sum_p[0], 0.0, 1.0e-6);
        ASSERT_NEAR(sum_p[1], 0.0, 1.0e-6);
        ASSERT_NEAR(sum_p[2], 500.0, 1.0e-6);
    }

    osh_nuclear_fermi_breakup_free(&model);
}

static void test_g4_multiplicity_anchors(void) {
    /* Anchor points from the canonical G4FermiBreakUp (Geant4 9.1 fixed),
     * extracted from I. Pshenichnov's standalone FermiTest (used with
     * permission; see examples/05_fermi_breakup_validation/).  Both lie in
     * the low-excitation region populated by the abrasion stage, where the
     * sequential-binary approximation reproduces the canonical model:
     *   C-12 at E* = 12.6 MeV (1.05 MeV/nucleon): <mult> = 3.015 (3-alpha)
     *   C-12 at E* = 22.2 MeV (1.85 MeV/nucleon): <mult> = 2.356
     * The high-excitation region (above ~3 MeV/nucleon) is a known deficit
     * of the binary scheme and is deliberately not pinned here. */
    struct osh_nuclear_fermi_breakup model;
    struct osh_rng rng;
    struct osh_nuclear_event ev;
    double p0[3];
    double mult_sum;
    double mean;
    int i;
    int n_events;

    memset(&model, 0, sizeof(model));
    ASSERT_TRUE(osh_nuclear_fermi_breakup_compile(&model) == OSH_OK);
    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, 314u, 0u);

    p0[0] = 0.0;
    p0[1] = 0.0;
    p0[2] = 0.0;
    n_events = 4000;

    mult_sum = 0.0;
    for (i = 0; i < n_events; ++i) {
        init_abrasion_event(&ev, 6u, 12u, 12.6, p0);
        osh_nuclear_fermi_breakup_step(&model, &ev.fragments[0], &rng, &ev);
        mult_sum += (double) (ev.n_secondaries + ev.n_fragments);
    }
    mean = mult_sum / n_events;
    ASSERT_NEAR(mean, 3.015, 3.015 * 0.05);

    mult_sum = 0.0;
    for (i = 0; i < n_events; ++i) {
        init_abrasion_event(&ev, 6u, 12u, 22.2, p0);
        osh_nuclear_fermi_breakup_step(&model, &ev.fragments[0], &rng, &ev);
        mult_sum += (double) (ev.n_secondaries + ev.n_fragments);
    }
    mean = mult_sum / n_events;
    ASSERT_NEAR(mean, 2.356, 2.356 * 0.05);

    osh_nuclear_fermi_breakup_free(&model);
}

int main(void) {
    test_compile_channels();
    test_be8_breaks_to_two_alphas();
    test_o16_below_threshold_untouched();
    test_conservation_c12();
    test_whitelist_only();
    test_c12_three_alpha_channel();
    test_capacity_truncation();
    test_moving_parent();
    test_g4_multiplicity_anchors();
    return 0;
}
