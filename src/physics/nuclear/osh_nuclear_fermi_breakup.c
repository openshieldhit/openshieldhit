#include "physics/nuclear/osh_nuclear_fermi_breakup.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "particle/osh_particle.h"
#include "particle/osh_particle_const.h"
#include "particle/osh_particle_pdg.h"
#include "physics/nuclear/osh_nuclear_handler.h"
#include "physics/osh_kinematics.h"
#include "random/osh_rng.h"

/* Dense (z, a) table size: z in [0, ZMAX], a in [0, AMAX]. */
#define FBU_NDENSE ((size_t) (OSH_FERMI_BREAKUP_ZMAX + 1) * (size_t) (OSH_FERMI_BREAKUP_AMAX + 1))

/* Final-product species indices into model->species. */
#define FBU_SPECIES_NEUTRON 0
#define FBU_SPECIES_PROTON 1
#define FBU_SPECIES_DEUTERON 2
#define FBU_SPECIES_TRITON 3
#define FBU_SPECIES_HE3 4
#define FBU_SPECIES_HE4 5
#define FBU_SPECIES_COUNT 6

/** One pending nuclide on the de-excitation work stack. */
struct fbu_work_item {
    double p[3];   /* lab momentum [MeV/c] */
    double e_star; /* excitation energy [MeV]; 0 for split products */
    uint8_t z;
    uint8_t a;
};

/** Ground-state spin degeneracy entry g = 2J+1 (ENSDF/TUNL adopted spins). */
struct fbu_spin_entry {
    uint8_t z;
    uint8_t a;
    uint8_t g;
};

/* Nuclides not listed default to g = 1 (conservative; refine as needed). */
static struct fbu_spin_entry const s_spin_table[] = {
    {0u, 1u, 2u},  /* n     1/2  */
    {1u, 1u, 2u},  /* p     1/2  */
    {1u, 2u, 3u},  /* d     1    */
    {1u, 3u, 2u},  /* t     1/2  */
    {2u, 3u, 2u},  /* He-3  1/2  */
    {2u, 4u, 1u},  /* He-4  0    */
    {2u, 5u, 4u},  /* He-5  3/2  */
    {2u, 6u, 1u},  /* He-6  0    */
    {3u, 5u, 4u},  /* Li-5  3/2  */
    {3u, 6u, 3u},  /* Li-6  1    */
    {3u, 7u, 4u},  /* Li-7  3/2  */
    {3u, 8u, 5u},  /* Li-8  2    */
    {3u, 9u, 4u},  /* Li-9  3/2  */
    {4u, 7u, 4u},  /* Be-7  3/2  */
    {4u, 8u, 1u},  /* Be-8  0    */
    {4u, 9u, 4u},  /* Be-9  3/2  */
    {4u, 10u, 1u}, /* Be-10 0    */
    {5u, 8u, 5u},  /* B-8   2    */
    {5u, 9u, 4u},  /* B-9   3/2  */
    {5u, 10u, 7u}, /* B-10  3    */
    {5u, 11u, 4u}, /* B-11  3/2  */
    {5u, 12u, 3u}, /* B-12  1    */
    {6u, 10u, 1u}, /* C-10  0    */
    {6u, 11u, 4u}, /* C-11  3/2  */
    {6u, 12u, 1u}, /* C-12  0    */
    {6u, 13u, 2u}, /* C-13  1/2  */
    {6u, 14u, 1u}, /* C-14  0    */
    {7u, 12u, 3u}, /* N-12  1    */
    {7u, 13u, 2u}, /* N-13  1/2  */
    {7u, 14u, 3u}, /* N-14  1    */
    {7u, 15u, 2u}, /* N-15  1/2  */
    {8u, 14u, 1u}, /* O-14  0    */
    {8u, 15u, 2u}, /* O-15  1/2  */
    {8u, 16u, 1u}, /* O-16  0    */
};

/** Dense (z, a) index into the per-nuclide tables. */
static inline size_t _za_idx(unsigned int z, unsigned int a) {
    return (size_t) z * (size_t) (OSH_FERMI_BREAKUP_AMAX + 1) + (size_t) a;
}

/**
 * Ground-state spin degeneracy 2J+1 of nuclide (z, a); 1 if not tabulated.
 * Cold path only (channel-table compilation).
 */
static unsigned int spin_degeneracy(unsigned int z, unsigned int a) {
    size_t i;

    for (i = 0u; i < sizeof(s_spin_table) / sizeof(s_spin_table[0]); ++i) {
        if (s_spin_table[i].z == z && s_spin_table[i].a == a) {
            return s_spin_table[i].g;
        }
    }
    return 1u;
}

/**
 * Index into model->species if (z, a) is a whitelisted transportable product
 * (n, p, d, t, He-3, He-4); -1 otherwise.
 */
static int final_species_index(unsigned int z, unsigned int a) {
    if (a == 1u && z == 0u) {
        return FBU_SPECIES_NEUTRON;
    }
    if (z == 1u) {
        if (a == 1u) {
            return FBU_SPECIES_PROTON;
        }
        if (a == 2u) {
            return FBU_SPECIES_DEUTERON;
        }
        if (a == 3u) {
            return FBU_SPECIES_TRITON;
        }
        return -1;
    }
    if (z == 2u) {
        if (a == 3u) {
            return FBU_SPECIES_HE3;
        }
        if (a == 4u) {
            return FBU_SPECIES_HE4;
        }
        return -1;
    }
    return -1;
}

/**
 * Fill the dense ground-state nuclear mass table [MeV/c²] from the isotope
 * database; entries absent from the database stay 0 and close any channel
 * that would produce them.
 */
static void fill_mass_table(double *mass_mev) {
    unsigned int z;
    unsigned int a;
    double m;

    /* The neutron is not an isotope-database entry; use the registry constant
     * (same NIST source as the database-derived masses, consistent to
     * sub-keV with the Q-values computed below). */
    mass_mev[_za_idx(0u, 1u)] = OSH_PART_MASS_NEUTRON;

    for (z = 1u; z <= OSH_FERMI_BREAKUP_ZMAX; ++z) {
        for (a = z; a <= OSH_FERMI_BREAKUP_AMAX; ++a) {
            if (osh_particle_nuclear_mass_mev_from_za(z, a, &m)) {
                mass_mev[_za_idx(z, a)] = m;
            }
        }
    }
}

/**
 * Enumerate all binary decay channels parent(z, a) -> (z1, a1) + (z2, a2)
 * with both products present in the mass table.  Unordered pairs are listed
 * once (a1 <= a2; z1 <= z2 when a1 == a2).
 *
 * When @p pool is NULL only the total channel count is returned (sizing
 * pass); otherwise the channels and the per-parent offset/count tables are
 * written.
 */
static size_t enumerate_channels(double const *mass_mev,
                                 struct osh_fermi_channel *pool,
                                 uint16_t *chan_offset,
                                 uint16_t *chan_count) {
    size_t total;
    size_t idx;
    unsigned int z;
    unsigned int a;
    unsigned int z1;
    unsigned int a1;
    unsigned int z2;
    unsigned int a2;
    unsigned int nch;
    double mp;
    double m1;
    double m2;
    double mu;
    double pref;

    total = 0u;
    for (z = 1u; z <= OSH_FERMI_BREAKUP_ZMAX; ++z) {
        for (a = 2u; a <= OSH_FERMI_BREAKUP_AMAX; ++a) {
            idx = _za_idx(z, a);
            nch = 0u;
            mp = mass_mev[idx];
            if (mp > 0.0) {
                for (a1 = 1u; a1 <= a / 2u; ++a1) {
                    a2 = a - a1;
                    for (z1 = 0u; z1 <= z; ++z1) {
                        z2 = z - z1;
                        if (a1 == a2 && z1 > z2) {
                            continue; /* unordered pair already listed */
                        }
                        m1 = mass_mev[_za_idx(z1, a1)];
                        m2 = mass_mev[_za_idx(z2, a2)];
                        if (m1 <= 0.0 || m2 <= 0.0) {
                            continue; /* product nuclide does not exist */
                        }
                        if (pool != NULL) {
                            mu = m1 * m2 / (m1 + m2);
                            pref = (double) spin_degeneracy(z1, a1) * (double) spin_degeneracy(z2, a2) * mu * sqrt(mu);
                            if (z1 == z2 && a1 == a2) {
                                pref *= 0.5; /* identical products */
                            }
                            pool[total].q_mev = (float) (mp - m1 - m2);
                            pool[total].weight_prefactor = (float) pref;
                            pool[total].z1 = (uint8_t) z1;
                            pool[total].a1 = (uint8_t) a1;
                            pool[total].z2 = (uint8_t) z2;
                            pool[total].a2 = (uint8_t) a2;
                        }
                        ++total;
                        ++nch;
                    }
                }
            }
            if (chan_offset != NULL) {
                chan_offset[idx] = (uint16_t) (total - nch);
                chan_count[idx] = (uint16_t) nch;
            }
        }
    }
    return total;
}

/**
 * Sum of open-channel weights w = prefactor * sqrt(e_star + Q) for parent at
 * dense index @p idx and excitation @p e_star; 0 if no channel is open.
 */
static double open_channel_weight_sum(struct osh_nuclear_fermi_breakup const *model, size_t idx, double e_star) {
    double wsum;
    double ekin;
    uint16_t off;
    uint16_t cnt;
    uint16_t c;

    wsum = 0.0;
    off = model->chan_offset[idx];
    cnt = model->chan_count[idx];
    for (c = 0u; c < cnt; ++c) {
        ekin = e_star + (double) model->channel_pool[off + c].q_mev;
        if (ekin > 0.0) {
            wsum += (double) model->channel_pool[off + c].weight_prefactor * sqrt(ekin);
        }
    }
    return wsum;
}

/**
 * Append @p node to event_out->fragments[] as an unprocessed fragment.
 * If the fragment array is full the product is dropped; the associated rest
 * of the energy budget is lost — accepted residual of the truncation policy.
 */
static void append_unprocessed_fragment(struct fbu_work_item const *node, struct osh_nuclear_event *event_out) {
    struct osh_nuclear_fragment *f;

    if (event_out->n_fragments >= OSH_NUCLEAR_MAX_FRAGMENTS) {
        return;
    }
    f = &event_out->fragments[event_out->n_fragments];
    f->excitation_energy = node->e_star;
    f->p[0] = node->p[0];
    f->p[1] = node->p[1];
    f->p[2] = node->p[2];
    f->z = node->z;
    f->a = node->a;
    event_out->n_fragments += 1u;
}

/* ---- Public API ----------------------------------------------------------- */

enum osh_status osh_nuclear_fermi_breakup_compile(struct osh_nuclear_fermi_breakup *out) {
    static int const s_pdgs[FBU_SPECIES_COUNT] = {OSH_PART_PDG_NEUTRON,
                                                  OSH_PART_PDG_PROTON,
                                                  OSH_PART_PDG_DEUTERON,
                                                  OSH_PART_PDG_TRITON,
                                                  OSH_PART_PDG_HE3,
                                                  OSH_PART_PDG_HE4};
    size_t total;
    int i;

    if (out == NULL) {
        return OSH_EINVAL;
    }

    out->mass_mev = calloc(FBU_NDENSE, sizeof(*out->mass_mev));
    out->chan_offset = calloc(FBU_NDENSE, sizeof(*out->chan_offset));
    out->chan_count = calloc(FBU_NDENSE, sizeof(*out->chan_count));
    out->species = calloc(FBU_SPECIES_COUNT, sizeof(*out->species));
    if (out->mass_mev == NULL || out->chan_offset == NULL || out->chan_count == NULL || out->species == NULL) {
        osh_nuclear_fermi_breakup_free(out);
        return OSH_ENOMEM;
    }

    for (i = 0; i < FBU_SPECIES_COUNT; ++i) {
        if (!osh_particle_from_pdg(&out->species[i], s_pdgs[i])) {
            osh_nuclear_fermi_breakup_free(out);
            return OSH_ESTATE;
        }
    }

    fill_mass_table(out->mass_mev);

    total = enumerate_channels(out->mass_mev, NULL, NULL, NULL);
    if (total > (size_t) UINT16_MAX) {
        osh_nuclear_fermi_breakup_free(out);
        return OSH_ESTATE; /* chan_offset is uint16_t; cannot index the pool */
    }
    out->channel_pool = calloc(total > 0u ? total : 1u, sizeof(*out->channel_pool));
    if (out->channel_pool == NULL) {
        osh_nuclear_fermi_breakup_free(out);
        return OSH_ENOMEM;
    }
    out->nchannels = enumerate_channels(out->mass_mev, out->channel_pool, out->chan_offset, out->chan_count);

    return OSH_OK;
}

void osh_nuclear_fermi_breakup_free(struct osh_nuclear_fermi_breakup *m) {
    if (m == NULL) {
        return;
    }
    free(m->channel_pool);
    free(m->mass_mev);
    free(m->species);
    free(m->chan_offset);
    free(m->chan_count);
    m->channel_pool = NULL;
    m->mass_mev = NULL;
    m->species = NULL;
    m->chan_offset = NULL;
    m->chan_count = NULL;
    m->nchannels = 0u;
}

void osh_nuclear_fermi_breakup_step(struct osh_nuclear_fermi_breakup const *model,
                                    struct osh_nuclear_fragment const *fragment,
                                    struct osh_rng *rng,
                                    struct osh_nuclear_event *event_out) {
    struct fbu_work_item stack[OSH_FERMI_BREAKUP_AMAX];
    struct fbu_work_item node;
    struct osh_nuclear_fragment frag;
    struct osh_fermi_channel const *ch;
    struct osh_nuclear_secondary *sec;
    size_t idx;
    size_t n_emitted;
    int sp;
    int spec_idx;
    int room;
    uint16_t off;
    uint16_t cnt;
    uint16_t c;
    double wsum;
    double ekin;
    double threshold;
    double cumulative;
    double m_eff;
    double m1;
    double m2;
    double p_star;
    double cos_theta;
    double sin_theta;
    double cos_phi;
    double sin_phi;
    double p_cm[3];
    double e_cm;
    double e_lab;
    double m_node;
    double p_norm;

    /* Local copy first: the fragment usually aliases event_out->fragments[0],
     * which is overwritten once the break-up commits. */
    frag = *fragment;

    if (frag.a < 2u || frag.a > OSH_FERMI_BREAKUP_AMAX || frag.z < 1u || frag.z > OSH_FERMI_BREAKUP_ZMAX) {
        return; /* outside model domain: leave for the fragment pool */
    }
    idx = _za_idx(frag.z, frag.a);
    if (model->mass_mev[idx] <= 0.0) {
        return; /* unknown nuclide */
    }

    /* Top-level gate: no open channel means the prefragment de-excites by
     * gamma emission, which is not modelled — the event stays untouched and
     * E* is dropped (documented limitation). */
    if (open_channel_weight_sum(model, idx, frag.excitation_energy) <= 0.0) {
        return;
    }

    /* Commit: the prefragment slot is consumed; unprocessed residues are
     * appended back as they are encountered. */
    event_out->n_fragments = 0u;
    n_emitted = 0u;

    stack[0].p[0] = frag.p[0];
    stack[0].p[1] = frag.p[1];
    stack[0].p[2] = frag.p[2];
    stack[0].e_star = frag.excitation_energy;
    stack[0].z = (uint8_t) frag.z;
    stack[0].a = (uint8_t) frag.a;
    sp = 1;

    /* Sequential binary de-excitation.  The stack never exceeds AMAX entries:
     * the summed mass number over pending items is conserved at frag.a <= AMAX
     * and every item carries a >= 1. */
    while (sp > 0) {
        --sp;
        node = stack[sp];
        idx = _za_idx(node.z, node.a);

        /* Truncation policy: once the secondary array is full, stop splitting
         * and drain the remaining work to the unprocessed-fragment list. */
        room = (event_out->n_secondaries < OSH_NUCLEAR_MAX_SECONDARIES);

        wsum = 0.0;
        if (room) {
            wsum = open_channel_weight_sum(model, idx, node.e_star);
        }

        if (wsum <= 0.0) {
            spec_idx = final_species_index(node.z, node.a);
            if (room && spec_idx >= 0) {
                /* Whitelisted transportable product: emit as secondary.
                 * Kinematics use the table mass for consistency with the
                 * Q-values; the registry descriptor is what transport sees. */
                m_node = model->mass_mev[idx];
                p_norm = sqrt(node.p[0] * node.p[0] + node.p[1] * node.p[1] + node.p[2] * node.p[2]);
                sec = &event_out->secondaries[event_out->n_secondaries];
                sec->energy = sqrt(m_node * m_node + p_norm * p_norm) - m_node;
                if (p_norm > 0.0) {
                    sec->dir[0] = node.p[0] / p_norm;
                    sec->dir[1] = node.p[1] / p_norm;
                    sec->dir[2] = node.p[2] / p_norm;
                } else {
                    /* Product exactly at rest (vanishing measure): emit along
                     * +z; the zero kinetic energy makes the direction moot. */
                    sec->dir[0] = 0.0;
                    sec->dir[1] = 0.0;
                    sec->dir[2] = 1.0;
                }
                sec->species = &model->species[spec_idx];
                event_out->n_secondaries += 1u;
                ++n_emitted;
            } else {
                append_unprocessed_fragment(&node, event_out);
            }
            continue;
        }

        /* Select the decay channel by cumulative open-channel weight. */
        off = model->chan_offset[idx];
        cnt = model->chan_count[idx];
        ch = NULL;
        threshold = osh_rng_double(rng) * wsum;
        cumulative = 0.0;
        for (c = 0u; c < cnt; ++c) {
            ekin = node.e_star + (double) model->channel_pool[off + c].q_mev;
            if (ekin > 0.0) {
                ch = &model->channel_pool[off + c];
                cumulative += (double) ch->weight_prefactor * sqrt(ekin);
                if (threshold <= cumulative) {
                    break;
                }
            }
        }
        if (ch == NULL) {
            /* Unreachable when wsum > 0; defensive against table corruption. */
            append_unprocessed_fragment(&node, event_out);
            continue;
        }

        /* Two-body decay of the excited parent (effective mass M + E*),
         * isotropic in the parent rest frame, boosted to the lab. */
        m1 = model->mass_mev[_za_idx(ch->z1, ch->a1)];
        m2 = model->mass_mev[_za_idx(ch->z2, ch->a2)];
        m_eff = model->mass_mev[idx] + node.e_star;
        p_star = osh_kinematics_two_body_decay_p(m_eff, m1, m2);

        cos_theta = 2.0 * osh_rng_double(rng) - 1.0;
        sin_theta = sqrt(fmax(0.0, 1.0 - cos_theta * cos_theta));
        osh_kinematics_azimuth(rng, &cos_phi, &sin_phi);
        p_cm[0] = p_star * sin_theta * cos_phi;
        p_cm[1] = p_star * sin_theta * sin_phi;
        p_cm[2] = p_star * cos_theta;

        if (sp + 2 > (int) OSH_FERMI_BREAKUP_AMAX) {
            /* Unreachable by mass-number conservation; defensive guard. */
            append_unprocessed_fragment(&node, event_out);
            continue;
        }

        e_cm = sqrt(m1 * m1 + p_star * p_star);
        osh_kinematics_boost_to_lab(m_eff, node.p, e_cm, p_cm, &e_lab, stack[sp].p);
        stack[sp].e_star = 0.0;
        stack[sp].z = ch->z1;
        stack[sp].a = ch->a1;
        ++sp;

        p_cm[0] = -p_cm[0];
        p_cm[1] = -p_cm[1];
        p_cm[2] = -p_cm[2];
        e_cm = sqrt(m2 * m2 + p_star * p_star);
        osh_kinematics_boost_to_lab(m_eff, node.p, e_cm, p_cm, &e_lab, stack[sp].p);
        stack[sp].e_star = 0.0;
        stack[sp].z = ch->z2;
        stack[sp].a = ch->a2;
        ++sp;
    }

    if (n_emitted > 0u) {
        event_out->kind = OSH_NUCLEAR_EVENT_FRAGMENTATION;
    }
}
