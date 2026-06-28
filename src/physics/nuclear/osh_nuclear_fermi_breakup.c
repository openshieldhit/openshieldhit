#include "physics/nuclear/osh_nuclear_fermi_breakup.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "particle/osh_particle.h"
#include "particle/osh_particle_const.h"
#include "particle/osh_particle_pdg.h"
#include "physics/nuclear/osh_nuclear_handler.h"
#include "physics/osh_kinematics.h"
#include "random/osh_rng.h"

/* Dense (z, a) table size: z in [0, ZMAX], a in [0, AMAX]. */
#define FBU_NDENSE ((size_t)(OSH_FERMI_BREAKUP_ZMAX + 1) * (size_t)(OSH_FERMI_BREAKUP_AMAX + 1))

/* Maximum N-body partition size we enumerate.
 * Capped at N=3: the 3-body channels (e.g. 3-alpha from C-12) fix the main
 * high-E* multiplicity deficit without producing the T^(3N/2-5/2) runaway
 * that occurs for N>=4.  A full N>=4 model requires nuclear excited-state
 * data to give smooth interpolation (as in G4FermiBreakUp 9.1); that is
 * deferred to a future enhancement. */
#define FBU_NMAX 3

/* Final-product species indices into model->species. */
#define FBU_SPECIES_NEUTRON  0
#define FBU_SPECIES_PROTON   1
#define FBU_SPECIES_DEUTERON 2
#define FBU_SPECIES_TRITON   3
#define FBU_SPECIES_HE3      4
#define FBU_SPECIES_HE4      5
#define FBU_SPECIES_COUNT    6

/* One pending nuclide on the de-excitation work stack. */
struct fbu_work_item {
    double p[3];   /* lab momentum [MeV/c] */
    double e_star; /* excitation energy [MeV]; 0 for break-up products */
    uint8_t z;
    uint8_t a;
};

/* Ground-state spin degeneracy entry g = 2J+1 (ENSDF/TUNL adopted spins). */
struct fbu_spin_entry {
    uint8_t z;
    uint8_t a;
    uint8_t g;
};

/* Nuclides not listed default to g = 1. */
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

/* Whitelist species (z, a) for N>=3 partitions, in canonical enumeration order
 * (ascending species index = ascending a, then z). */
static uint8_t const s_wl_z[FBU_SPECIES_COUNT] = {0u, 1u, 1u, 1u, 2u, 2u};
static uint8_t const s_wl_a[FBU_SPECIES_COUNT] = {1u, 1u, 2u, 3u, 3u, 4u};

/* Gamma(3N/2 - 3/2) for N = 2..FBU_NMAX.
 * N=2: Gamma(3/2) = sqrt(pi)/2 ~ 0.8862
 * N=3: Gamma(3)   = 2 */
static double const s_gamma[FBU_NMAX + 1] = {
    0.0, 0.0,
    0.88622692545275801, /* N=2 */
    2.0,                 /* N=3 */
};

/* Particle-unstable nuclear excited states included as binary-channel products.
 * Each entry (z, a, exc_kev) is added to the N=2 channel table with effective
 * mass M_gs + exc_kev/1000.  When selected at runtime, the product is pushed
 * to the work stack with e_star = exc_kev/1000 MeV, enabling further decay.
 *
 * Selected levels (ENSDF/TUNL):
 *  Li-7*(4.63 MeV): decays to t+He-4 (Q=-0.925, T=3.7 MeV).  Key for
 *    chains C-12→p+B-11*(14)→p+He-4+Li-7*(4.63)→p+He-4+t+He-4 (mult=4).
 *  B-11*(8.92 MeV): decays to He-4+Li-7 (T≈0.18 MeV, just open).
 *  B-11*(9.27 MeV): decays to He-4+Li-7 (T≈0.53 MeV, stronger).
 *  B-11*(14.0 MeV): decays to He-4+Li-7*(4.63) (T≈0.6 MeV), enabling the
 *    mult=4 chain.  Approximates several real B-11 levels near 14 MeV.
 *  Be-9*(2.43 MeV): decays to n+Be-8(unstable)→n+He-4+He-4 (mult=3). */
struct fbu_exc_entry {
    uint8_t z;
    uint8_t a;
    uint16_t exc_kev;
};
static struct fbu_exc_entry const s_exc_states[] = {
    {3u,  7u,  4630u},
    {5u, 11u,  8920u},
    {5u, 11u,  9270u},
    {5u, 11u, 14000u},
    {4u,  9u,  2430u},
};

/* Dense (z, a) index. */
static inline size_t _za_idx(unsigned int z, unsigned int a) {
    return (size_t)z * (size_t)(OSH_FERMI_BREAKUP_AMAX + 1) + (size_t)a;
}

/* Ground-state spin degeneracy 2J+1; 1 if not tabulated. */
static unsigned int spin_degeneracy(unsigned int z, unsigned int a) {
    size_t i;
    for (i = 0u; i < sizeof(s_spin_table) / sizeof(s_spin_table[0]); ++i) {
        if (s_spin_table[i].z == z && s_spin_table[i].a == a) {
            return s_spin_table[i].g;
        }
    }
    return 1u;
}

/* Index into model->species if (z, a) is a whitelisted product; -1 otherwise. */
static int final_species_index(unsigned int z, unsigned int a) {
    if (a == 1u && z == 0u) return FBU_SPECIES_NEUTRON;
    if (z == 1u) {
        if (a == 1u) return FBU_SPECIES_PROTON;
        if (a == 2u) return FBU_SPECIES_DEUTERON;
        if (a == 3u) return FBU_SPECIES_TRITON;
        return -1;
    }
    if (z == 2u) {
        if (a == 3u) return FBU_SPECIES_HE3;
        if (a == 4u) return FBU_SPECIES_HE4;
        return -1;
    }
    return -1;
}

/* Fill ground-state nuclear mass table [MeV/c^2] from the isotope database. */
static void fill_mass_table(double *mass_mev) {
    unsigned int z, a;
    double m;

    mass_mev[_za_idx(0u, 1u)] = OSH_PART_MASS_NEUTRON;
    for (z = 1u; z <= OSH_FERMI_BREAKUP_ZMAX; ++z) {
        for (a = z; a <= OSH_FERMI_BREAKUP_AMAX; ++a) {
            if (osh_particle_nuclear_mass_mev_from_za(z, a, &m)) {
                mass_mev[_za_idx(z, a)] = m;
            }
        }
    }
}

/* Append an unprocessed fragment.  Drops silently if the array is full. */
static void append_unprocessed_fragment(struct fbu_work_item const *node,
                                         struct osh_nuclear_event *event_out) {
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

/* ---- Partition enumeration ------------------------------------------------ */

/* Context for the recursive N>=3 whitelist-only partition enumeration. */
struct enum_ctx {
    double const *mass_mev;
    double m_parent;
    int n_total;       /* target N for this recursion subtree */
    uint8_t buf_z[FBU_NMAX];
    uint8_t buf_a[FBU_NMAX];
    /* output pointers (NULL during size pass) */
    struct osh_fermi_partition *part_pool;
    struct osh_fermi_frag_spec *fspec_pool;
    size_t *part_total;
    size_t *fspec_total;
};

/* Compute weight prefactor for the partition currently assembled in ctx->buf_[]. */
static float compute_nbody_prefactor(struct enum_ctx const *ctx, int n) {
    double q;
    double prod_m32;
    double prod_g;
    int counts[FBU_SPECIES_COUNT];
    int i, wsp, c;
    double id_denom;
    double mi;

    memset(counts, 0, sizeof(counts));
    q = ctx->m_parent;
    prod_m32 = 1.0;
    prod_g = 1.0;
    for (i = 0; i < n; ++i) {
        mi = ctx->mass_mev[_za_idx(ctx->buf_z[i], ctx->buf_a[i])];
        q -= mi;
        prod_m32 *= mi * sqrt(mi);
        prod_g *= spin_degeneracy(ctx->buf_z[i], ctx->buf_a[i]);
        for (wsp = 0; wsp < FBU_SPECIES_COUNT; ++wsp) {
            if (s_wl_z[wsp] == ctx->buf_z[i] && s_wl_a[wsp] == ctx->buf_a[i]) {
                counts[wsp]++;
                break;
            }
        }
    }
    id_denom = 1.0;
    for (wsp = 0; wsp < FBU_SPECIES_COUNT; ++wsp) {
        for (c = counts[wsp]; c > 1; --c) {
            id_denom *= c;
        }
    }
    (void)q; /* q is used by caller via part_pool write */
    return (float)(prod_g / id_denom * prod_m32 / s_gamma[n]);
}

/* Recursive enumeration of whitelist-only N-body partitions.
 * k_rem: number of fragments still to place.
 * z_rem, a_rem: remaining charge and mass to distribute.
 * min_sp: minimum species index (canonical non-decreasing order).
 * buf_pos: next write position in ctx->buf_[]. */
static void enum_wl_rec(struct enum_ctx *ctx, int k_rem, int z_rem, int a_rem,
                         int min_sp, int buf_pos) {
    int sp, zi, ai, n, i;
    double q;

    if (k_rem == 0) {
        if (z_rem != 0 || a_rem != 0) {
            return; /* charge/mass not exhausted: invalid */
        }
        n = ctx->n_total;
        if (ctx->part_pool != NULL) {
            q = ctx->m_parent;
            for (i = 0; i < n; ++i) {
                q -= ctx->mass_mev[_za_idx(ctx->buf_z[i], ctx->buf_a[i])];
            }
            ctx->part_pool[*ctx->part_total].q_mev = (float)q;
            ctx->part_pool[*ctx->part_total].weight_prefactor = compute_nbody_prefactor(ctx, n);
            ctx->part_pool[*ctx->part_total].fspec_offset = (uint32_t)*ctx->fspec_total;
            ctx->part_pool[*ctx->part_total].n_frags = (uint8_t)n;
            ctx->part_pool[*ctx->part_total]._pad[0] = 0;
            ctx->part_pool[*ctx->part_total]._pad[1] = 0;
            ctx->part_pool[*ctx->part_total]._pad[2] = 0;
            for (i = 0; i < n; ++i) {
                ctx->fspec_pool[*ctx->fspec_total].z = ctx->buf_z[i];
                ctx->fspec_pool[*ctx->fspec_total].a = ctx->buf_a[i];
                (*ctx->fspec_total)++;
            }
        } else {
            *ctx->fspec_total += (size_t)ctx->n_total;
        }
        (*ctx->part_total)++;
        return;
    }

    /* Try each whitelist species at or after min_sp. */
    for (sp = min_sp; sp < FBU_SPECIES_COUNT; ++sp) {
        zi = s_wl_z[sp];
        ai = s_wl_a[sp];
        if (zi > z_rem || ai > a_rem) {
            continue;
        }
        /* Pruning: remaining k_rem fragments each need at least ai mass. */
        if (ai * k_rem > a_rem) {
            continue;
        }
        ctx->buf_z[buf_pos] = (uint8_t)zi;
        ctx->buf_a[buf_pos] = (uint8_t)ai;
        enum_wl_rec(ctx, k_rem - 1, z_rem - zi, a_rem - ai, sp, buf_pos + 1);
    }
}

/* Enumerate all N-body partitions for every valid parent (z, a).
 * N=2: full mass-table binary pairs.
 * N=3..NMAX: whitelist-only.
 * Pass pool==NULL, fspec==NULL for the sizing pass. */
static size_t enumerate_all_partitions(double const *mass_mev,
                                        struct osh_fermi_partition *part_pool,
                                        struct osh_fermi_frag_spec *fspec_pool,
                                        uint32_t *part_offset,
                                        uint16_t *part_count,
                                        size_t *nfspecs_out) {
    size_t part_total;
    size_t fspec_total;
    unsigned int z, a, z1, a1, z2, a2;
    double mp, m1, m2, pref;
    size_t part_start;
    size_t idx;
    int n;

    part_total = 0u;
    fspec_total = 0u;
    for (z = 1u; z <= OSH_FERMI_BREAKUP_ZMAX; ++z) {
        for (a = 2u; a <= OSH_FERMI_BREAKUP_AMAX; ++a) {
            idx = _za_idx(z, a);
            mp = mass_mev[idx];
            part_start = part_total;

            if (mp > 0.0) {
                /* N=2 binary pairs from the full mass table. */
                for (a1 = 1u; a1 <= a / 2u; ++a1) {
                    a2 = a - a1;
                    for (z1 = 0u; z1 <= z; ++z1) {
                        z2 = z - z1;
                        if (a1 == a2 && z1 > z2) {
                            continue; /* already listed in the other order */
                        }
                        m1 = mass_mev[_za_idx(z1, a1)];
                        m2 = mass_mev[_za_idx(z2, a2)];
                        if (m1 <= 0.0 || m2 <= 0.0) {
                            continue;
                        }
                        if (part_pool != NULL) {
                            /* Weight prefactor uses prod(m_i^3/2)/Gamma(3/2)
                             * for consistency with the N>=3 formula. */
                            pref = (double)spin_degeneracy(z1, a1)
                                        * (double)spin_degeneracy(z2, a2)
                                        * m1 * sqrt(m1) * m2 * sqrt(m2)
                                        / s_gamma[2];
                            if (z1 == z2 && a1 == a2) {
                                pref *= 0.5; /* identical products: 1/2! */
                            }
                            part_pool[part_total].q_mev = (float)(mp - m1 - m2);
                            part_pool[part_total].weight_prefactor = (float)pref;
                            part_pool[part_total].fspec_offset = (uint32_t)fspec_total;
                            part_pool[part_total].n_frags = 2u;
                            part_pool[part_total]._pad[0] = 0;
                            part_pool[part_total]._pad[1] = 0;
                            part_pool[part_total]._pad[2] = 0;
                            fspec_pool[fspec_total].z = (uint8_t)z1;
                            fspec_pool[fspec_total].a = (uint8_t)a1;
                            fspec_pool[fspec_total].exc_kev = 0u;
                            fspec_pool[fspec_total + 1u].z = (uint8_t)z2;
                            fspec_pool[fspec_total + 1u].a = (uint8_t)a2;
                            fspec_pool[fspec_total + 1u].exc_kev = 0u;
                        }
                        part_total++;
                        fspec_total += 2u;
                    }
                }

                /* N=2 channels with particle-unstable excited-state products.
                 * Each s_exc_states entry (z_exc, a_exc, exc_kev) is the
                 * excited fragment; the other product (z_other, a_other) is
                 * ground-state.  Q uses the excited mass M_gs + exc_mev.
                 * When the excited product is later popped from the work
                 * stack its e_star = exc_kev/1000 enables further decay. */
                {
                    size_t ei;
                    size_t nexc;
                    unsigned int z_exc, a_exc, z_oth, a_oth;
                    uint16_t exc_kev;
                    double exc_mev;
                    double m_exc, m_oth, m_exc_eff;

                    nexc = sizeof(s_exc_states) / sizeof(s_exc_states[0]);
                    for (ei = 0u; ei < nexc; ++ei) {
                        z_exc = s_exc_states[ei].z;
                        a_exc = s_exc_states[ei].a;
                        exc_kev = s_exc_states[ei].exc_kev;
                        exc_mev = (double)exc_kev / 1000.0;

                        if (a_exc >= a || z_exc > z) {
                            continue;
                        }
                        a_oth = a - a_exc;
                        z_oth = z - z_exc;
                        if (a_oth < 1u) {
                            continue;
                        }
                        m_exc = mass_mev[_za_idx(z_exc, a_exc)];
                        m_oth = mass_mev[_za_idx(z_oth, a_oth)];
                        if (m_exc <= 0.0 || m_oth <= 0.0) {
                            continue;
                        }
                        m_exc_eff = m_exc + exc_mev;
                        if (part_pool != NULL) {
                            pref = (double)spin_degeneracy(z_oth, a_oth)
                                 * (double)spin_degeneracy(z_exc, a_exc)
                                 * m_oth * sqrt(m_oth)
                                 * m_exc_eff * sqrt(m_exc_eff)
                                 / s_gamma[2];
                            part_pool[part_total].q_mev = (float)(mp - m_oth - m_exc_eff);
                            part_pool[part_total].weight_prefactor = (float)pref;
                            part_pool[part_total].fspec_offset = (uint32_t)fspec_total;
                            part_pool[part_total].n_frags = 2u;
                            part_pool[part_total]._pad[0] = 0;
                            part_pool[part_total]._pad[1] = 0;
                            part_pool[part_total]._pad[2] = 0;
                            fspec_pool[fspec_total].z = (uint8_t)z_oth;
                            fspec_pool[fspec_total].a = (uint8_t)a_oth;
                            fspec_pool[fspec_total].exc_kev = 0u;
                            fspec_pool[fspec_total + 1u].z = (uint8_t)z_exc;
                            fspec_pool[fspec_total + 1u].a = (uint8_t)a_exc;
                            fspec_pool[fspec_total + 1u].exc_kev = exc_kev;
                        }
                        part_total++;
                        fspec_total += 2u;
                    }
                }

                /* N=3..NMAX whitelist-only partitions. */
                {
                    struct enum_ctx ctx;
                    ctx.mass_mev = mass_mev;
                    ctx.m_parent = mp;
                    ctx.part_pool = part_pool;
                    ctx.fspec_pool = fspec_pool;
                    ctx.part_total = &part_total;
                    ctx.fspec_total = &fspec_total;
                    for (n = 3; n <= FBU_NMAX; ++n) {
                        if (n > (int)a) {
                            break; /* more fragments than nucleons */
                        }
                        ctx.n_total = n;
                        enum_wl_rec(&ctx, n, (int)z, (int)a, 0, 0);
                    }
                }
            }

            if (part_offset != NULL) {
                part_offset[idx] = (uint32_t)part_start;
                part_count[idx] = (uint16_t)(part_total - part_start);
            }
        }
    }

    if (nfspecs_out != NULL) {
        *nfspecs_out = fspec_total;
    }
    return part_total;
}

/* ---- Kopylov N-body kinematics ------------------------------------------- */

/* Sample chi from pdf proportional to chi^(3K-5) * (1-chi) via rejection.
 * Used by the Kopylov algorithm for K >= 2. */
static double beta_kopylov(int k, struct osh_rng *rng) {
    int npow;
    double xn;
    double fmax;
    double chi, f;

    npow = 3 * k - 5;
    xn = (double)npow;
    fmax = sqrt(pow(xn / (xn + 1.0), npow) / (xn + 1.0));
    do {
        chi = osh_rng_double(rng);
        f = sqrt(pow(chi, npow) * (1.0 - chi));
    } while (fmax * osh_rng_double(rng) > f);
    return chi;
}

/* Generate N-body phase space momenta (Kopylov algorithm, cf. G4FermiPhaseSpaceDecay).
 * m[0..n-1]: fragment rest masses [MeV].
 * m_eff: parent effective mass (m_gs + E*) [MeV].
 * p_parent[3]: parent lab momentum [MeV/c].
 * Output: p_out[i][3] lab momenta; e_tot[i] total lab energies (includes rest mass). */
static void kopylov_nbody(int n, double const m[], double m_eff,
                           double const p_parent[3], struct osh_rng *rng,
                           double p_out[][3], double e_tot[]) {
    double mu;          /* running remaining mass (unplaced fragments) */
    double mass;        /* current subsystem effective mass */
    double t;           /* available kinetic energy for current subsystem */
    double rest_mass;   /* effective mass of the "rest" pseudo-particle */
    double p_rest[3];   /* lab momentum of the rest pseudo-particle */
    double rest_m_cur;  /* rest mass of the rest pseudo-particle (for boost) */
    double p_cm[3];
    double p_cm_mag;
    double cos_theta, sin_theta, cos_phi, sin_phi;
    double e_frag_cm, e_rest_cm;
    double e_frag_lab, e_rest_lab;
    double p_rest_new[3];
    double p2;
    int k;

    mu = 0.0;
    for (k = 0; k < n; ++k) {
        mu += m[k];
    }

    mass = m_eff;
    t = fmax(0.0, m_eff - mu);

    p_rest[0] = p_parent[0];
    p_rest[1] = p_parent[1];
    p_rest[2] = p_parent[2];
    rest_m_cur = m_eff;

    for (k = n - 1; k > 0; --k) {
        mu -= m[k];
        if (k > 1) {
            t *= beta_kopylov(k, rng);
        } else {
            t = 0.0; /* last step: no remaining kinetic energy for fragment 0 */
        }
        rest_mass = mu + t;

        p_cm_mag = osh_kinematics_two_body_decay_p(mass, m[k], rest_mass);

        cos_theta = 2.0 * osh_rng_double(rng) - 1.0;
        sin_theta = sqrt(fmax(0.0, 1.0 - cos_theta * cos_theta));
        osh_kinematics_azimuth(rng, &cos_phi, &sin_phi);
        p_cm[0] = p_cm_mag * sin_theta * cos_phi;
        p_cm[1] = p_cm_mag * sin_theta * sin_phi;
        p_cm[2] = p_cm_mag * cos_theta;

        e_frag_cm = sqrt(m[k] * m[k] + p_cm_mag * p_cm_mag);
        osh_kinematics_boost_to_lab(rest_m_cur, p_rest, e_frag_cm, p_cm,
                                     &e_frag_lab, p_out[k]);
        e_tot[k] = e_frag_lab;

        /* Boost rest pseudo-particle (opposite direction) to lab. */
        p_cm[0] = -p_cm[0];
        p_cm[1] = -p_cm[1];
        p_cm[2] = -p_cm[2];
        e_rest_cm = sqrt(rest_mass * rest_mass + p_cm_mag * p_cm_mag);
        osh_kinematics_boost_to_lab(rest_m_cur, p_rest, e_rest_cm, p_cm,
                                     &e_rest_lab, p_rest_new);
        p_rest[0] = p_rest_new[0];
        p_rest[1] = p_rest_new[1];
        p_rest[2] = p_rest_new[2];
        rest_m_cur = rest_mass;
        mass = rest_mass;
    }

    /* Fragment 0: whatever remains in the rest pseudo-particle. */
    p_out[0][0] = p_rest[0];
    p_out[0][1] = p_rest[1];
    p_out[0][2] = p_rest[2];
    p2 = p_rest[0] * p_rest[0] + p_rest[1] * p_rest[1] + p_rest[2] * p_rest[2];
    e_tot[0] = sqrt(m[0] * m[0] + p2);
}

/* ---- Weight selection helper --------------------------------------------- */

/* Compute total open-partition weight for parent at dense index idx and the
 * given excitation energy.  Returns 0 if no channel is open. */
static double partition_weight_sum(struct osh_nuclear_fermi_breakup const *model,
                                    size_t idx, double e_star) {
    double wsum;
    uint32_t off;
    uint16_t cnt;
    uint16_t i;
    double t;
    int nf;

    wsum = 0.0;
    off = model->part_offset[idx];
    cnt = model->part_count[idx];
    for (i = 0u; i < cnt; ++i) {
        t = e_star + (double)model->part_pool[off + i].q_mev;
        if (t > 0.0) {
            nf = model->part_pool[off + i].n_frags;
            wsum += (double)model->part_pool[off + i].weight_prefactor
                  * pow(t, 1.5 * nf - 2.5);
        }
    }
    return wsum;
}

/* ---- Public API ----------------------------------------------------------- */

enum osh_status osh_nuclear_fermi_breakup_compile(struct osh_nuclear_fermi_breakup *out) {
    static int const s_pdgs[FBU_SPECIES_COUNT] = {
        OSH_PART_PDG_NEUTRON, OSH_PART_PDG_PROTON,
        OSH_PART_PDG_DEUTERON, OSH_PART_PDG_TRITON,
        OSH_PART_PDG_HE3, OSH_PART_PDG_HE4};
    size_t nfspecs;
    size_t npart;
    int i;

    if (out == NULL) {
        return OSH_EINVAL;
    }

    out->mass_mev = calloc(FBU_NDENSE, sizeof(*out->mass_mev));
    out->part_offset = calloc(FBU_NDENSE, sizeof(*out->part_offset));
    out->part_count = calloc(FBU_NDENSE, sizeof(*out->part_count));
    out->species = calloc(FBU_SPECIES_COUNT, sizeof(*out->species));
    if (!out->mass_mev || !out->part_offset || !out->part_count || !out->species) {
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

    /* Sizing pass. */
    npart = enumerate_all_partitions(out->mass_mev, NULL, NULL, NULL, NULL, &nfspecs);

    out->part_pool = calloc(npart > 0u ? npart : 1u, sizeof(*out->part_pool));
    out->fspec_pool = calloc(nfspecs > 0u ? nfspecs : 1u, sizeof(*out->fspec_pool));
    if (!out->part_pool || !out->fspec_pool) {
        osh_nuclear_fermi_breakup_free(out);
        return OSH_ENOMEM;
    }

    /* Fill pass. */
    out->npartitions = enumerate_all_partitions(out->mass_mev, out->part_pool,
                                                 out->fspec_pool, out->part_offset,
                                                 out->part_count, &out->nfspecs);

    return OSH_OK;
}

void osh_nuclear_fermi_breakup_free(struct osh_nuclear_fermi_breakup *m) {
    if (m == NULL) {
        return;
    }
    free(m->mass_mev);
    free(m->species);
    free(m->part_pool);
    free(m->fspec_pool);
    free(m->part_offset);
    free(m->part_count);
    memset(m, 0, sizeof(*m));
}

void osh_nuclear_fermi_breakup_step(struct osh_nuclear_fermi_breakup const *model,
                                     struct osh_nuclear_fragment const *fragment,
                                     struct osh_rng *rng,
                                     struct osh_nuclear_event *event_out) {
    /* Work stack for N=2 chain decays (He-5→αn, Li-5→αp, Be-8→2α). */
    struct fbu_work_item stack[OSH_FERMI_BREAKUP_AMAX];
    struct fbu_work_item node;
    struct fbu_work_item tmp;
    struct osh_nuclear_fragment frag;
    struct osh_nuclear_secondary *sec;
    double m_arr[FBU_NMAX];
    double p_out[FBU_NMAX][3];
    double e_tot[FBU_NMAX];
    size_t idx;
    size_t n_emitted;
    int sp;
    int spec_idx;
    int room;
    int fi;
    int nf;
    uint32_t off;
    uint32_t fspec_off;
    uint16_t cnt;
    uint16_t exc1, exc2;
    uint16_t i, c;
    uint8_t z1, a1, z2, a2;
    uint8_t zi, ai;
    double wsum, threshold, cumulative, t, w;
    double m_eff, m1, m2;
    double mi;
    double p_star, cos_theta, sin_theta, cos_phi, sin_phi;
    double p_cm[3], e_cm, e_lab;
    double m_node, p_norm, p2;

    frag = *fragment; /* local copy: fragment may alias event_out->fragments[0] */

    if (frag.a < 2u || frag.a > OSH_FERMI_BREAKUP_AMAX
        || frag.z < 1u || frag.z > OSH_FERMI_BREAKUP_ZMAX) {
        return;
    }
    idx = _za_idx(frag.z, frag.a);
    if (model->mass_mev[idx] <= 0.0) {
        return;
    }

    /* Gate: any open partition? */
    if (partition_weight_sum(model, idx, frag.excitation_energy) <= 0.0) {
        return;
    }

    /* Commit: consume the prefragment slot. */
    event_out->n_fragments = 0u;
    n_emitted = 0u;

    /* Push the parent as the initial work item. */
    stack[0].p[0] = frag.p[0];
    stack[0].p[1] = frag.p[1];
    stack[0].p[2] = frag.p[2];
    stack[0].e_star = frag.excitation_energy;
    stack[0].z = (uint8_t)frag.z;
    stack[0].a = (uint8_t)frag.a;
    sp = 1;

    while (sp > 0) {
        --sp;
        node = stack[sp];
        idx = _za_idx(node.z, node.a);

        room = (event_out->n_secondaries < OSH_NUCLEAR_MAX_SECONDARIES);
        wsum = 0.0;
        if (room) {
            wsum = partition_weight_sum(model, idx, node.e_star);
        }

        if (wsum <= 0.0) {
            /* No open channel (or secondaries full): route to terminal. */
            spec_idx = final_species_index(node.z, node.a);
            if (room && spec_idx >= 0) {
                /* Whitelist species: emit as transportable secondary. */
                m_node = model->mass_mev[idx];
                p2 = node.p[0]*node.p[0] + node.p[1]*node.p[1] + node.p[2]*node.p[2];
                p_norm = sqrt(p2);
                sec = &event_out->secondaries[event_out->n_secondaries];
                sec->energy = sqrt(m_node * m_node + p2) - m_node;
                if (p_norm > 0.0) {
                    sec->dir[0] = node.p[0] / p_norm;
                    sec->dir[1] = node.p[1] / p_norm;
                    sec->dir[2] = node.p[2] / p_norm;
                } else {
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

        /* Select one partition by cumulative weight. */
        off = model->part_offset[idx];
        cnt = model->part_count[idx];
        threshold = osh_rng_double(rng) * wsum;
        cumulative = 0.0;
        c = 0;
        for (i = 0u; i < cnt; ++i) {
            t = node.e_star + (double)model->part_pool[off + i].q_mev;
            if (t > 0.0) {
                nf = model->part_pool[off + i].n_frags;
                w = (double)model->part_pool[off + i].weight_prefactor
                  * pow(t, 1.5 * nf - 2.5);
                cumulative += w;
                c = i;
                if (threshold <= cumulative) {
                    break;
                }
            }
        }

        nf = model->part_pool[off + c].n_frags;
        fspec_off = model->part_pool[off + c].fspec_offset;

        if (nf == 2) {
            /* Two-body decay: push both products on the work stack. */
            z1 = model->fspec_pool[fspec_off].z;
            a1 = model->fspec_pool[fspec_off].a;
            exc1 = model->fspec_pool[fspec_off].exc_kev;
            z2 = model->fspec_pool[fspec_off + 1u].z;
            a2 = model->fspec_pool[fspec_off + 1u].a;
            exc2 = model->fspec_pool[fspec_off + 1u].exc_kev;

            /* Use excited effective masses for kinematics when exc > 0. */
            m1 = model->mass_mev[_za_idx(z1, a1)] + (double)exc1 / 1000.0;
            m2 = model->mass_mev[_za_idx(z2, a2)] + (double)exc2 / 1000.0;
            m_eff = model->mass_mev[idx] + node.e_star;
            p_star = osh_kinematics_two_body_decay_p(m_eff, m1, m2);

            cos_theta = 2.0 * osh_rng_double(rng) - 1.0;
            sin_theta = sqrt(fmax(0.0, 1.0 - cos_theta * cos_theta));
            osh_kinematics_azimuth(rng, &cos_phi, &sin_phi);
            p_cm[0] = p_star * sin_theta * cos_phi;
            p_cm[1] = p_star * sin_theta * sin_phi;
            p_cm[2] = p_star * cos_theta;

            if (sp + 2 > (int)OSH_FERMI_BREAKUP_AMAX) {
                /* Guard: mass conservation ensures this never fires. */
                append_unprocessed_fragment(&node, event_out);
                continue;
            }

            e_cm = sqrt(m1 * m1 + p_star * p_star);
            osh_kinematics_boost_to_lab(m_eff, node.p, e_cm, p_cm, &e_lab, stack[sp].p);
            stack[sp].e_star = (double)exc1 / 1000.0;
            stack[sp].z = z1;
            stack[sp].a = a1;
            ++sp;

            p_cm[0] = -p_cm[0];
            p_cm[1] = -p_cm[1];
            p_cm[2] = -p_cm[2];
            e_cm = sqrt(m2 * m2 + p_star * p_star);
            osh_kinematics_boost_to_lab(m_eff, node.p, e_cm, p_cm, &e_lab, stack[sp].p);
            stack[sp].e_star = (double)exc2 / 1000.0;
            stack[sp].z = z2;
            stack[sp].a = a2;
            ++sp;

        } else {
            /* N>=3 partition: Kopylov kinematics, emit all directly as secondaries. */
            m_eff = model->mass_mev[idx] + node.e_star;
            for (fi = 0; fi < nf; ++fi) {
                m_arr[fi] = model->mass_mev[_za_idx(
                    model->fspec_pool[fspec_off + fi].z,
                    model->fspec_pool[fspec_off + fi].a)];
            }
            kopylov_nbody(nf, m_arr, m_eff, node.p, rng, p_out, e_tot);

            for (fi = 0; fi < nf && room; ++fi) {
                zi = model->fspec_pool[fspec_off + fi].z;
                ai = model->fspec_pool[fspec_off + fi].a;
                spec_idx = final_species_index(zi, ai);
                if (spec_idx < 0) {
                    /* Shouldn't happen for whitelist-only N>=3 partitions. */
                    tmp.z = zi;
                    tmp.a = ai;
                    tmp.e_star = 0.0;
                    tmp.p[0] = p_out[fi][0];
                    tmp.p[1] = p_out[fi][1];
                    tmp.p[2] = p_out[fi][2];
                    append_unprocessed_fragment(&tmp, event_out);
                    continue;
                }
                mi = m_arr[fi];
                p2 = p_out[fi][0]*p_out[fi][0] + p_out[fi][1]*p_out[fi][1]
                   + p_out[fi][2]*p_out[fi][2];
                p_norm = sqrt(p2);
                sec = &event_out->secondaries[event_out->n_secondaries];
                sec->energy = sqrt(mi * mi + p2) - mi;
                if (p_norm > 0.0) {
                    sec->dir[0] = p_out[fi][0] / p_norm;
                    sec->dir[1] = p_out[fi][1] / p_norm;
                    sec->dir[2] = p_out[fi][2] / p_norm;
                } else {
                    sec->dir[0] = 0.0;
                    sec->dir[1] = 0.0;
                    sec->dir[2] = 1.0;
                }
                sec->species = &model->species[spec_idx];
                event_out->n_secondaries += 1u;
                ++n_emitted;
                room = (event_out->n_secondaries < OSH_NUCLEAR_MAX_SECONDARIES);
            }
        }
    }

    if (n_emitted > 0u) {
        event_out->kind = OSH_NUCLEAR_EVENT_FRAGMENTATION;
    }
}
