#include "physics/nuclear/osh_nuclear_fermi_breakup.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "openshieldhit/const.h"
#include "particle/osh_particle.h"
#include "particle/osh_particle_const.h"
#include "particle/osh_particle_pdg.h"
#include "physics/nuclear/osh_nuclear_handler.h"
#include "physics/osh_kinematics.h"
#include "random/osh_rng.h"

/* Dense (z, a) table size: z in [0, ZMAX], a in [0, AMAX]. */
#define FBU_NDENSE ((size_t) (OSH_FERMI_BREAKUP_ZMAX + 1) * (size_t) (OSH_FERMI_BREAKUP_AMAX + 1))

/* Maximum N-body partition size we enumerate (ground-state fragments). */
#define FBU_NMAX OSH_FERMI_BREAKUP_NMAX

/* Work-stack depth: a partition pushes up to FBU_NMAX products and the sum of
 * nucleons across all pending items can never exceed the parent A (<= AMAX),
 * so AMAX + FBU_NMAX is a safe upper bound. */
#define FBU_WORK_STACK (OSH_FERMI_BREAKUP_AMAX + FBU_NMAX)

/* Final-product species indices into model->species. */
enum fbu_species {
    FBU_SPECIES_NEUTRON = 0,
    FBU_SPECIES_PROTON,
    FBU_SPECIES_DEUTERON,
    FBU_SPECIES_TRITON,
    FBU_SPECIES_HE3,
    FBU_SPECIES_HE4,
    FBU_SPECIES_COUNT
};

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

/* Gamma(3N/2 - 3/2) for N = 2..FBU_NMAX.
 * N=2: Gamma(3/2)  = sqrt(pi)/2
 * N=3: Gamma(3)
 * N=4: Gamma(9/2)
 * N=5: Gamma(6)
 * N=6: Gamma(15/2) */
static double const s_gamma[FBU_NMAX + 1] = {
    0.0,
    0.0,
    0.88622692545275801,  /* N=2 */
    2.0,                  /* N=3 */
    11.63172839656744790, /* N=4 */
    120.0,                /* N=5 */
    1871.25430573160160,  /* N=6 */
};

/* Dense (z, a) index. */
static inline size_t _za_idx(unsigned int z, unsigned int a) {
    return ((size_t) z * (size_t) (OSH_FERMI_BREAKUP_AMAX + 1)) + (size_t) a;
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

/* Index into model->species if (z, a) is a whitelisted transportable product;
 * -1 otherwise. */
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

/* Break-up free-volume / density-of-states coefficient C such that the
 * N-fragment phase-space weight scales as C^(N-1).
 *   C = V / ((2 pi)^(3/2) * (hbar c)^3),   V = (4/3) pi r0^3 A   [MeV^-3]
 * A is the parent mass number (the break-up volume is a property of the
 * decaying nucleus). */
static double volume_coeff(unsigned int a_parent) {
    double r0 = OSH_FERMI_BREAKUP_R0_FM;
    double v = (4.0 / 3.0) * OSH_M_PI * r0 * r0 * r0 * (double) a_parent;
    double denom = pow(2.0 * OSH_M_PI, 1.5) * OSH_HBARC * OSH_HBARC * OSH_HBARC;
    return v / denom;
}

/* Fill ground-state nuclear mass table [MeV/c^2] from the isotope database. */
static void fill_mass_table(double *mass_mev) {
    unsigned int z; /* candidate nuclide charge */
    unsigned int a; /* candidate nuclide mass number */
    double m;       /* nuclear mass from the isotope database [MeV/c^2] */

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

/* ---- Partition enumeration ------------------------------------------------ */

/* Candidate ground-state fragment species (the enumeration alphabet), built
 * once from the mass table and shared by the recursion. */
struct fbu_alphabet {
    uint8_t z[FBU_NDENSE];
    uint8_t a[FBU_NDENSE];
    int n;
};

/* Context for the recursive N-body partition enumeration. */
struct enum_ctx {
    double const *mass_mev;
    struct fbu_alphabet const *alpha;
    double m_parent;
    double c_parent; /* volume coefficient C for this parent */
    int n_total;     /* target N for this recursion subtree */
    uint8_t buf_z[FBU_NMAX];
    uint8_t buf_a[FBU_NMAX];
    /* output pointers (NULL during the size pass) */
    struct osh_fermi_partition *part_pool;
    struct osh_fermi_frag_spec *fspec_pool;
    size_t *part_total;
    size_t *fspec_total;
};

/* Build the enumeration alphabet: every (z, a) with a tabulated ground-state
 * mass, in canonical non-decreasing order (a, then z). */
static void build_alphabet(double const *mass_mev, struct fbu_alphabet *alpha) {
    unsigned int z; /* charge of the candidate species */
    unsigned int a; /* mass number of the candidate species */
    alpha->n = 0;
    for (a = 1u; a <= OSH_FERMI_BREAKUP_AMAX; ++a) {
        for (z = 0u; z <= OSH_FERMI_BREAKUP_ZMAX; ++z) {
            if (mass_mev[_za_idx(z, a)] > 0.0) {
                alpha->z[alpha->n] = (uint8_t) z;
                alpha->a[alpha->n] = (uint8_t) a;
                alpha->n += 1;
            }
        }
    }
}

/* Weight prefactor for the partition currently assembled in ctx->buf_[]:
 *   C^(N-1) * prod(g_i) / prod(n_k!) * prod(m_i^3/2) / Gamma(3N/2-3/2).
 * buf_[] is in non-decreasing alphabet order, so identical species are
 * adjacent and the identical-particle factorial is a product over run lengths. */
static double compute_prefactor(struct enum_ctx const *ctx, int n) {
    double prod_m32; /* running product of m_i^(3/2) over fragments */
    double prod_g;   /* running product of spin degeneracies (2J+1) */
    double id_denom; /* identical-particle factorial prod(n_k!) */
    double mi;       /* current fragment mass [MeV/c^2] */
    int i;           /* start index of the current identical-species run */
    int run;         /* length of the identical-species run at i */
    int f;           /* factorial counter while folding the run into id_denom */

    prod_m32 = 1.0;
    prod_g = 1.0;
    id_denom = 1.0;
    i = 0;
    while (i < n) {
        mi = ctx->mass_mev[_za_idx(ctx->buf_z[i], ctx->buf_a[i])];
        prod_m32 *= mi * sqrt(mi);
        prod_g *= (double) spin_degeneracy(ctx->buf_z[i], ctx->buf_a[i]);
        /* Count the run of identical species starting at i. */
        run = 1;
        while (i + run < n && ctx->buf_z[i + run] == ctx->buf_z[i] && ctx->buf_a[i + run] == ctx->buf_a[i]) {
            mi = ctx->mass_mev[_za_idx(ctx->buf_z[i + run], ctx->buf_a[i + run])];
            prod_m32 *= mi * sqrt(mi);
            prod_g *= (double) spin_degeneracy(ctx->buf_z[i + run], ctx->buf_a[i + run]);
            ++run;
        }
        for (f = run; f > 1; --f) {
            id_denom *= (double) f;
        }
        i += run;
    }

    return pow(ctx->c_parent, n - 1) * prod_g / id_denom * prod_m32 / s_gamma[n];
}

/* Recursive enumeration of N-body partitions over the full alphabet.
 * k_rem: number of fragments still to place.
 * z_rem, a_rem: remaining charge and mass to distribute.
 * min_sp: minimum alphabet index (canonical non-decreasing order).
 * buf_pos: next write position in ctx->buf_[]. */
static void enum_rec(struct enum_ctx *ctx, int k_rem, int z_rem, int a_rem, int min_sp, int buf_pos) {
    struct fbu_alphabet const *alpha = ctx->alpha; /* candidate species list */
    struct osh_fermi_partition *part;              /* partition slot written on the fill pass */
    int sp;                                        /* current alphabet index being tried */
    int zi;                                        /* charge of species sp */
    int ai;                                        /* mass number of species sp */
    int n;                                         /* fragment multiplicity of the completed partition */
    int i;                                         /* fragment loop index */
    double q;                                      /* Q-value = M_parent - sum(m_i) [MeV] */

    if (k_rem == 0) {
        if (z_rem != 0 || a_rem != 0) {
            return; /* charge/mass not exhausted: invalid */
        }
        n = ctx->n_total;
        if (ctx->part_pool != NULL) {
            part = &ctx->part_pool[*ctx->part_total];
            q = ctx->m_parent;
            for (i = 0; i < n; ++i) {
                q -= ctx->mass_mev[_za_idx(ctx->buf_z[i], ctx->buf_a[i])];
            }
            part->weight_prefactor = compute_prefactor(ctx, n);
            part->fspec_offset = (uint32_t) *ctx->fspec_total;
            part->q_mev = (float) q;
            part->n_frags = (uint8_t) n;
            part->_pad[0] = 0;
            part->_pad[1] = 0;
            part->_pad[2] = 0;
            for (i = 0; i < n; ++i) {
                ctx->fspec_pool[*ctx->fspec_total].z = ctx->buf_z[i];
                ctx->fspec_pool[*ctx->fspec_total].a = ctx->buf_a[i];
                (*ctx->fspec_total)++;
            }
        } else {
            *ctx->fspec_total += (size_t) n;
        }
        (*ctx->part_total)++;
        return;
    }

    for (sp = min_sp; sp < alpha->n; ++sp) {
        zi = alpha->z[sp];
        ai = alpha->a[sp];
        if (zi > z_rem || ai > a_rem) {
            continue;
        }
        /* Pruning: the remaining k_rem fragments each have at least mass ai
         * (species are non-decreasing in a, so ai is the current minimum). */
        if (ai * k_rem > a_rem) {
            continue;
        }
        ctx->buf_z[buf_pos] = (uint8_t) zi;
        ctx->buf_a[buf_pos] = (uint8_t) ai;
        enum_rec(ctx, k_rem - 1, z_rem - zi, a_rem - ai, sp, buf_pos + 1);
    }
}

/* Enumerate all N=2..FBU_NMAX partitions for every valid parent (z, a).
 * Pass pool==NULL, fspec==NULL for the sizing pass. */
static size_t enumerate_all_partitions(double const *mass_mev,
                                       struct fbu_alphabet const *alpha,
                                       struct osh_fermi_partition *part_pool,
                                       struct osh_fermi_frag_spec *fspec_pool,
                                       uint32_t *part_offset,
                                       uint16_t *part_count,
                                       size_t *nfspecs_out) {
    size_t part_total;  /* running count of partitions emitted so far */
    size_t fspec_total; /* running count of fragment-spec (z,a) entries */
    unsigned int z;     /* parent charge */
    unsigned int a;     /* parent mass number */
    double mp;          /* parent ground-state mass [MeV/c^2] */
    size_t part_start;  /* first partition index belonging to this parent */
    size_t idx;         /* dense (z,a) table index of the parent */
    int n;              /* target fragment multiplicity for this pass */
    struct enum_ctx ctx;

    part_total = 0u;
    fspec_total = 0u;
    for (z = 1u; z <= OSH_FERMI_BREAKUP_ZMAX; ++z) {
        for (a = 2u; a <= OSH_FERMI_BREAKUP_AMAX; ++a) {
            idx = _za_idx(z, a);
            mp = mass_mev[idx];
            part_start = part_total;

            if (mp > 0.0) {
                ctx.mass_mev = mass_mev;
                ctx.alpha = alpha;
                ctx.m_parent = mp;
                ctx.c_parent = volume_coeff(a);
                ctx.part_pool = part_pool;
                ctx.fspec_pool = fspec_pool;
                ctx.part_total = &part_total;
                ctx.fspec_total = &fspec_total;
                for (n = 2; n <= FBU_NMAX; ++n) {
                    if (n > (int) a) {
                        break; /* more fragments than nucleons */
                    }
                    ctx.n_total = n;
                    enum_rec(&ctx, n, (int) z, (int) a, 0, 0);
                }
            }

            if (part_offset != NULL) {
                part_offset[idx] = (uint32_t) part_start;
                part_count[idx] = (uint16_t) (part_total - part_start);
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
    int npow;    /* exponent 3K-5 of the chi^(3K-5)*(1-chi) pdf */
    double xn;   /* npow as a double */
    double fmax; /* pdf maximum, used as the rejection envelope */
    double chi;  /* trial sample in [0,1) */
    double f;    /* pdf value at chi */

    npow = (3 * k) - 5;
    xn = (double) npow;
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
static void kopylov_nbody(int n,
                          double const m[],
                          double m_eff,
                          double const p_parent[3],
                          struct osh_rng *rng,
                          double p_out[][3],
                          double e_tot[]) {
    double mu;            /* running remaining mass (unplaced fragments) */
    double mass;          /* current subsystem effective mass */
    double t;             /* available kinetic energy for current subsystem */
    double rest_mass;     /* effective mass of the "rest" pseudo-particle */
    double p_rest[3];     /* lab momentum of the rest pseudo-particle */
    double rest_m_cur;    /* rest mass of the rest pseudo-particle (for boost) */
    double p_cm[3];       /* fragment momentum in the current subsystem CM frame */
    double p_cm_mag;      /* magnitude of p_cm */
    double cos_theta;     /* polar cosine of the CM emission direction */
    double sin_theta;     /* polar sine of the CM emission direction */
    double cos_phi;       /* azimuth cosine */
    double sin_phi;       /* azimuth sine */
    double e_frag_cm;     /* fragment k total energy in the subsystem CM */
    double e_rest_cm;     /* rest pseudo-particle total energy in the subsystem CM */
    double e_frag_lab;    /* fragment k total energy in the lab */
    double e_rest_lab;    /* rest pseudo-particle total energy in the lab */
    double p_rest_new[3]; /* updated lab momentum of the rest pseudo-particle */
    double p2;            /* squared lab momentum of fragment 0 */
    int k;                /* fragment index, processed from n-1 down to 1 */

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

        cos_theta = (2.0 * osh_rng_double(rng)) - 1.0;
        sin_theta = sqrt(fmax(0.0, 1.0 - (cos_theta * cos_theta)));
        osh_kinematics_azimuth(rng, &cos_phi, &sin_phi);
        p_cm[0] = p_cm_mag * sin_theta * cos_phi;
        p_cm[1] = p_cm_mag * sin_theta * sin_phi;
        p_cm[2] = p_cm_mag * cos_theta;

        e_frag_cm = sqrt((m[k] * m[k]) + (p_cm_mag * p_cm_mag));
        osh_kinematics_boost_to_lab(rest_m_cur, p_rest, e_frag_cm, p_cm, &e_frag_lab, p_out[k]);
        e_tot[k] = e_frag_lab;

        /* Boost rest pseudo-particle (opposite direction) to lab. */
        p_cm[0] = -p_cm[0];
        p_cm[1] = -p_cm[1];
        p_cm[2] = -p_cm[2];
        e_rest_cm = sqrt((rest_mass * rest_mass) + (p_cm_mag * p_cm_mag));
        osh_kinematics_boost_to_lab(rest_m_cur, p_rest, e_rest_cm, p_cm, &e_rest_lab, p_rest_new);
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
    p2 = (p_rest[0] * p_rest[0]) + (p_rest[1] * p_rest[1]) + (p_rest[2] * p_rest[2]);
    e_tot[0] = sqrt((m[0] * m[0]) + p2);
}

/* ---- Weight selection helper --------------------------------------------- */

/* Compute total open-partition weight for parent at dense index idx and the
 * given excitation energy.  Returns 0 if no channel is open. */
static double partition_weight_sum(struct osh_nuclear_fermi_breakup const *model, size_t idx, double e_star) {
    double wsum;  /* accumulated weight of all open partitions */
    uint32_t off; /* first partition index for this parent */
    uint16_t cnt; /* number of partitions for this parent */
    uint16_t i;   /* partition loop index */
    double t;     /* released kinetic energy E*+Q of partition i [MeV] */
    int nf;       /* fragment multiplicity of partition i */

    wsum = 0.0;
    off = model->part_offset[idx];
    cnt = model->part_count[idx];
    for (i = 0u; i < cnt; ++i) {
        t = e_star + (double) model->part_pool[off + i].q_mev;
        if (t > 0.0) {
            nf = model->part_pool[off + i].n_frags;
            wsum += model->part_pool[off + i].weight_prefactor * pow(t, (1.5 * nf) - 2.5);
        }
    }
    return wsum;
}

/* ---- Public API ----------------------------------------------------------- */

enum osh_status osh_nuclear_fermi_breakup_compile(struct osh_nuclear_fermi_breakup *out) {
    static int const s_pdgs[FBU_SPECIES_COUNT] = {OSH_PART_PDG_NEUTRON,
                                                  OSH_PART_PDG_PROTON,
                                                  OSH_PART_PDG_DEUTERON,
                                                  OSH_PART_PDG_TRITON,
                                                  OSH_PART_PDG_HE3,
                                                  OSH_PART_PDG_HE4};
    struct fbu_alphabet *alpha; /* candidate species list (freed before return) */
    size_t nfspecs;             /* total fragment-spec entries (sizing pass) */
    size_t npart;               /* total partitions (sizing pass) */
    int i;                      /* species loop index */

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

    alpha = calloc(1u, sizeof(*alpha));
    if (!alpha) {
        osh_nuclear_fermi_breakup_free(out);
        return OSH_ENOMEM;
    }
    build_alphabet(out->mass_mev, alpha);

    /* Sizing pass. */
    npart = enumerate_all_partitions(out->mass_mev, alpha, NULL, NULL, NULL, NULL, &nfspecs);

    out->part_pool = calloc(npart > 0u ? npart : 1u, sizeof(*out->part_pool));
    out->fspec_pool = calloc(nfspecs > 0u ? nfspecs : 1u, sizeof(*out->fspec_pool));
    if (!out->part_pool || !out->fspec_pool) {
        free(alpha);
        osh_nuclear_fermi_breakup_free(out);
        return OSH_ENOMEM;
    }

    /* Fill pass. */
    out->npartitions = enumerate_all_partitions(
        out->mass_mev, alpha, out->part_pool, out->fspec_pool, out->part_offset, out->part_count, &out->nfspecs);

    free(alpha);
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

/* Emit one work-stack node as a transportable secondary (whitelist) or as an
 * unprocessed residual fragment.  Returns 1 if a transportable secondary was
 * emitted. */
static int emit_terminal(struct osh_nuclear_fermi_breakup const *model,
                         struct fbu_work_item const *node,
                         int room,
                         struct osh_nuclear_event *event_out) {
    struct osh_nuclear_secondary *sec; /* output secondary slot */
    int spec_idx;                      /* whitelist species index, or -1 if not transportable */
    double m_node;                     /* rest mass of this fragment [MeV/c^2] */
    double p2;                         /* squared lab momentum [MeV^2/c^2] */
    double p_norm;                     /* lab momentum magnitude [MeV/c] */

    spec_idx = final_species_index(node->z, node->a);
    if (!room || spec_idx < 0) {
        append_unprocessed_fragment(node, event_out);
        return 0;
    }

    m_node = model->mass_mev[_za_idx(node->z, node->a)];
    p2 = (node->p[0] * node->p[0]) + (node->p[1] * node->p[1]) + (node->p[2] * node->p[2]);
    p_norm = sqrt(p2);
    sec = &event_out->secondaries[event_out->n_secondaries];
    sec->energy = sqrt((m_node * m_node) + p2) - m_node;
    if (p_norm > 0.0) {
        sec->dir[0] = node->p[0] / p_norm;
        sec->dir[1] = node->p[1] / p_norm;
        sec->dir[2] = node->p[2] / p_norm;
    } else {
        sec->dir[0] = 0.0;
        sec->dir[1] = 0.0;
        sec->dir[2] = 1.0;
    }
    sec->species = &model->species[spec_idx];
    event_out->n_secondaries += 1u;
    return 1;
}

void osh_nuclear_fermi_breakup_step(struct osh_nuclear_fermi_breakup const *model,
                                    struct osh_nuclear_fragment const *fragment,
                                    struct osh_rng *rng,
                                    struct osh_nuclear_event *event_out) {
    struct fbu_work_item stack[FBU_WORK_STACK]; /* pending nuclides still to de-excite */
    struct fbu_work_item node;                  /* current item popped off the stack */
    struct fbu_work_item child;                 /* one decay product being pushed back */
    struct osh_nuclear_fragment frag;           /* local copy of the input prefragment */
    double m_arr[FBU_NMAX];                     /* selected partition's fragment masses [MeV/c^2] */
    double p_out[FBU_NMAX][3];                  /* per-fragment lab momenta from the decay [MeV/c] */
    double e_tot[FBU_NMAX];                     /* per-fragment lab energies (Kopylov scratch; unused here) */
    size_t idx;                                 /* dense (z,a) index of the current node */
    size_t n_emitted;                           /* number of transportable secondaries produced */
    int sp;                                     /* work-stack height (next free slot) */
    int room;                                   /* 1 while the secondaries array still has space */
    int fi;                                     /* fragment loop index */
    int nf;                                     /* fragment multiplicity of the selected partition */
    uint32_t off;                               /* first partition index for the current node */
    uint32_t fspec_off;                         /* fragment-spec offset of the selected partition */
    uint16_t cnt;                               /* partition count for the current node */
    uint16_t i;                                 /* partition loop index */
    uint16_t c;                                 /* index of the selected partition */
    double wsum;                                /* total weight of the current node's open partitions */
    double threshold;                           /* RNG draw scaled into [0,wsum) for selection */
    double cumulative;                          /* running weight sum during selection */
    double t;                                   /* released kinetic energy E*+Q of a partition [MeV] */
    double w;                                   /* weight of a single partition */
    double m_eff;                               /* effective parent mass m_gs + E* [MeV/c^2] */
    double p_star;                              /* two-body CM momentum magnitude [MeV/c] */
    double cos_theta;                           /* polar cosine of the two-body CM direction */
    double sin_theta;                           /* polar sine of the two-body CM direction */
    double cos_phi;                             /* azimuth cosine */
    double sin_phi;                             /* azimuth sine */
    double p_cm[3];                             /* product momentum in the parent rest frame [MeV/c] */
    double e_cm;                                /* product total energy in the parent rest frame [MeV] */
    double e_lab;                               /* product total lab energy [MeV] (discarded; momentum is kept) */

    frag = *fragment; /* local copy: fragment may alias event_out->fragments[0] */

    if (frag.a < 2u || frag.a > OSH_FERMI_BREAKUP_AMAX || frag.z < 1u || frag.z > OSH_FERMI_BREAKUP_ZMAX) {
        return;
    }
    idx = _za_idx(frag.z, frag.a);
    if (model->mass_mev[idx] <= 0.0) {
        return;
    }

    /* Gate: any open partition?  If not, leave the event untouched. */
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
    stack[0].z = (uint8_t) frag.z;
    stack[0].a = (uint8_t) frag.a;
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
            n_emitted += (size_t) emit_terminal(model, &node, room, event_out);
            continue;
        }

        /* Select one partition by cumulative weight. */
        off = model->part_offset[idx];
        cnt = model->part_count[idx];
        threshold = osh_rng_double(rng) * wsum;
        cumulative = 0.0;
        c = 0;
        for (i = 0u; i < cnt; ++i) {
            t = node.e_star + (double) model->part_pool[off + i].q_mev;
            if (t > 0.0) {
                nf = model->part_pool[off + i].n_frags;
                w = model->part_pool[off + i].weight_prefactor * pow(t, (1.5 * nf) - 2.5);
                cumulative += w;
                c = i;
                if (threshold <= cumulative) {
                    break;
                }
            }
        }

        nf = model->part_pool[off + c].n_frags;
        fspec_off = model->part_pool[off + c].fspec_offset;
        m_eff = model->mass_mev[idx] + node.e_star;

        /* Invariant: every compiled partition has 2..FBU_NMAX fragments.  Guard
         * it explicitly so a corrupt table can never index m_arr[]/p_out[] out
         * of bounds (and so the static analyzer can prove m_arr[] is fully
         * written before kopylov_nbody() reads m_arr[0]). */
        if (nf < 2 || nf > FBU_NMAX) {
            append_unprocessed_fragment(&node, event_out);
            continue;
        }

        /* Capacity guard: a partition can push up to nf products. */
        if (sp + nf > FBU_WORK_STACK) {
            append_unprocessed_fragment(&node, event_out);
            continue;
        }

        for (fi = 0; fi < nf; ++fi) {
            m_arr[fi] =
                model->mass_mev[_za_idx(model->fspec_pool[fspec_off + fi].z, model->fspec_pool[fspec_off + fi].a)];
        }

        if (nf == 2) {
            /* Two-body decay (exact). */
            p_star = osh_kinematics_two_body_decay_p(m_eff, m_arr[0], m_arr[1]);
            cos_theta = (2.0 * osh_rng_double(rng)) - 1.0;
            sin_theta = sqrt(fmax(0.0, 1.0 - (cos_theta * cos_theta)));
            osh_kinematics_azimuth(rng, &cos_phi, &sin_phi);
            p_cm[0] = p_star * sin_theta * cos_phi;
            p_cm[1] = p_star * sin_theta * sin_phi;
            p_cm[2] = p_star * cos_theta;

            e_cm = sqrt((m_arr[0] * m_arr[0]) + (p_star * p_star));
            osh_kinematics_boost_to_lab(m_eff, node.p, e_cm, p_cm, &e_lab, p_out[0]);
            p_cm[0] = -p_cm[0];
            p_cm[1] = -p_cm[1];
            p_cm[2] = -p_cm[2];
            e_cm = sqrt((m_arr[1] * m_arr[1]) + (p_star * p_star));
            osh_kinematics_boost_to_lab(m_eff, node.p, e_cm, p_cm, &e_lab, p_out[1]);
        } else {
            /* N>=3: Kopylov phase space. */
            kopylov_nbody(nf, m_arr, m_eff, node.p, rng, p_out, e_tot);
        }

        /* Push all products back onto the work stack (ground state). */
        for (fi = 0; fi < nf; ++fi) {
            child.p[0] = p_out[fi][0];
            child.p[1] = p_out[fi][1];
            child.p[2] = p_out[fi][2];
            child.e_star = 0.0;
            child.z = model->fspec_pool[fspec_off + fi].z;
            child.a = model->fspec_pool[fspec_off + fi].a;
            stack[sp] = child;
            ++sp;
        }
    }

    if (n_emitted > 0u) {
        event_out->kind = OSH_NUCLEAR_EVENT_FRAGMENTATION;
    }
}
