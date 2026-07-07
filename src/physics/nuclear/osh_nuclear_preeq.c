#include "physics/nuclear/osh_nuclear_preeq.h"

#include <math.h> /* sqrt, exp, log, cbrt, lgamma, ceil */
#include <string.h>

#include "openshieldhit/const.h"
#include "particle/osh_particle.h"
#include "particle/osh_particle_const.h"
#include "particle/osh_particle_pdg.h"
#include "physics/nuclear/osh_nuclear_handler.h"
#include "physics/osh_kinematics.h"
#include "random/osh_rng.h"

/* Ejectile species indices (order matches the species tables below). */
enum preeq_species {
    PREEQ_SPECIES_NEUTRON = 0,
    PREEQ_SPECIES_PROTON,
    PREEQ_SPECIES_DEUTERON,
    PREEQ_SPECIES_TRITON,
    PREEQ_SPECIES_HE3,
    PREEQ_SPECIES_HE4
};

static int const s_pdg[OSH_PREEQ_NSPECIES] = {
    OSH_PART_PDG_NEUTRON,
    OSH_PART_PDG_PROTON,
    OSH_PART_PDG_DEUTERON,
    OSH_PART_PDG_TRITON,
    OSH_PART_PDG_HE3,
    OSH_PART_PDG_HE4,
};

static unsigned int const s_z[OSH_PREEQ_NSPECIES] = {0u, 1u, 1u, 1u, 2u, 2u};
static unsigned int const s_a[OSH_PREEQ_NSPECIES] = {1u, 1u, 2u, 3u, 3u, 4u};

/* Spin multiplicity 2s+1 of the ejectile ground state. */
static double const s_spin_mult[OSH_PREEQ_NSPECIES] = {2.0, 2.0, 3.0, 2.0, 2.0, 1.0};

/* Binomial C(a_j, z_j) of the ejectile composition (for the R_j factor). */
static double const s_comp_binom[OSH_PREEQ_NSPECIES] = {1.0, 1.0, 2.0, 3.0, 3.0, 6.0};

/* Emission energy grid: bins of at most ~1 MeV, capped for stack safety. */
#define PREEQ_GRID_MAX 512

/* Hard bound on emission + transition steps per fragment; physically the
 * loop ends at n_eq = sqrt(2 g E*) <~ 25 for any residue this table holds. */
#define PREEQ_MAX_ITER 128

/* Dense (z, a) index into the mass table. */
static inline size_t _za_idx(unsigned int z, unsigned int a) {
    return ((size_t) z * (size_t) (OSH_PREEQ_AMAX + 1)) + (size_t) a;
}

/* Condensation probability gamma_j for cluster ejectiles (Geant4 closed
 * forms: 16/A, 243/A^2, 4096/A^3), times the per-species calibration scale;
 * 1 for nucleons. */
static double condensation_gamma(int spec, double a_parent) {
    if (spec == PREEQ_SPECIES_DEUTERON) {
        return OSH_PREEQ_GAMMA_SCALE_D * 16.0 / a_parent;
    }
    if (spec == PREEQ_SPECIES_TRITON) {
        return OSH_PREEQ_GAMMA_SCALE_T * 243.0 / (a_parent * a_parent);
    }
    if (spec == PREEQ_SPECIES_HE3) {
        return OSH_PREEQ_GAMMA_SCALE_HE3 * 243.0 / (a_parent * a_parent);
    }
    if (spec == PREEQ_SPECIES_HE4) {
        return OSH_PREEQ_GAMMA_SCALE_ALPHA * 4096.0 / (a_parent * a_parent * a_parent);
    }
    return 1.0;
}

/* Inverse cross section [fm^2], Dostrovsky-style geometric estimate.
 * Neutrons: pi R^2 alpha (1 + beta / T); charged: pi R^2 (1 - V_C / T),
 * clamped at zero below the barrier. */
static double inverse_sigma_fm2(int spec, double t_mev, double a_res, double z_res) {
    double r_fm;
    double sigma_geom;
    double alpha;
    double beta;
    double v_coul;

    r_fm = OSH_PREEQ_EMISSION_R0_FM * (cbrt(a_res) + cbrt((double) s_a[spec]));
    sigma_geom = OSH_M_PI * r_fm * r_fm;

    if (spec == PREEQ_SPECIES_NEUTRON) {
        alpha = 0.76 + 2.2 / cbrt(a_res);
        beta = (2.12 / (cbrt(a_res) * cbrt(a_res)) - 0.05) / alpha;
        return sigma_geom * alpha * (1.0 + beta / t_mev);
    }

    v_coul = OSH_E2_MEV_FM * (double) s_z[spec] * z_res / r_fm;
    if (t_mev <= v_coul) {
        return 0.0;
    }
    return sigma_geom * (1.0 - v_coul / t_mev);
}

/*
 * Emission width of species spec from configuration (p, h, E*) of the parent
 * (z, a), integrated over the ejectile energy grid [Betak-Dobes detailed
 * balance with Williams state densities; see the header].  Returns the total
 * width [1/fm].  When grid_out is non-NULL it receives the per-bin
 * contributions (nbins_out entries of width dt_out), for CDF sampling.
 */
static double emission_width(struct osh_nuclear_preeq const *model,
                             int spec,
                             unsigned int z,
                             unsigned int a,
                             double e_star,
                             unsigned int exc_p,
                             unsigned int exc_h,
                             double *grid_out,
                             int *nbins_out,
                             double *dt_out,
                             double *b_sep_out) {
    double m_parent;
    double m_res;
    double m_ej;
    double b_sep;
    double t_max;
    double g_level;
    double mu_red;
    double prefactor;
    double gamma_j;
    double r_comp;
    double comb;
    double log_ge;
    double f_p;
    double f_n;
    double total;
    double t_mid;
    double u_res;
    double ratio;
    double w;
    double dt;
    unsigned int z_res;
    unsigned int a_res;
    unsigned int n_exc;
    int nbins;
    int i;
    int n_minus;

    if (grid_out != NULL) {
        *nbins_out = 0;
        *dt_out = 0.0;
    }
    if (b_sep_out != NULL) {
        *b_sep_out = 0.0;
    }

    /* Channel-closure guards: enough particle excitons of any kind, a real
     * residual nucleus left behind (A_res >= 2, so no bare-nucleon residues
     * enter the fragment pool), and both nuclides in the mass table. */
    n_exc = exc_p + exc_h;
    if (exc_p < s_a[spec] || n_exc < s_a[spec] + 1u) {
        return 0.0;
    }
    if (a < s_a[spec] + 2u || z < s_z[spec] || (a - s_a[spec]) < (z - s_z[spec])) {
        return 0.0;
    }
    z_res = z - s_z[spec];
    a_res = a - s_a[spec];
    if (z > OSH_PREEQ_ZMAX || a > OSH_PREEQ_AMAX) {
        return 0.0;
    }
    m_parent = model->mass_mev[_za_idx(z, a)];
    m_res = model->mass_mev[_za_idx(z_res, a_res)];
    m_ej = model->species[spec].mass;
    if (m_parent <= 0.0 || m_res <= 0.0) {
        return 0.0;
    }

    b_sep = m_res + m_ej - m_parent;
    if (b_sep_out != NULL) {
        *b_sep_out = b_sep;
    }
    t_max = e_star - b_sep;
    if (t_max <= 0.0) {
        return 0.0;
    }

    /* Williams-ratio combinatorics (independent of T):
     *   omega(p - a_j, h, U) / omega(p, h, E) =
     *     p!/(p-a_j)! * (n-1)!/(n-a_j-1)! * (gU)^(n-a_j-1) / (gE)^(n-1) */
    comb = exp(lgamma((double) exc_p + 1.0) - lgamma((double) (exc_p - s_a[spec]) + 1.0) + lgamma((double) n_exc)
               - lgamma((double) (n_exc - s_a[spec])));

    g_level = OSH_PREEQ_LEVEL_DENSITY_PER_A * (double) a;
    log_ge = log(g_level * e_star);
    n_minus = (int) n_exc - (int) s_a[spec] - 1;

    mu_red = m_ej * m_res / (m_ej + m_res);
    prefactor = s_spin_mult[spec] * mu_red / (OSH_M_PI * OSH_M_PI * OSH_HBARC * OSH_HBARC * OSH_HBARC);
    gamma_j = condensation_gamma(spec, (double) a);

    f_p = (double) z / (double) a;
    f_n = ((double) a - (double) z) / (double) a;
    r_comp = s_comp_binom[spec] * pow(f_p, (double) s_z[spec]) * pow(f_n, (double) (s_a[spec] - s_z[spec]));

    nbins = (int) ceil(t_max);
    if (nbins < 1) {
        nbins = 1;
    }
    if (nbins > PREEQ_GRID_MAX) {
        nbins = PREEQ_GRID_MAX;
    }
    dt = t_max / (double) nbins;

    total = 0.0;
    for (i = 0; i < nbins; ++i) {
        t_mid = ((double) i + 0.5) * dt;
        u_res = t_max - t_mid; /* = E* - B_j - T, residual excitation */
        ratio = comb * exp(((double) n_minus * log(g_level * u_res)) - (((double) n_exc - 1.0) * log_ge));
        w = prefactor * t_mid * inverse_sigma_fm2(spec, t_mid, (double) a_res, (double) z_res) * gamma_j * r_comp
            * ratio * dt;
        if (w < 0.0) {
            w = 0.0;
        }
        if (grid_out != NULL) {
            grid_out[i] = w;
        }
        total += w;
    }
    if (grid_out != NULL) {
        *nbins_out = nbins;
        *dt_out = dt;
    }
    return total;
}

/*
 * Intranuclear transition rate lambda_plus (Delta-n = +2) [1/fm], CEM form
 * (Gudima-Mashnik-Toneev; same expressions as the Geant4
 * G4PreCompoundTransitions CEM branch): averaged in-medium NN cross section
 * at the mean exciton relative energy, Pauli-blocking factor, interaction
 * volume from the transitions radius.
 */
static double lambda_plus(unsigned int z, unsigned int a, double e_star, unsigned int n_exc) {
    double e_rel;
    double v_rel;
    double v_sq;
    double sigma_pp;
    double sigma_np;
    double sigma_avg;
    double f_p;
    double f_n;
    double pauli;
    double pauli_x;
    double ratio;
    double xx;
    double v_int;
    double za;
    double na;

    if (n_exc == 0u || a < 2u) {
        return 0.0;
    }

    e_rel = 1.6 * OSH_PREEQ_FERMI_ENERGY_MEV + e_star / (double) n_exc;
    v_sq = 2.0 * e_rel / OSH_PART_MASS_PROTON;
    v_rel = sqrt(v_sq);

    /* In-medium NN cross sections [fm^2]; parametrisation in beta = v/c
     * (1 mb = 0.1 fm^2). */
    sigma_pp = (10.63 / v_sq - 29.92 / v_rel + 42.9) * 0.1;
    sigma_np = (34.10 / v_sq - 82.20 / v_rel + 82.2) * 0.1;

    /* Average over the (untracked) charge of the fast exciton with the
     * residue proton/neutron fractions as weights. */
    za = (double) z;
    na = (double) (a - z);
    f_p = za / (double) a;
    f_n = na / (double) a;
    sigma_avg = f_p * ((za - 1.0) * sigma_pp + na * sigma_np) / ((double) a - 1.0)
                + f_n * ((na - 1.0) * sigma_pp + za * sigma_np) / ((double) a - 1.0);

    ratio = OSH_PREEQ_FERMI_ENERGY_MEV / e_rel;
    pauli = 1.0 - 1.4 * ratio;
    if (ratio > 0.5) {
        pauli_x = 2.0 - 1.0 / ratio;
        pauli += 0.4 * ratio * pauli_x * pauli_x * sqrt(pauli_x);
    }

    xx = 2.0 * OSH_PREEQ_TRANSITIONS_R0_FM + OSH_HBARC / (OSH_PART_MASS_PROTON * v_rel);
    v_int = (4.0 / 3.0) * OSH_M_PI * xx * xx * xx;

    if (sigma_avg < 0.0) {
        return 0.0;
    }
    if (pauli < 0.0) {
        return 0.0;
    }
    return sigma_avg * pauli * v_rel / v_int;
}

enum osh_status osh_nuclear_preeq_compile(struct osh_nuclear_preeq *model) {
    unsigned int z;
    unsigned int a;
    double m;
    int i;

    if (model == NULL) {
        return OSH_EINVAL;
    }
    memset(model, 0, sizeof(*model));

    model->mass_mev[_za_idx(0u, 1u)] = OSH_PART_MASS_NEUTRON;
    for (z = 1u; z <= OSH_PREEQ_ZMAX; ++z) {
        for (a = z; a <= OSH_PREEQ_AMAX; ++a) {
            if (osh_particle_nuclear_mass_mev_from_za(z, a, &m)) {
                model->mass_mev[_za_idx(z, a)] = m;
            }
        }
    }

    for (i = 0; i < OSH_PREEQ_NSPECIES; ++i) {
        if (!osh_particle_from_pdg(&model->species[i], s_pdg[i])) {
            return OSH_ESTATE;
        }
    }

    model->compiled = 1;
    return OSH_OK;
}

void osh_nuclear_preeq_step(struct osh_nuclear_preeq const *model,
                            struct osh_nuclear_fragment *fragment,
                            struct osh_rng *rng,
                            struct osh_nuclear_event *event_out) {
    double grid[PREEQ_GRID_MAX];
    double width[OSH_PREEQ_NSPECIES];
    double b_sep[OSH_PREEQ_NSPECIES];
    double width_total;
    double lp;
    double g_level;
    double n_eq;
    double u;
    double cum;
    double t_emit;
    double dt;
    double p_ej;
    double cos_theta;
    double sin_theta;
    double cos_phi;
    double sin_phi;
    double dir[3];
    unsigned int n_exc;
    int nbins;
    int iter;
    int spec;
    int i;
    size_t slot;

    if (model == NULL || !model->compiled || fragment == NULL || rng == NULL || event_out == NULL) {
        return;
    }

    for (iter = 0; iter < PREEQ_MAX_ITER; ++iter) {
        /* Thermalized, empty, or out-of-domain residues are left untouched:
         * the compound/neutron path and break-up residues carry (0, 0). */
        if (fragment->excitons_p == 0u || fragment->excitation_energy <= 0.0) {
            return;
        }
        if (fragment->z > OSH_PREEQ_ZMAX || fragment->a > OSH_PREEQ_AMAX || fragment->a < 2u) {
            return;
        }

        /* Equilibrium exit: n_eq = sqrt(2 g E*) (CEM03.03 Eq. 34). */
        n_exc = fragment->excitons_p + fragment->excitons_h;
        g_level = OSH_PREEQ_LEVEL_DENSITY_PER_A * (double) fragment->a;
        n_eq = sqrt(2.0 * g_level * fragment->excitation_energy);
        if ((double) n_exc >= n_eq) {
            return;
        }

        width_total = 0.0;
        for (spec = 0; spec < OSH_PREEQ_NSPECIES; ++spec) {
            width[spec] = emission_width(model,
                                         spec,
                                         fragment->z,
                                         fragment->a,
                                         fragment->excitation_energy,
                                         fragment->excitons_p,
                                         fragment->excitons_h,
                                         NULL,
                                         NULL,
                                         NULL,
                                         &b_sep[spec]);
            width_total += width[spec];
        }
        lp = lambda_plus(fragment->z, fragment->a, fragment->excitation_energy, n_exc);

        if (width_total + lp <= 0.0) {
            return;
        }

        u = osh_rng_double(rng) * (width_total + lp);
        if (u < lp) {
            /* Intranuclear collision: one more particle-hole pair; the
             * configuration random-walks toward equilibrium. */
            fragment->excitons_p += 1u;
            fragment->excitons_h += 1u;
            continue;
        }

        /* Emission: pick the species proportional to its width. */
        u -= lp;
        cum = 0.0;
        spec = OSH_PREEQ_NSPECIES - 1;
        for (i = 0; i < OSH_PREEQ_NSPECIES; ++i) {
            cum += width[i];
            if (u <= cum) {
                spec = i;
                break;
            }
        }

        if (event_out->n_secondaries >= OSH_NUCLEAR_MAX_SECONDARIES) {
            return; /* out of slots: leave the rest to equilibrium break-up */
        }

        /* Sample the ejectile energy from the per-bin width grid. */
        emission_width(model,
                       spec,
                       fragment->z,
                       fragment->a,
                       fragment->excitation_energy,
                       fragment->excitons_p,
                       fragment->excitons_h,
                       grid,
                       &nbins,
                       &dt,
                       &b_sep[spec]);
        if (nbins < 1) {
            return;
        }
        /* CDF walk over the bins; falls through to the last bin on float
         * round-off. */
        u = osh_rng_double(rng) * width[spec];
        cum = 0.0;
        for (i = 0; i < nbins - 1; ++i) {
            cum += grid[i];
            if (u <= cum) {
                break;
            }
        }
        t_emit = ((double) i + osh_rng_double(rng)) * dt;

        /* Isotropic lab emission (Kalbach systematics: planned refinement). */
        cos_theta = 2.0 * osh_rng_double(rng) - 1.0;
        sin_theta = sqrt(fmax(0.0, 1.0 - (cos_theta * cos_theta)));
        osh_kinematics_azimuth(rng, &cos_phi, &sin_phi);
        dir[0] = sin_theta * cos_phi;
        dir[1] = sin_theta * sin_phi;
        dir[2] = cos_theta;

        slot = event_out->n_secondaries;
        event_out->secondaries[slot].dir[0] = dir[0];
        event_out->secondaries[slot].dir[1] = dir[1];
        event_out->secondaries[slot].dir[2] = dir[2];
        event_out->secondaries[slot].energy = t_emit;
        event_out->secondaries[slot].species = &model->species[spec];
        event_out->n_secondaries = slot + 1u;

        /* Momentum balance: the ejectile momentum leaves the residue. */
        p_ej = sqrt(t_emit * (t_emit + 2.0 * model->species[spec].mass));
        fragment->p[0] -= p_ej * dir[0];
        fragment->p[1] -= p_ej * dir[1];
        fragment->p[2] -= p_ej * dir[2];

        /* Exact bookkeeping: E* loses the separation energy and the
         * ejectile kinetic energy; the residue drops a_j particle excitons
         * and (z_j, a_j) nucleons. */
        fragment->excitation_energy -= b_sep[spec] + t_emit;
        if (fragment->excitation_energy < 0.0) {
            fragment->excitation_energy = 0.0;
        }
        fragment->a -= s_a[spec];
        fragment->z -= s_z[spec];
        fragment->excitons_p -= s_a[spec];

        event_out->kind = OSH_NUCLEAR_EVENT_FRAGMENTATION;
    }
}
