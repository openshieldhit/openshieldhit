#include "physics/nuclear/osh_nuclear_abrasion.h"

#include <math.h> /* sqrt, exp, cbrt */

#include "openshieldhit/const.h"
#include "particle/osh_particle.h"
#include "particle/osh_particle_const.h"
#include "particle/osh_particle_pdg.h"
#include "physics/nuclear/osh_nuclear_handler.h"
#include "physics/osh_kinematics.h"
#include "random/osh_rng.h"

/* Compile-time constants; values match osh_particle_from_pdg for these PDG codes. */
static struct particle const s_proton = {OSH_PART_MASS_PROTON, OSH_PART_PDG_PROTON, 1, 1u, 1u, 0u};
static struct particle const s_neutron = {OSH_PART_MASS_NEUTRON, OSH_PART_PDG_NEUTRON, 0, 0u, 1u, 0u};

/**
 * Retention probability for a knocked-out nucleon of lab kinetic energy
 * e_mev (issue #263): a smooth Fermi-function turn-on below the threshold —
 * the lean analog of INCL4.6's "back to spectator" re-absorption.  The
 * threshold for protons includes a Coulomb-barrier term supplied by the
 * caller; width 0 degrades gracefully to a hard step.
 */
static double retention_probability(double e_mev, double e_thr_mev) {
    double const w = OSH_ABRASION_RETENTION_WIDTH_MEV;

    if (e_thr_mev <= 0.0) {
        return 0.0;
    }
    if (w <= 0.0) {
        if (e_mev <= e_thr_mev) {
            return 1.0;
        }
        return 0.0;
    }
    return 1.0 / (1.0 + exp((e_mev - e_thr_mev) / w));
}

/**
 * Zero-truncated Poisson sample (nu >= 1) by CDF inversion.
 *
 * The reaction cross-section has already fired when this is called, so the
 * event must wound at least one nucleon; sampling nu = 0 would produce a
 * null event that silently weakens the effective inelastic cross-section.
 * The loop is bounded by OSH_NUCLEAR_MAX_SECONDARIES.
 */
static int sample_nu_truncated(struct osh_rng *rng, double lambda) {
    double u;
    double term;
    double cum;
    int k;

    if (lambda <= 0.0) {
        return 1;
    }

    /* Map a uniform deviate into the truncated CDF range (P(0), 1]. */
    term = exp(-lambda);
    u = term + osh_rng_double(rng) * (1.0 - term);
    cum = term;
    k = 0;
    while (cum < u && k < OSH_NUCLEAR_MAX_SECONDARIES) {
        ++k;
        term *= lambda / (double) k;
        cum += term;
    }
    if (k < 1) {
        return 1;
    }
    return k;
}

void osh_nuclear_abrasion_step(double T_lab_mev,
                               double const incident_dir[3],
                               double a_eff,
                               double z_eff,
                               double sigma_pa_cm2,
                               struct osh_rng *rng,
                               struct osh_nuclear_event *event_out) {
    double sigma_pn_cm2;
    double nu_mean;
    double n_over_a;
    double T_current;
    double cos_cm;
    double cos_phi;
    double sin_phi;
    double cos_prim; /* deflection cos of the continuing cascade proton after collision j */
    double sin_prim;
    double e_prim;  /* KE of the continuing cascade proton [MeV] */
    double cos_sec; /* deflection cos of the knocked-out nucleon */
    double e_sec;   /* lab KE of the knocked-out nucleon [MeV] */
    double sin_sec;
    double p_in;        /* incident lab momentum magnitude [MeV/c] */
    double p_sec;       /* emitted nucleon lab momentum magnitude [MeV/c] */
    double w_curr[3];   /* current cascade-proton direction (composed per collision) */
    double e_star;      /* prefragment excitation energy: hole charges + retained KE [MeV] */
    double coulomb_mev; /* ~2/3 of the residue Coulomb barrier, added to the proton retention threshold */
    double e_thr;       /* retention threshold for the current knockout species [MeV] */
    int is_neutron;
    int n_knockout_p;
    int n_knockout_n;
    unsigned int n_excitons_p;
    unsigned int n_excitons_h;
    int nu; /* sampled number of participant nucleons */
    int j;
    size_t s;
    struct particle const *species;

    /* The primary slot is always terminated in an inelastic event; the
     * degraded cascade proton continues as a secondary when it escapes. */
    event_out->kind = OSH_NUCLEAR_EVENT_ABRASION;
    event_out->primary_energy = 0.0;
    event_out->n_secondaries = 0u;
    event_out->n_fragments = 0u;

    if (sigma_pa_cm2 <= 0.0 || a_eff <= 0.0) {
        return;
    }

    /* The residual target prefragment is handed to the Fermi break-up stage.
     * Fast nucleon emission below will update A/Z, E*, and momentum before
     * return.  Momentum balance: the prefragment starts with the full
     * incident momentum and the momentum of every escaping particle
     * (knock-out nucleons and the cascade proton) is subtracted; an absorbed
     * cascade proton's momentum thus stays with the fragment. */
    event_out->n_fragments = 1u;
    event_out->fragments[0].a = (unsigned int) floor(a_eff + 0.5);
    event_out->fragments[0].z = (unsigned int) floor(z_eff + 0.5);
    event_out->fragments[0].excitation_energy = 0.0;
    event_out->fragments[0].excitons_p = 0u;
    event_out->fragments[0].excitons_h = 0u;
    p_in = sqrt(T_lab_mev * (T_lab_mev + 2.0 * OSH_PART_MASS_PROTON));
    event_out->fragments[0].p[0] = p_in * incident_dir[0];
    event_out->fragments[0].p[1] = p_in * incident_dir[1];
    event_out->fragments[0].p[2] = p_in * incident_dir[2];

    /* Wounded-nucleon mean (Glauber 1970):
     *   <nu> = sigma_pN * A / sigma_pA
     * sigma_pN is a tunable constant (OSH_ABRASION_SIGMA_PN_MB); calibrate
     * against public secondary-nucleon yield benchmark data.  nu is sampled
     * zero-truncated: the reaction already fired, so at least one nucleon
     * participates. */
    sigma_pn_cm2 = OSH_ABRASION_SIGMA_PN_MB * OSH_MB_TO_CM2;
    nu_mean = sigma_pn_cm2 * a_eff / sigma_pa_cm2;
    nu = sample_nu_truncated(rng, nu_mean);
    if (nu > OSH_NUCLEAR_MAX_SECONDARIES - 1) {
        nu = OSH_NUCLEAR_MAX_SECONDARIES - 1; /* reserve one slot for the cascade proton */
    }
    if (nu > (int) a_eff) {
        nu = (int) a_eff; /* cannot wound more nucleons than the target holds */
    }

    /* N/A = neutron fraction of the target material.
     * Each participant is independently drawn as n (prob N/A) or p (prob Z/A).
     * Clamp handles edge cases where z_eff or a_eff come from a pure-hydrogen
     * material (N/A = 0) or a hypothetical neutron-only material. */
    n_over_a = (a_eff - z_eff) / a_eff;
    if (n_over_a < 0.0) {
        n_over_a = 0.0;
    }
    if (n_over_a > 1.0) {
        n_over_a = 1.0;
    }

    /* Intranuclear-cascade picture: the proton undergoes nu quasi-elastic
     * collisions, each knocking out one nucleon and deflecting the continuing
     * proton (directions composed collision by collision).  Each abraded
     * nucleon leaves a hole; its excitation cost
     * (OSH_ABRASION_EXCITATION_PER_HOLE_MEV) is charged to the cascade proton
     * and booked as prefragment excitation energy.  A knockout near or below
     * the retention threshold is re-absorbed (issue #263): its kinetic energy
     * funds E* and it becomes a particle exciton.  The event conserves
     * kinetic energy exactly:  T_in = sum(KE_escaped) + E*. */
    T_current = T_lab_mev;
    e_star = 0.0;
    n_knockout_p = 0;
    n_knockout_n = 0;
    n_excitons_p = 0u;
    n_excitons_h = 0u;
    w_curr[0] = incident_dir[0];
    w_curr[1] = incident_dir[1];
    w_curr[2] = incident_dir[2];

    /* Proton retention threshold includes ~2/3 of the residue Coulomb barrier
     * (INCL4.6 recipe); V_C = e^2 Z / (r0 A^(1/3)), evaluated once for the
     * initial target — the per-event change of the residue is beyond this
     * model's accuracy. */
    coulomb_mev = (2.0 / 3.0) * OSH_E2_MEV_FM * z_eff / (OSH_ABRASION_COULOMB_R0_FM * cbrt(a_eff));

    for (j = 0; j < nu; ++j) {
        /* Stop if the remaining proton energy is too low for meaningful
         * kinematics (below ~1 MeV the equal-mass formula loses accuracy
         * and the nucleon would range out immediately anyway). */
        if (T_current < 1.0) {
            break;
        }

        /* Nucleon type: neutron with prob N/A, proton with prob Z/A. */
        is_neutron = (osh_rng_double(rng) < n_over_a);
        if (is_neutron) {
            species = &s_neutron;
        } else {
            species = &s_proton;
        }

        /* Isotropic CM scattering — simplest physical model.
         * A forward-peaked angular distribution (e.g. from a diffraction
         * model) would be a refinement for future work. */
        cos_cm = 2.0 * osh_rng_double(rng) - 1.0;
        osh_kinematics_azimuth(rng, &cos_phi, &sin_phi);

        /* Full relativistic equal-mass kinematics (proton mass used for both
         * p and n since mp ≈ mn to better than 0.1%). */
        osh_kinematics_elastic_equal_mass_lab(
            T_current, OSH_PART_MASS_PROTON, cos_cm, &cos_prim, &e_prim, &cos_sec, &e_sec);

        /* Every collision vacates one orbital: a hole exciton. */
        ++n_excitons_h;

        e_thr = OSH_ABRASION_RETENTION_THRESHOLD_MEV;
        if (e_thr > 0.0 && !is_neutron) {
            e_thr += coulomb_mev;
        }
        if (osh_rng_double(rng) < retention_probability(e_sec, e_thr)) {
            /* Retained: the struck nucleon stays in the residue as a particle
             * exciton; its kinetic energy funds the excitation.  A/Z and the
             * fragment momentum are left untouched. */
            e_star += e_sec;
            ++n_excitons_p;
        } else {
            if (is_neutron) {
                ++n_knockout_n;
            } else {
                ++n_knockout_p;
            }

            /* Rotate the knocked-out nucleon from the current cascade
             * direction.  Convention matches pp elastic: the recoil uses the
             * opposite azimuth (-cos_phi, -sin_phi) so that momentum is
             * conserved transversely. */
            s = event_out->n_secondaries;
            sin_sec = sqrt(fmax(0.0, 1.0 - (cos_sec * cos_sec)));
            osh_kinematics_rotate_dir_cos(w_curr, event_out->secondaries[s].dir, cos_sec, sin_sec, -cos_phi, -sin_phi);

            event_out->secondaries[s].energy = e_sec;
            event_out->secondaries[s].species = species;
            event_out->n_secondaries = s + 1u;

            /* Momentum balance: subtract the emitted nucleon from the fragment. */
            p_sec = sqrt(e_sec * (e_sec + 2.0 * OSH_PART_MASS_PROTON));
            event_out->fragments[0].p[0] -= p_sec * event_out->secondaries[s].dir[0];
            event_out->fragments[0].p[1] -= p_sec * event_out->secondaries[s].dir[1];
            event_out->fragments[0].p[2] -= p_sec * event_out->secondaries[s].dir[2];
        }

        /* Deflect the continuing cascade proton. */
        sin_prim = sqrt(fmax(0.0, 1.0 - (cos_prim * cos_prim)));
        osh_kinematics_rotate_dir_cos(w_curr, w_curr, cos_prim, sin_prim, cos_phi, sin_phi);

        /* Charge the hole excitation cost to the cascade proton; if it cannot
         * pay, it is absorbed and its remaining energy funds the excitation. */
        if (e_prim > OSH_ABRASION_EXCITATION_PER_HOLE_MEV) {
            e_star += OSH_ABRASION_EXCITATION_PER_HOLE_MEV;
            T_current = e_prim - OSH_ABRASION_EXCITATION_PER_HOLE_MEV;
        } else {
            e_star += e_prim;
            T_current = 0.0;
            break;
        }
    }

    /* The degraded cascade proton escapes the nucleus and continues as a
     * transportable secondary (intranuclear-cascade picture).  Below ~1 MeV
     * it is absorbed instead and its energy adds to the excitation. */
    if (T_current >= 1.0 && event_out->n_secondaries < OSH_NUCLEAR_MAX_SECONDARIES) {
        s = event_out->n_secondaries;
        event_out->secondaries[s].dir[0] = w_curr[0];
        event_out->secondaries[s].dir[1] = w_curr[1];
        event_out->secondaries[s].dir[2] = w_curr[2];
        event_out->secondaries[s].energy = T_current;
        event_out->secondaries[s].species = &s_proton;
        event_out->n_secondaries = s + 1u;

        p_sec = sqrt(T_current * (T_current + 2.0 * OSH_PART_MASS_PROTON));
        event_out->fragments[0].p[0] -= p_sec * w_curr[0];
        event_out->fragments[0].p[1] -= p_sec * w_curr[1];
        event_out->fragments[0].p[2] -= p_sec * w_curr[2];
    } else {
        /* Absorbed cascade proton: its remaining energy funds the excitation
         * and it joins the residue as a particle exciton (it entered from
         * outside, so no hole is associated with it). */
        e_star += T_current;
        ++n_excitons_p;
    }

    if (event_out->n_fragments > 0u) {
        unsigned int da = (unsigned int) (n_knockout_p + n_knockout_n);
        unsigned int dz = (unsigned int) n_knockout_p;

        if (event_out->fragments[0].a > da) {
            event_out->fragments[0].a -= da;
        } else {
            event_out->fragments[0].a = 0u;
        }
        if (event_out->fragments[0].z > dz) {
            event_out->fragments[0].z -= dz;
        } else {
            event_out->fragments[0].z = 0u;
        }
        if (event_out->fragments[0].a == 0u || event_out->fragments[0].z == 0u) {
            /* Fully disassembled residual (light targets): any accumulated
             * excitation is dropped with it — rare edge, accepted. */
            event_out->n_fragments = 0u;
        }
    }

    if (event_out->n_fragments > 0u) {
        event_out->fragments[0].excitation_energy = e_star;
        event_out->fragments[0].excitons_p = n_excitons_p;
        event_out->fragments[0].excitons_h = n_excitons_h;
    }
}
