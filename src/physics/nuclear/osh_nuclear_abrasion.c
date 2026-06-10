#include "physics/nuclear/osh_nuclear_abrasion.h"

#include <math.h> /* sqrt */

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
    double cos_prim; /* lab cos of the continuing (degraded) proton after collision j */
    double e_prim;   /* lab KE of the continuing proton [MeV] */
    double cos_sec;  /* lab cos of the knocked-out nucleon */
    double e_sec;    /* lab KE of the knocked-out nucleon [MeV] */
    double sin_sec;
    int is_neutron;
    int n_knockout_p;
    int n_knockout_n;
    int nu; /* sampled number of participant nucleons */
    int j;
    struct particle const *species;

    /* Primary is always absorbed in an inelastic event.  Secondaries are filled
     * in below; n_secondaries stays 0 if nu=0 (Poisson fluctuation). */
    event_out->kind = OSH_NUCLEAR_EVENT_ABRASION;
    event_out->primary_energy = 0.0;
    event_out->n_secondaries = 0u;
    event_out->n_fragments = 0u;

    if (sigma_pa_cm2 <= 0.0 || a_eff <= 0.0) {
        return;
    }

    /* The residual target nucleus is stored for the future Fermi-breakup
     * stage.  Fast nucleon emission below will update A/Z before return. */
    event_out->n_fragments = 1u;
    event_out->fragments[0].a = (unsigned int) floor(a_eff + 0.5);
    event_out->fragments[0].z = (unsigned int) floor(z_eff + 0.5);
    event_out->fragments[0].excitation_energy = 0.0;

    /* Wounded-nucleon mean (Glauber 1970):
     *   <nu> = sigma_pN * A / sigma_pA
     * sigma_pN is a tunable constant (OSH_ABRASION_SIGMA_PN_MB); calibrate
     * against public secondary-nucleon yield benchmark data. */
    sigma_pn_cm2 = OSH_ABRASION_SIGMA_PN_MB * OSH_MB_TO_CM2;
    nu_mean = sigma_pn_cm2 * a_eff / sigma_pa_cm2;
    nu = osh_rng_poisson(rng, nu_mean);
    if (nu <= 0) {
        /* Poisson fluctuation: no participants this event.  Primary still
         * absorbed; caller sees ABRASION with n_secondaries = 0. */
        return;
    }
    if (nu > OSH_NUCLEAR_MAX_SECONDARIES) {
        nu = OSH_NUCLEAR_MAX_SECONDARIES;
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

    /* Minimal fast-emission model: each participant is sampled against the
     * original beam axis, while the available proton energy is degraded after
     * each emission so total emitted kinetic energy remains bounded. */
    T_current = T_lab_mev;
    n_knockout_p = 0;
    n_knockout_n = 0;

    for (j = 0; j < nu; ++j) {
        /* Stop if the remaining proton energy is too low for meaningful
         * kinematics (below ~1 MeV the equal-mass formula loses accuracy
         * and the nucleon would range out immediately anyway). */
        if (T_current < 1.0) {
            break;
        }

        /* Nucleon type: neutron with prob N/A, proton with prob Z/A. */
        is_neutron = (osh_rng_double(rng) < n_over_a);
        species = is_neutron ? &s_neutron : &s_proton;
        if (is_neutron) {
            ++n_knockout_n;
        } else {
            ++n_knockout_p;
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

        /* Rotate the recoil nucleon into the lab frame.
         * Convention matches pp elastic: secondary uses the opposite azimuth
         * (-cos_phi, -sin_phi) so that momentum is conserved transversely. */
        sin_sec = sqrt(fmax(0.0, 1.0 - cos_sec * cos_sec));
        osh_kinematics_rotate_dir_cos(
            incident_dir, event_out->secondaries[j].dir, cos_sec, sin_sec, -cos_phi, -sin_phi);

        event_out->secondaries[j].energy = e_sec;
        event_out->secondaries[j].species = species;
        event_out->n_secondaries = (size_t) j + 1u;

        /* Degrade the proton for the next collision. */
        T_current = e_prim;
    }

    if (event_out->n_fragments > 0u) {
        unsigned int da = (unsigned int) (n_knockout_p + n_knockout_n);
        unsigned int dz = (unsigned int) n_knockout_p;

        event_out->fragments[0].a = (event_out->fragments[0].a > da) ? event_out->fragments[0].a - da : 0u;
        event_out->fragments[0].z = (event_out->fragments[0].z > dz) ? event_out->fragments[0].z - dz : 0u;
        if (event_out->fragments[0].a == 0u || event_out->fragments[0].z == 0u) {
            event_out->n_fragments = 0u;
        }
    }
}
