#include "physics/nuclear/osh_nuclear_abrasion.h"
#include "physics/nuclear/osh_nuclear_handler.h"

#include <math.h> /* sqrt */

#include "particle/osh_particle.h"
#include "particle/osh_particle_const.h"
#include "particle/osh_particle_pdg.h"
#include "physics/osh_kinematics.h"
#include "random/osh_rng.h"

#define OSH_MB_TO_CM2 1.0e-27

/* Species descriptors initialised on first call; read-only thereafter. */
static struct particle s_proton;
static struct particle s_neutron;
static int s_species_valid = 0;

static void ensure_species(void) {
    if (!s_species_valid) {
        osh_particle_from_pdg(&s_proton, OSH_PART_PDG_PROTON);
        osh_particle_from_pdg(&s_neutron, OSH_PART_PDG_NEUTRON);
        s_species_valid = 1;
    }
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
    double cos_prim; /* lab cos of the continuing (degraded) proton after collision j */
    double e_prim;   /* lab KE of the continuing proton [MeV] */
    double cos_sec;  /* lab cos of the knocked-out nucleon */
    double e_sec;    /* lab KE of the knocked-out nucleon [MeV] */
    double sin_sec;
    int is_neutron;
    int nu;          /* sampled number of participant nucleons */
    int j;
    struct particle const *species;

    ensure_species();

    /* Primary is always absorbed in an inelastic event.  Secondaries are filled
     * in below; n_secondaries stays 0 if nu=0 (Poisson fluctuation). */
    event_out->kind = OSH_NUCLEAR_EVENT_ABRASION;
    event_out->primary_energy = 0.0;
    event_out->n_secondaries = 0u;

    if (sigma_pa_cm2 <= 0.0 || a_eff <= 0.0) {
        return;
    }

    /* Wounded-nucleon mean (Glauber 1970):
     *   <nu> = sigma_pN * A / sigma_pA
     * sigma_pN is a tunable constant (OSH_ABRASION_SIGMA_PN_MB); calibrate
     * against FLUKA / SH12A secondary proton yield data. */
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

    /* Sequential collision model: each participant scatters off the
     * incoming proton in turn.  The proton energy degrades after each
     * collision so that energy is conserved across all nu steps.
     * (This is equivalent to a one-dimensional intranuclear cascade with
     * nu steps and isotropic CM scattering.) */
    T_current = T_lab_mev;

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

        /* Isotropic CM scattering — simplest physical model.
         * A forward-peaked angular distribution (e.g. from a diffraction
         * model) would be a refinement for future work. */
        cos_cm = 2.0 * osh_rng_double(rng) - 1.0;
        osh_kinematics_azimuth(rng, &cos_phi, &sin_phi);

        /* Full relativistic equal-mass kinematics (proton mass used for both
         * p and n since mp ≈ mn to better than 0.1%). */
        osh_kinematics_elastic_equal_mass_lab(T_current, OSH_PART_MASS_PROTON,
                                              cos_cm,
                                              &cos_prim, &e_prim,
                                              &cos_sec, &e_sec);

        /* Rotate the recoil nucleon into the lab frame.
         * Convention matches pp elastic: secondary uses the opposite azimuth
         * (-cos_phi, -sin_phi) so that momentum is conserved transversely. */
        sin_sec = sqrt(1.0 - cos_sec * cos_sec);
        osh_kinematics_rotate_dir_cos(incident_dir,
                                      event_out->secondaries[j].dir,
                                      cos_sec, sin_sec,
                                      -cos_phi, -sin_phi);

        event_out->secondaries[j].energy = e_sec;
        event_out->secondaries[j].species = species;
        event_out->n_secondaries = (size_t)(j + 1);

        /* Degrade the proton for the next collision. */
        T_current = e_prim;
    }
}
