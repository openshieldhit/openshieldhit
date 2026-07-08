#include "physics/atomic/osh_physics_strag.h"

#include <math.h>

#include "openshieldhit/const.h"

double osh_physics_strag_xi(double z_eff, double z_over_a, double thickness_gcm2, double beta2) {
    if (z_eff <= 0.0 || z_over_a <= 0.0 || thickness_gcm2 <= 0.0 || beta2 <= 0.0) {
        return 0.0;
    }
    /* ξ = (K/2)·(Z/A)·(z_eff²/β²)·d   [MeV] */
    return 0.5 * OSH_K_BETHE * z_over_a * (z_eff * z_eff / beta2) * thickness_gcm2;
}

double osh_physics_strag_emax(double t_kin_mev, double mass_mev) {
    double gamma; /* Lorentz factor γ = 1 + T/M                    */
    double beta2; /* projectile β² = 1 − 1/γ²                      */

    if (t_kin_mev <= 0.0 || mass_mev <= 0.0) {
        return 0.0;
    }
    gamma = 1.0 + t_kin_mev / mass_mev;
    beta2 = 1.0 - 1.0 / (gamma * gamma);
    return osh_physics_strag_emax_beta2(t_kin_mev, mass_mev, beta2);
}

double osh_physics_strag_emax_beta2(double t_kin_mev, double mass_mev, double beta2) {
    double gamma; /* Lorentz factor γ = 1 + T/M                    */
    double ratio; /* electron-to-projectile mass ratio mₑ/M       */

    if (t_kin_mev <= 0.0 || mass_mev <= 0.0 || beta2 <= 0.0 || beta2 >= 1.0) {
        return 0.0;
    }
    gamma = 1.0 + t_kin_mev / mass_mev;
    ratio = OSH_ELECTRON_MASS_MEV / mass_mev;
    /* E_max = 2·mₑc²·β²γ² / (1 + 2γ·mₑ/M + (mₑ/M)²) */
    return 2.0 * OSH_ELECTRON_MASS_MEV * beta2 * gamma * gamma / (1.0 + 2.0 * gamma * ratio + ratio * ratio);
}

double osh_physics_strag_kappa(double xi, double e_max) {
    if (e_max <= 0.0 || xi <= 0.0) {
        return 0.0;
    }
    return xi / e_max;
}

double osh_physics_strag_lambda_bar(double kappa, double beta2) {
    if (kappa <= 0.0) {
        /* NAN, not 0.0: unlike ξ/E_max/κ (positive magnitudes), 0.0 is a
         * legitimate λ̄, so a distinct sentinel is needed for invalid input. */
        return NAN;
    }
    return osh_physics_strag_lambda_bar_log(log(kappa), beta2);
}

double osh_physics_strag_lambda_bar_log(double log_kappa, double beta2) {
    if (!isfinite(log_kappa)) {
        return NAN;
    }
    /* λ̄ = −β² − ln κ + (γ_Euler − 1) */
    return -beta2 - log_kappa + (OSH_EULER_GAMMA - 1.0);
}
