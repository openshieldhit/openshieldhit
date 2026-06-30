#include "physics/atomic/osh_physics_scat.h"

#include "physics/atomic/osh_physics_scat_highland.h"
#include "physics/atomic/osh_physics_scat_moliere.h"

int osh_physics_mcs_scatter(enum osh_mcs_model model,
                            double const v[3],
                            double w[3],
                            double t_kin_mev,
                            double mass_mev,
                            double z_eff,
                            double d_gcm2,
                            double path_scale_gcm2,
                            double x0_gcm2,
                            double chic2_coeff,
                            double screen_z,
                            struct osh_rng *rng) {
    double theta0;

    w[0] = v[0];
    w[1] = v[1];
    w[2] = v[2];

    /* Highland θ₀ is the validated scattering width; both Gaussian and Molière
     * use it for magnitude (Molière adds only the distribution shape/tail). */
    theta0 = osh_physics_highland_theta0(t_kin_mev, mass_mev, z_eff, d_gcm2, path_scale_gcm2, x0_gcm2);
    if (theta0 <= 0.0) {
        return 0;
    }

    if (model == OSH_MCS_GAUSSIAN) {
        osh_physics_highland_scatter(v, w, theta0, rng);
        return 1;
    }
    if (model == OSH_MCS_MOLIERE) {
        return osh_physics_moliere_scatter(
            v, w, theta0, t_kin_mev, mass_mev, z_eff, d_gcm2, path_scale_gcm2, chic2_coeff, screen_z, rng);
    }
    return 0;
}
