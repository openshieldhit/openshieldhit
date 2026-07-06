#include "physics/atomic/osh_physics_scat_highland.h"

#include <math.h>

#include "common/osh_vect.h"
#include "random/osh_rng.h"
#include "physics/atomic/osh_physics_scat_highland_hd.h"

double osh_physics_highland_theta0(
    double t_total_mev, double mass_mev, double z_eff, double thickness_gcm2, double path_scale_gcm2, double x0_gcm2) {
    return _osh_physics_highland_theta0_hd(
        t_total_mev, mass_mev, z_eff, thickness_gcm2, path_scale_gcm2, x0_gcm2);
}

double osh_physics_highland_s_theta(
    double t_total_mev, double mass_mev, double z_eff, double rho_gcm3, double x0_gcm2, double theta_max_rad) {
    return _osh_physics_highland_s_theta_hd(
        t_total_mev, mass_mev, z_eff, rho_gcm3, x0_gcm2, theta_max_rad);
}

void osh_physics_highland_scatter(double const v[3], double w[3], double theta0, struct osh_rng *rng) {
    _osh_physics_highland_scatter_hd(v, w, theta0, rng);
}
