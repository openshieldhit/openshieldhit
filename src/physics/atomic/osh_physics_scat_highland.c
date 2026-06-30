#include "physics/atomic/osh_physics_scat_highland.h"

#include <math.h>

#include "common/osh_vect.h"
#include "random/osh_rng.h"

double osh_physics_highland_theta0(
    double t_total_mev, double mass_mev, double z_eff, double thickness_gcm2, double path_scale_gcm2, double x0_gcm2) {
    double momentum; /* p = sqrt(T*(T + 2m)) [MeV/c] */
    double e_total;  /* E = T + m [MeV] */
    double pv;       /* beta*c*p = p^2/E [MeV] */
    double d_over_x0;
    double s_over_x0;
    double log_arg;
    double theta0;

    if (x0_gcm2 <= 0.0 || thickness_gcm2 <= 0.0 || z_eff <= 0.0) {
        return 0.0;
    }
    /* The log correction is evaluated at the macroscopic path scale, not the
     * substep thickness, so that summing the per-substep variances reproduces
     * the full-path Highland value.  Falling back to the substep thickness
     * keeps the legacy single-step behaviour. */
    if (path_scale_gcm2 <= 0.0) {
        path_scale_gcm2 = thickness_gcm2;
    }

    momentum = sqrt(t_total_mev * (t_total_mev + 2.0 * mass_mev));
    e_total = t_total_mev + mass_mev;
    pv = momentum * momentum / e_total;

    if (pv <= 0.0) {
        return 0.0;
    }

    d_over_x0 = thickness_gcm2 / x0_gcm2;
    s_over_x0 = path_scale_gcm2 / x0_gcm2;
    log_arg = z_eff * z_eff * s_over_x0;
    theta0 = (13.6 / pv) * z_eff * sqrt(d_over_x0) * (1.0 + 0.038 * log(log_arg));

    return (theta0 > 0.0) ? theta0 : 0.0;
}

double osh_physics_highland_s_theta(
    double t_total_mev, double mass_mev, double z_eff, double rho_gcm3, double x0_gcm2, double theta_max_rad) {
    double pv; /* beta*c*p = p^2/E [MeV] */
    double momentum;
    double e_total;
    double ratio;

    if (rho_gcm3 <= 0.0 || z_eff <= 0.0 || x0_gcm2 <= 0.0 || theta_max_rad <= 0.0) {
        return 0.0;
    }

    momentum = sqrt(t_total_mev * (t_total_mev + 2.0 * mass_mev));
    e_total = t_total_mev + mass_mev;
    pv = momentum * momentum / e_total;

    if (pv <= 0.0) {
        return 0.0;
    }

    ratio = theta_max_rad * pv / (13.6 * z_eff);
    return (x0_gcm2 / rho_gcm3) * ratio * ratio;
}

void osh_physics_highland_scatter(double const v[3], double w[3], double theta0, struct osh_rng *rng) {
    double u1[3];
    double u2[3];
    double tx;
    double ty;
    double norm;

    if (theta0 < 1e-9) {
        w[0] = v[0];
        w[1] = v[1];
        w[2] = v[2];
        return;
    }

    osh_vect_orthogonal_basis_norm(v, u1, u2);

    tx = osh_rng_gauss01(rng) * theta0;
    ty = osh_rng_gauss01(rng) * theta0;

    w[0] = v[0] + tx * u1[0] + ty * u2[0];
    w[1] = v[1] + tx * u1[1] + ty * u2[1];
    w[2] = v[2] + tx * u1[2] + ty * u2[2];

    norm = sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
    if (norm > 0.0) {
        w[0] /= norm;
        w[1] /= norm;
        w[2] /= norm;
    } else {
        w[0] = v[0];
        w[1] = v[1];
        w[2] = v[2];
    }
}
