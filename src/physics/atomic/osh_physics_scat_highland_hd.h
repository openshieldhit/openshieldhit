#ifndef OSH_PHYSICS_SCAT_HIGHLAND_HD_H
#define OSH_PHYSICS_SCAT_HIGHLAND_HD_H

/*
 * osh_physics_scat_highland_hd.h — device-compilable Highland MCS.
 *
 * Bodies are marked OSH_HD static inline so they compile both as host
 * functions (via plain C compilation) and as device functions (via nvcc
 * with __host__ __device__).  The original .c file includes this header
 * and re-exports each with its unchanged public signature.
 */

#include "common/osh_hd.h"
#include "physics/atomic/osh_physics_scat_highland.h"

#include <math.h>

#include "common/osh_vect.h"
#include "common/osh_vect_hd.h"
#include "random/osh_rng.h"
#include "random/osh_rng_hd.h"

OSH_HD static inline double _osh_physics_highland_theta0_hd(
    double t_total_mev, double mass_mev, double z_eff, double thickness_gcm2, double path_scale_gcm2, double x0_gcm2) {
    double momentum;
    double e_total;
    double pv;
    double d_over_x0;
    double s_over_x0;
    double log_arg;
    double theta0;

    if (x0_gcm2 <= 0.0 || thickness_gcm2 <= 0.0 || z_eff <= 0.0) {
        return 0.0;
    }
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

OSH_HD static inline double _osh_physics_highland_s_theta_hd(
    double t_total_mev, double mass_mev, double z_eff, double rho_gcm3, double x0_gcm2, double theta_max_rad) {
    double pv;
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

OSH_HD static inline void _osh_physics_highland_scatter_hd(
    double const v[3], double w[3], double theta0, struct osh_rng *rng) {
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

    _osh_vect_orthogonal_basis_norm_hd(v, u1, u2);

    tx = _osh_rng_gauss01_hd(rng) * theta0;
    ty = _osh_rng_gauss01_hd(rng) * theta0;

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

#endif /* OSH_PHYSICS_SCAT_HIGHLAND_HD_H */
