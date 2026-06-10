#include "physics/atomic/osh_physics_moliere.h"

#include <math.h>

#include "common/osh_vect.h"
#include "random/osh_rng.h"

/* ---- Highland formula ----------------------------------------------------- */

double
osh_physics_moliere_theta0(double t_total_mev, double mass_mev, double z_eff, double thickness_gcm2, double x0_gcm2) {
    double momentum; /* p = sqrt(T*(T + 2m))                      [MeV/c] */
    double e_total;  /* E = T + m                                  [MeV]   */
    double pv;       /* βcp = p²/E = T(T+2m)/(T+m)               [MeV]   */
    double d_over_x0;
    double log_arg;
    double theta0;

    if (x0_gcm2 <= 0.0 || thickness_gcm2 <= 0.0 || z_eff <= 0.0) {
        return 0.0;
    }

    /*
     * Relativistic kinematics.
     *
     * p = sqrt(T*(T + 2m)),  E = T + m,  β = p/E.
     * βcp = β × c × p = p²/E = T*(T+2m)/(T+m)  [MeV] in natural units (c=1).
     *
     * This is the denominator of the Highland formula as given in PDG eq. 34.15.
     */
    momentum = sqrt(t_total_mev * (t_total_mev + 2.0 * mass_mev));
    e_total = t_total_mev + mass_mev;
    pv = momentum * momentum / e_total; /* = βcp */

    if (pv <= 0.0) {
        return 0.0;
    }

    d_over_x0 = thickness_gcm2 / x0_gcm2;

    /*
     * Highland formula [Hig75, PDG34.3]:
     *
     *   θ₀ = (13.6 MeV / βcp) × z_eff × √(d/X₀) × [1 + 0.038 ln(z_eff² d/X₀)]
     *
     * The log argument is z_eff² × d/X₀; for fully stripped ions z_eff = Z_proj
     * at high energies, and the log correction is typically small (< 10%).
     *
     * If the correction factor becomes negative (very thin slab), clamp to 0.
     */
    log_arg = z_eff * z_eff * d_over_x0;
    theta0 = (13.6 / pv) * z_eff * sqrt(d_over_x0) * (1.0 + 0.038 * log(log_arg));

    /*
     * For ultra-thin slabs (log_arg << 1) the log correction term is a large
     * negative number and the correction factor (1 + 0.038*ln(...)) can become
     * negative before the Highland approximation has any physical meaning.
     * Clamp to zero in that regime rather than returning a negative angle.
     */
    return (theta0 > 0.0) ? theta0 : 0.0;
}

/* ---- MCS substep criterion ----------------------------------------------- */

double osh_physics_moliere_s_theta(
    double t_total_mev, double mass_mev, double z_eff, double rho_gcm3, double x0_gcm2, double theta_max_rad) {
    double pv; /* βcp = p²/E = T(T+2m)/(T+m)  [MeV] */
    double momentum;
    double e_total;
    double ratio;

    if (rho_gcm3 <= 0.0 || z_eff <= 0.0 || x0_gcm2 <= 0.0 || theta_max_rad <= 0.0) {
        return 0.0;
    }

    /*
     * βcp = p²/E, identical to the kinematics in theta0().
     * For a stopped particle (T → 0) pv → 0 and s_theta → 0 as well —
     * the criterion becomes infinitely tight, but the particle is dead.
     */
    momentum = sqrt(t_total_mev * (t_total_mev + 2.0 * mass_mev));
    e_total = t_total_mev + mass_mev;
    pv = momentum * momentum / e_total;

    if (pv <= 0.0) {
        return 0.0;
    }

    /*
     * Invert θ₀ ≈ (13.6/βcp) × z_eff × √(d/X₀) for d, ignoring the log
     * correction.  Then convert areal density to geometric length:
     *
     *   s_theta = (x0_gcm2 / rho_gcm3) × (θ_max × βcp / (13.6 × z_eff))²
     */
    ratio = theta_max_rad * pv / (13.6 * z_eff);
    return (x0_gcm2 / rho_gcm3) * ratio * ratio;
}

/* ---- Direction scatter ---------------------------------------------------- */

void osh_physics_moliere_scatter(double const v[3], double w[3], double theta0, struct osh_rng *rng) {
    double u1[3]; /* first transverse basis vector  */
    double u2[3]; /* second transverse basis vector  */
    double tx;    /* projected angle in the u1 plane [rad] */
    double ty;    /* projected angle in the u2 plane [rad] */
    double norm;

    if (theta0 < 1e-9) {
        w[0] = v[0];
        w[1] = v[1];
        w[2] = v[2];
        return;
    }

    osh_vect_orthogonal_basis_norm(v, u1, u2);

    /*
     * Sample two independent projected scattering angles from N(0, θ₀).
     *
     * The Highland distribution is Gaussian in each of the two transverse
     * projected planes independently [Hig75].  osh_rng_gauss01() uses
     * Box-Muller with a cached spare, so the pair (tx, ty) costs one
     * log+sqrt and one stored value — essentially one call's worth of work.
     */
    tx = osh_rng_gauss01(rng) * theta0;
    ty = osh_rng_gauss01(rng) * theta0;

    /*
     * Scattered direction in the caller's coordinate frame:
     *   w_raw = v + tx·u1 + ty·u2
     *
     * For |θ| << 1 this is the small-angle approximation v + δv, accurate to
     * O(θ²).  Normalisation recovers a unit vector even for large tail draws,
     * keeping the direction physically valid regardless of the Gaussian sample.
     */
    w[0] = v[0] + tx * u1[0] + ty * u2[0];
    w[1] = v[1] + tx * u1[1] + ty * u2[1];
    w[2] = v[2] + tx * u1[2] + ty * u2[2];

    norm = sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
    if (norm > 0.0) {
        w[0] /= norm;
        w[1] /= norm;
        w[2] /= norm;
    } else {
        /* Degenerate (tx = ty = 0 and v = 0 — cannot happen in practice) */
        w[0] = v[0];
        w[1] = v[1];
        w[2] = v[2];
    }
}
