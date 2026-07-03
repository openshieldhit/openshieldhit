#include "physics/osh_kinematics.h"

#include <math.h>
#include <stddef.h>

#include "common/osh_vect.h"
#include "random/osh_rng.h"

/* ---- Direction rotation -------------------------------------------------- */

void osh_kinematics_rotate_dir_cos(
    double const v[3], double w[3], double cos_theta, double sin_theta, double a, double b) {
    double u1[3];
    double u2[3];
    double norm;
    double tmp[3];

    if (sin_theta < 1e-9) {
        w[0] = v[0];
        w[1] = v[1];
        w[2] = v[2];
        return;
    }

    osh_vect_orthogonal_basis_norm(v, u1, u2);

    /* w_raw = cos_theta·v + sin_theta·(a·u1 + b·u2) */
    tmp[0] = cos_theta * v[0] + sin_theta * (a * u1[0] + b * u2[0]);
    tmp[1] = cos_theta * v[1] + sin_theta * (a * u1[1] + b * u2[1]);
    tmp[2] = cos_theta * v[2] + sin_theta * (a * u1[2] + b * u2[2]);

    /* Renormalise for numerical robustness; result may alias v */
    norm = sqrt(tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2]);
    if (norm > 0.0) {
        w[0] = tmp[0] / norm;
        w[1] = tmp[1] / norm;
        w[2] = tmp[2] / norm;
    } else {
        w[0] = v[0];
        w[1] = v[1];
        w[2] = v[2];
    }
}

/* ---- Azimuthal direction ------------------------------------------------- */

void osh_kinematics_azimuth(struct osh_rng *rng, double *cos_phi, double *sin_phi) {
    double Q, R, QQ, RR, QR;

    /*
     * Knuth disk-rejection (TAOCP Vol. 2, §3.4.1):
     *   Q = 2·U₁ − 1  ∈ [-1, 1)
     *   R = U₂         ∈ [0,  1)   (upper half of unit square)
     * Accept if Q²+R² ≤ 1.  Average ~1.27 uniform deviates consumed.
     *
     * Weierstrass substitution:
     *   cos φ = (Q²−R²)/(Q²+R²)
     *   sin φ = 2QR/(Q²+R²)
     * No sqrt, no trig.
     */
    do {
        Q = 2.0 * osh_rng_double(rng) - 1.0;
        R = osh_rng_double(rng);
        QQ = Q * Q;
        RR = R * R;
        QR = QQ + RR;
    } while (QR > 1.0 || QR == 0.0);

    *cos_phi = (QQ - RR) / QR;
    *sin_phi = (2.0 * Q * R) / QR;
}

/* ---- Relativistic 2-body elastic kinematics ----------------------------- */

void osh_kinematics_elastic_equal_mass_lab(double T_lab,
                                           double mp,
                                           double cos_theta_cm,
                                           double *cos_theta1_lab,
                                           double *e1_lab_mev,
                                           double *cos_theta2_lab,
                                           double *e2_lab_mev) {
    double p_lab;        /* lab momentum of incident particle [MeV/c] */
    double W;            /* invariant mass sqrt(s) [MeV] */
    double p_cm;         /* CM momentum magnitude [MeV/c] */
    double E_cm;         /* CM energy of each particle = W/2 [MeV] */
    double beta_cm;      /* CM velocity β = p_lab / (T + 2mp) */
    double gamma_cm;     /* CM Lorentz factor γ = (T + 2mp) / W */
    double sin_theta_cm; /* sin(θ_CM) derived from sqrt(1 - cos²) */
    double p1z_cm;       /* z-component of particle 1 momentum in CM */
    double p1t_cm;       /* transverse momentum of particle 1 in CM */
    double p1z_lab;      /* z-component of particle 1 momentum in lab */
    double p2z_lab;      /* z-component of particle 2 momentum in lab */
    double E1_lab;       /* total energy of particle 1 in lab */
    double E2_lab;       /* total energy of particle 2 in lab */
    double p1_tot_lab;   /* total momentum magnitude of particle 1 in lab */
    double p2_tot_lab;   /* total momentum magnitude of particle 2 in lab */

    /*
     * Lab kinematics: incident particle at energy T_lab, target at rest.
     *   p_lab = sqrt(T(T + 2m))
     *   s     = 2m(T + 2m)     → W = sqrt(2m(T + 2m))
     *
     * CM frame:
     *   p_cm   = m·p_lab/W
     *   E_cm   = W/2  (symmetric for equal masses)
     *   β_cm   = p_lab / (T + 2m)
     *   γ_cm   = (T + 2m) / W
     */
    p_lab = sqrt(T_lab * (T_lab + 2.0 * mp));
    W = sqrt(2.0 * mp * (T_lab + 2.0 * mp));
    p_cm = (W > 0.0) ? (p_lab * mp / W) : 0.0;
    E_cm = W * 0.5;
    beta_cm = p_lab / (T_lab + 2.0 * mp);
    gamma_cm = (T_lab + 2.0 * mp) / W;

    sin_theta_cm = sqrt(1.0 - cos_theta_cm * cos_theta_cm);

    p1z_cm = p_cm * cos_theta_cm;
    p1t_cm = p_cm * sin_theta_cm;

    /*
     * Lorentz boost along z from CM to lab.
     *
     * Particle 1 (scattered incident):
     *   E1_lab  = γ·(E_cm + β·p1z_cm)
     *   p1z_lab = γ·(p1z_cm + β·E_cm)
     *   p1t_lab = p1t_cm             (transverse unchanged)
     *
     * Particle 2 (recoil, back-to-back in CM: p2z_cm = -p1z_cm):
     *   E2_lab  = γ·(E_cm − β·p1z_cm)
     *   p2z_lab = γ·(−p1z_cm + β·E_cm)
     *   p2t_lab = p1t_cm             (same transverse magnitude)
     */
    E1_lab = gamma_cm * (E_cm + beta_cm * p1z_cm);
    p1z_lab = gamma_cm * (p1z_cm + beta_cm * E_cm);

    E2_lab = gamma_cm * (E_cm - beta_cm * p1z_cm);
    /*
     * p2z_lab = γ(-p1z_cm + β·E_cm).  In exact arithmetic β·E_cm == p_cm,
     * but computing β·E_cm separately introduces floating-point error that can
     * make p2z_lab slightly negative near cos_theta_cm=1 (catastrophic
     * cancellation).  Using p_cm directly avoids this: the result is then
     * γ·p_cm·(1 − cos_theta_cm) which is manifestly non-negative.
     */
    p2z_lab = gamma_cm * (p_cm - p1z_cm);

    /* p1t_lab = p2t_lab = p1t_cm (transverse magnitude unchanged by boost) */

    p1_tot_lab = sqrt(p1z_lab * p1z_lab + p1t_cm * p1t_cm);
    p2_tot_lab = sqrt(p2z_lab * p2z_lab + p1t_cm * p1t_cm);

    *e1_lab_mev = E1_lab - mp;
    *e2_lab_mev = E2_lab - mp;

    *cos_theta1_lab = (p1_tot_lab > 0.0) ? (p1z_lab / p1_tot_lab) : 1.0;
    *cos_theta2_lab = (p2_tot_lab > 0.0) ? (p2z_lab / p2_tot_lab) : 0.0;

    /* Guard against tiny negative kinetic energies from floating-point rounding */
    if (*e1_lab_mev < 0.0)
        *e1_lab_mev = 0.0;
    if (*e2_lab_mev < 0.0)
        *e2_lab_mev = 0.0;
}

void osh_kinematics_elastic_lab(double T_lab,
                                double m1,
                                double m2,
                                double cos_theta_cm,
                                double *cos_theta1_lab,
                                double *e1_lab_mev,
                                double *cos_theta2_lab,
                                double *e2_lab_mev) {
    double p_lab;        /* incident projectile lab momentum [MeV/c]           */
    double E1_lab_in;    /* incident projectile total lab energy [MeV]         */
    double s;            /* Mandelstam s = W²                                  */
    double W;            /* invariant mass sqrt(s) [MeV]                       */
    double p_cm;         /* CM momentum of either particle (elastic) [MeV/c]   */
    double E1_cm;        /* projectile CM total energy [MeV]                   */
    double E2_cm;        /* target CM total energy [MeV]                       */
    double beta_cm;      /* CM velocity                                        */
    double gamma_cm;     /* CM Lorentz factor                                  */
    double sin_theta_cm; /* sin θ_CM                                           */
    double p1z_cm;       /* projectile CM longitudinal momentum               */
    double p1t;          /* transverse momentum (shared magnitude)            */
    double E1_lab;
    double E2_lab;
    double p1z_lab;
    double p2z_lab;
    double p1_tot;
    double p2_tot;

    p_lab = sqrt(T_lab * (T_lab + 2.0 * m1));
    E1_lab_in = T_lab + m1;
    s = (m1 * m1) + (m2 * m2) + (2.0 * m2 * E1_lab_in);
    W = sqrt(s);
    p_cm = (W > 0.0) ? (p_lab * m2 / W) : 0.0;
    E1_cm = sqrt((p_cm * p_cm) + (m1 * m1));
    E2_cm = sqrt((p_cm * p_cm) + (m2 * m2));
    beta_cm = p_lab / (E1_lab_in + m2);
    gamma_cm = (E1_lab_in + m2) / W;

    sin_theta_cm = sqrt(fmax(0.0, 1.0 - (cos_theta_cm * cos_theta_cm)));
    p1z_cm = p_cm * cos_theta_cm;
    p1t = p_cm * sin_theta_cm;

    /* Boost the scattered projectile (particle 1) and the recoil target
     * (particle 2, back-to-back in CM: p2z_cm = -p1z_cm) from CM to lab. */
    E1_lab = gamma_cm * (E1_cm + (beta_cm * p1z_cm));
    p1z_lab = gamma_cm * (p1z_cm + (beta_cm * E1_cm));
    E2_lab = gamma_cm * (E2_cm - (beta_cm * p1z_cm));
    p2z_lab = gamma_cm * ((beta_cm * E2_cm) - p1z_cm);

    p1_tot = sqrt((p1z_lab * p1z_lab) + (p1t * p1t));
    p2_tot = sqrt((p2z_lab * p2z_lab) + (p1t * p1t));

    *e1_lab_mev = E1_lab - m1;
    *e2_lab_mev = E2_lab - m2;
    *cos_theta1_lab = (p1_tot > 0.0) ? (p1z_lab / p1_tot) : 1.0;
    *cos_theta2_lab = (p2_tot > 0.0) ? (p2z_lab / p2_tot) : 1.0;

    if (*e1_lab_mev < 0.0)
        *e1_lab_mev = 0.0;
    if (*e2_lab_mev < 0.0)
        *e2_lab_mev = 0.0;
}

/* ---- Two-body decay ------------------------------------------------------- */

double osh_kinematics_two_body_decay_p(double m_parent, double m1, double m2) {
    double sum;
    double diff;
    double num;

    if (m_parent <= 0.0) {
        return 0.0;
    }

    sum = m1 + m2;
    diff = m1 - m2;
    num = (m_parent * m_parent - sum * sum) * (m_parent * m_parent - diff * diff);
    if (num <= 0.0) {
        return 0.0;
    }
    return sqrt(num) / (2.0 * m_parent);
}

/* ---- Lorentz boost to lab ------------------------------------------------- */

void osh_kinematics_boost_to_lab(double m_parent,
                                 double const p_parent[3],
                                 double e_cm,
                                 double const p_cm[3],
                                 double *e_lab_out,
                                 double p_lab_out[3]) {
    double p2;
    double e_parent;
    double dot;
    double coef;

    p2 = p_parent[0] * p_parent[0] + p_parent[1] * p_parent[1] + p_parent[2] * p_parent[2];
    e_parent = sqrt(m_parent * m_parent + p2);
    dot = p_cm[0] * p_parent[0] + p_cm[1] * p_parent[1] + p_cm[2] * p_parent[2];

    /* coef = (p_cm·P)/(M(E_P+M)) + e_cm/M; for P = 0 the whole P-term vanishes */
    coef = dot / (m_parent * (e_parent + m_parent)) + e_cm / m_parent;

    *e_lab_out = (e_parent * e_cm + dot) / m_parent;
    p_lab_out[0] = p_cm[0] + p_parent[0] * coef;
    p_lab_out[1] = p_cm[1] + p_parent[1] * coef;
    p_lab_out[2] = p_cm[2] + p_parent[2] * coef;
}
