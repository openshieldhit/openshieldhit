#include "physics/atomic/osh_physics_scat_moliere.h"

#include <math.h>

#include "openshieldhit/const.h"
#include "physics/osh_kinematics.h"
#include "random/osh_rng.h"

#define OSH_MOLIERE_NF 257    /* reduced-angle (ϑ) grid points for f⁽ⁿ⁾ */
#define OSH_MOLIERE_THMAX 12.0 /* maximum reduced angle ϑ */

/* Inverse-CDF support table: ϑ as a function of (B, u).  The reduced-angle
 * density depends on B only through f(ϑ) = f⁽⁰⁾ + f⁽¹⁾/B + f⁽²⁾/B², so ϑ(u) is a
 * smooth function of (B, u) that we tabulate once at init and bilinearly
 * interpolate on the hot path — no per-scatter CDF reconstruction.  Built from
 * our own f⁽ⁿ⁾ formulas (no external tables). */
#define OSH_MOLIERE_NB 64    /* B-grid nodes */
#define OSH_MOLIERE_NU 1024  /* u-grid nodes (cumulative probability) */
#define OSH_MOLIERE_B_MIN 1.6 /* ≈ solve_b threshold (ln Ω = 1.0816 ⇒ B ≈ 1.6) */
#define OSH_MOLIERE_B_MAX 40.0 /* above this the 1/B, 1/B² corrections are negligible */

static double osh_moliere_f0[OSH_MOLIERE_NF];
static double osh_moliere_f1[OSH_MOLIERE_NF];
static double osh_moliere_f2[OSH_MOLIERE_NF];
static double osh_moliere_dth;
static float osh_moliere_invcdf[OSH_MOLIERE_NB][OSH_MOLIERE_NU];
static double osh_moliere_inv_db; /* (NB-1) / (B_MAX - B_MIN) */
static int osh_moliere_inited;

/* Bessel J0(x), Abramowitz & Stegun 9.4.1 / 9.4.3 polynomial approximation. */
static double osh_bessel_j0(double x) {
    double ax;
    double y;
    double z;
    double f;
    double theta;

    ax = fabs(x);
    if (ax < 3.0) {
        y = (x / 3.0) * (x / 3.0);
        return 1.0
               + y
                     * (-2.2499997
                        + y * (1.2656208 + y * (-0.3163866 + y * (0.0444479 + y * (-0.0039444 + y * 0.0002100)))));
    }

    z = 3.0 / ax;
    f = 0.79788456
        + z
              * (-0.00000077
                 + z * (-0.00552740 + z * (-0.00009512 + z * (0.00137237 + z * (-0.00072805 + z * 0.00014476)))));
    theta = ax - 0.78539816
            + z
                  * (-0.04166397
                     + z * (-0.00003954 + z * (0.00262573 + z * (-0.00054125 + z * (-0.00029333 + z * 0.00013558)))));
    return f * cos(theta) / sqrt(ax);
}

/**
 * @brief Build the inverse CDF ϑ(u) for one B value onto a uniform u grid.
 *
 * @details
 * Forms the (negative-clamped) cumulative of the sampling density ϑ·f(ϑ;B) over
 * the reduced-angle grid, then inverts it onto NU uniform cumulative-probability
 * nodes.  Called only at init, once per B node; not on the hot path.
 *
 * @param[in]  b    Molière B parameter for this row.
 * @param[out] row  Destination of NU interpolated ϑ values.
 */
static void moliere_build_invcdf_row(double b, float *row) {
    double cdf[OSH_MOLIERE_NF];
    double inv_b;
    double inv_b2;
    double acc;
    double total;
    double th;
    double pdf;
    double target;
    double seg;
    int i;
    int k;
    int idx;

    inv_b = 1.0 / b;
    inv_b2 = inv_b * inv_b;
    acc = 0.0;
    for (i = 0; i < OSH_MOLIERE_NF; ++i) {
        th = (double) i * osh_moliere_dth;
        pdf = th * (osh_moliere_f0[i] + osh_moliere_f1[i] * inv_b + osh_moliere_f2[i] * inv_b2);
        if (pdf < 0.0) {
            pdf = 0.0; /* clamp rare negative dips so the CDF stays monotone */
        }
        acc += pdf * osh_moliere_dth;
        cdf[i] = acc;
    }

    total = cdf[OSH_MOLIERE_NF - 1];
    if (total <= 0.0) {
        for (k = 0; k < OSH_MOLIERE_NU; ++k) {
            row[k] = 0.0f;
        }
        return;
    }

    /* Invert: for each uniform u_k = k/(NU-1) find ϑ with cdf(ϑ) = u_k·total.
     * cdf is monotone non-decreasing, so advance idx forward as u_k rises. */
    row[0] = 0.0f;
    idx = 1;
    for (k = 1; k < OSH_MOLIERE_NU; ++k) {
        target = ((double) k / (double) (OSH_MOLIERE_NU - 1)) * total;
        while (idx < OSH_MOLIERE_NF - 1 && cdf[idx] < target) {
            ++idx;
        }
        seg = cdf[idx] - cdf[idx - 1];
        if (seg > 0.0) {
            row[k] = (float) (((double) (idx - 1) + (target - cdf[idx - 1]) / seg) * osh_moliere_dth);
        } else {
            row[k] = (float) ((double) idx * osh_moliere_dth);
        }
    }
}

void osh_physics_moliere_init(void) {
    double ymax;
    double dy;
    double th;
    double s0;
    double s1;
    double s2;
    double y;
    double e;
    double base;
    double q;
    double t1;
    double bnode;
    int my;
    int i;
    int j;
    int jb;

    if (osh_moliere_inited) {
        return;
    }

    ymax = 12.0;
    my = 4000;
    dy = ymax / (double) my;
    osh_moliere_dth = OSH_MOLIERE_THMAX / (double) (OSH_MOLIERE_NF - 1);

    for (i = 0; i < OSH_MOLIERE_NF; ++i) {
        th = (double) i * osh_moliere_dth;
        s0 = 0.0;
        s1 = 0.0;
        s2 = 0.0;
        for (j = 1; j <= my; ++j) {
            y = (double) j * dy;
            e = exp(-0.25 * y * y);
            base = y * osh_bessel_j0(th * y) * e;
            q = 0.25 * y * y;
            t1 = q * log(q);
            s0 += base;
            s1 += base * t1;
            s2 += base * t1 * t1;
        }
        osh_moliere_f0[i] = s0 * dy;
        osh_moliere_f1[i] = s1 * dy;
        osh_moliere_f2[i] = 0.5 * s2 * dy;
    }

    /* Precompute the inverse-CDF support table over a linear B grid. */
    osh_moliere_inv_db = (double) (OSH_MOLIERE_NB - 1) / (OSH_MOLIERE_B_MAX - OSH_MOLIERE_B_MIN);
    for (jb = 0; jb < OSH_MOLIERE_NB; ++jb) {
        bnode = OSH_MOLIERE_B_MIN + (double) jb / osh_moliere_inv_db;
        moliere_build_invcdf_row(bnode, osh_moliere_invcdf[jb]);
    }

    osh_moliere_inited = 1;
}

int osh_physics_moliere_solve_b(double omega, double *b_out) {
    double cnst;
    double b;
    double f;
    double fp;
    double db;
    int it;

    if (omega <= 0.0) {
        return 0;
    }

    cnst = log(omega);
    if (cnst <= 1.0816) {
        return 0;
    }

    b = cnst + log(cnst);
    for (it = 0; it < 50; ++it) {
        f = b - log(b) - cnst;
        fp = 1.0 - 1.0 / b;
        if (fp == 0.0) {
            break;
        }
        db = -f / fp;
        b += db;
        if (b < cnst) {
            b = cnst;
        }
        if (fabs(db) < 1e-8) {
            break;
        }
    }

    if (b_out) {
        *b_out = b;
    }
    return 1;
}

double osh_physics_moliere_reduced_f(int n, double theta_reduced) {
    double const *tab;
    double g;
    double frac;
    int i;

    if (!osh_moliere_inited) {
        osh_physics_moliere_init();
    }
    if (theta_reduced < 0.0 || theta_reduced > OSH_MOLIERE_THMAX) {
        return 0.0;
    }

    switch (n) {
    case 0:
        tab = osh_moliere_f0;
        break;
    case 1:
        tab = osh_moliere_f1;
        break;
    case 2:
        tab = osh_moliere_f2;
        break;
    default:
        return 0.0;
    }

    g = theta_reduced / osh_moliere_dth;
    i = (int) g;
    if (i >= OSH_MOLIERE_NF - 1) {
        return tab[OSH_MOLIERE_NF - 1];
    }
    frac = g - (double) i;
    return tab[i] * (1.0 - frac) + tab[i + 1] * frac;
}

double osh_physics_moliere_inv_cdf(double b, double u) {
    double gb;
    double fb;
    double gu;
    double fu;
    double r0;
    double r1;
    int jb;
    int ku;

    if (!osh_moliere_inited) {
        osh_physics_moliere_init();
    }
    if (b <= 0.0) {
        return 0.0;
    }

    /* Clamp B into the tabulated range and locate the B bracket (linear grid). */
    if (b < OSH_MOLIERE_B_MIN) {
        b = OSH_MOLIERE_B_MIN;
    } else if (b > OSH_MOLIERE_B_MAX) {
        b = OSH_MOLIERE_B_MAX;
    }
    gb = (b - OSH_MOLIERE_B_MIN) * osh_moliere_inv_db;
    jb = (int) gb;
    if (jb >= OSH_MOLIERE_NB - 1) {
        jb = OSH_MOLIERE_NB - 2;
    }
    fb = gb - (double) jb;

    /* Clamp u to [0,1] and locate the u bracket (uniform grid). */
    if (u < 0.0) {
        u = 0.0;
    } else if (u > 1.0) {
        u = 1.0;
    }
    gu = u * (double) (OSH_MOLIERE_NU - 1);
    ku = (int) gu;
    if (ku >= OSH_MOLIERE_NU - 1) {
        ku = OSH_MOLIERE_NU - 2;
    }
    fu = gu - (double) ku;

    /* Bilinear interpolation of the inverse-CDF table in (B, u). */
    r0 = (double) osh_moliere_invcdf[jb][ku] * (1.0 - fu) + (double) osh_moliere_invcdf[jb][ku + 1] * fu;
    r1 = (double) osh_moliere_invcdf[jb + 1][ku] * (1.0 - fu) + (double) osh_moliere_invcdf[jb + 1][ku + 1] * fu;
    return r0 * (1.0 - fb) + r1 * fb;
}

double osh_physics_moliere_sample_reduced(double b, struct osh_rng *rng) {
    return osh_physics_moliere_inv_cdf(b, osh_rng_double(rng));
}

int osh_physics_moliere_scatter(double const v[3],
                                double w[3],
                                double t_kin_mev,
                                double mass_mev,
                                double z_eff,
                                double d_gcm2,
                                double path_scale_gcm2,
                                double chic2_coeff,
                                double screen_z,
                                struct osh_rng *rng) {
    double p2;
    double e_total;
    double pv;
    double beta2;
    double beta;
    double chic2_step;
    double chic2_path;
    double alpha;
    double chi_a2;
    double omega;
    double b;
    double theta;
    double st;
    double ct;
    double tr;
    double th;
    double cos_phi;
    double sin_phi;
    int tries;

    w[0] = v[0];
    w[1] = v[1];
    w[2] = v[2];

    if (chic2_coeff <= 0.0 || z_eff <= 0.0 || screen_z <= 0.0 || d_gcm2 <= 0.0) {
        return 0;
    }
    if (path_scale_gcm2 <= 0.0) {
        path_scale_gcm2 = d_gcm2;
    }

    p2 = t_kin_mev * (t_kin_mev + 2.0 * mass_mev);
    e_total = t_kin_mev + mass_mev;
    if (p2 <= 0.0 || e_total <= 0.0) {
        return 0;
    }

    pv = p2 / e_total;
    beta2 = p2 / (e_total * e_total);
    beta = sqrt(beta2);
    if (beta <= 0.0) {
        return 0;
    }

    chic2_step = chic2_coeff * z_eff * z_eff * d_gcm2 / (pv * pv);
    chic2_path = chic2_coeff * z_eff * z_eff * path_scale_gcm2 / (pv * pv);

    alpha = z_eff * screen_z * OSH_ALPHA_FS / beta;
    chi_a2 = 2.007e-5 * pow(screen_z, 2.0 / 3.0) * (1.0 + 3.34 * alpha * alpha) / p2;
    if (chi_a2 <= 0.0) {
        return 0;
    }

    omega = chic2_path / (1.167 * chi_a2);
    if (!osh_physics_moliere_solve_b(omega, &b)) {
        return 0;
    }

    theta = 0.0;
    for (tries = 0; tries < 16; ++tries) {
        tr = osh_physics_moliere_sample_reduced(b, rng);
        th = sqrt(chic2_step * b) * tr;
        if (th >= OSH_M_PI) {
            continue;
        }
        if (th < 1e-6 || osh_rng_double(rng) * th <= sin(th)) {
            theta = th;
            break;
        }
    }
    if (theta <= 0.0) {
        return 0;
    }

    st = sin(theta);
    ct = cos(theta);
    osh_kinematics_azimuth(rng, &cos_phi, &sin_phi);
    osh_kinematics_rotate_dir_cos(v, w, ct, st, cos_phi, sin_phi);
    return 1;
}
