#include "physics/atomic/osh_physics_scat_moliere.h"

#include <math.h>

#include "openshieldhit/const.h"
#include "physics/osh_kinematics.h"
#include "random/osh_rng.h"

#define OSH_MOLIERE_NF 257     /* reduced-angle (ϑ) grid points for f⁽ⁿ⁾ */
#define OSH_MOLIERE_THMAX 12.0 /* maximum reduced angle ϑ */

/* Inverse-CDF support table: ϑ as a function of (B, u).  The reduced-angle
 * density depends on B only through f(ϑ) = f⁽⁰⁾ + f⁽¹⁾/B + f⁽²⁾/B², so ϑ(u) is a
 * smooth function of (B, u) that we tabulate once at init and bilinearly
 * interpolate on the hot path — no per-scatter CDF reconstruction.  Built from
 * our own f⁽ⁿ⁾ formulas (no external tables). */
#define OSH_MOLIERE_NB 64      /* B-grid nodes */
#define OSH_MOLIERE_NU 1024    /* u-grid nodes (cumulative probability) */
#define OSH_MOLIERE_B_MIN 1.6  /* ≈ solve_b threshold (ln Ω = 1.0816 ⇒ B ≈ 1.6) */
#define OSH_MOLIERE_B_MAX 40.0 /* above this the 1/B, 1/B² corrections are negligible */

static double osh_moliere_f0[OSH_MOLIERE_NF];
static double osh_moliere_f1[OSH_MOLIERE_NF];
static double osh_moliere_f2[OSH_MOLIERE_NF];
static double osh_moliere_dth;
static float osh_moliere_invcdf[OSH_MOLIERE_NB][OSH_MOLIERE_NU];
static double osh_moliere_m2[OSH_MOLIERE_NB]; /* ⟨ϑ²⟩ per B node (incl. Rutherford tail) */
static double osh_moliere_inv_db;             /* (NB-1) / (B_MAX - B_MIN) */
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
 * @param[in]  b      Molière B parameter for this row.
 * @param[out] row    Destination of NU interpolated ϑ values.
 * @param[out] m2_out ⟨ϑ²⟩ of the (clamped, truncated) distribution for this B.
 */
static void moliere_build_invcdf_row(double b, float *row, double *m2_out) {
    double cdf[OSH_MOLIERE_NF]; /* cumulative of ϑ·f(ϑ;B) over the ϑ grid       */
    double inv_b;               /* 1/B                                           */
    double inv_b2;              /* 1/B²                                          */
    double acc;                 /* running CDF accumulator (∫ ϑ·f dϑ)            */
    double acc_m2;              /* running ∫ ϑ²·(ϑ·f) dϑ for ⟨ϑ²⟩               */
    double total;               /* CDF normalisation = ∫ ϑ·f dϑ                  */
    double th;                  /* reduced angle ϑ at the current grid point     */
    double pdf;                 /* sampling density ϑ·f(ϑ;B) (clamped ≥ 0)       */
    double target;              /* u_k·total: CDF value to invert for u node k    */
    double seg;                 /* CDF span of the bracketing ϑ interval         */
    int i;                      /* ϑ-grid index                                  */
    int k;                      /* u-grid (cumulative-probability) index         */
    int idx;                    /* walking ϑ-grid cursor during inversion        */

    inv_b = 1.0 / b;
    inv_b2 = inv_b * inv_b;
    acc = 0.0;
    acc_m2 = 0.0;
    for (i = 0; i < OSH_MOLIERE_NF; ++i) {
        th = (double) i * osh_moliere_dth;
        pdf = th * (osh_moliere_f0[i] + osh_moliere_f1[i] * inv_b + osh_moliere_f2[i] * inv_b2);
        if (pdf < 0.0) {
            pdf = 0.0; /* clamp rare negative dips so the CDF stays monotone */
        }
        acc += pdf * osh_moliere_dth;
        acc_m2 += th * th * pdf * osh_moliere_dth;
        cdf[i] = acc;
    }

    total = cdf[OSH_MOLIERE_NF - 1];
    /* ⟨ϑ²⟩ over the same clamped/truncated density that is actually sampled. */
    *m2_out = (total > 0.0) ? (acc_m2 / total) : 1.0;
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
    double ymax;  /* upper limit of the y (Hankel) integration                 */
    double dy;    /* y quadrature step                                         */
    double th;    /* current reduced angle ϑ on the f⁽ⁿ⁾ grid                  */
    double s0;    /* running ∫ for f⁽⁰⁾(ϑ)                                     */
    double s1;    /* running ∫ for f⁽¹⁾(ϑ)                                     */
    double s2;    /* running ∫ for f⁽²⁾(ϑ) (before the 1/2! factor)            */
    double y;     /* Hankel integration variable                              */
    double e;     /* Gaussian damping e^(−y²/4)                               */
    double base;  /* y·J₀(ϑy)·e^(−y²/4): common factor of all three integrands */
    double q;     /* y²/4                                                      */
    double t1;    /* (y²/4)·ln(y²/4): the per-order weight                     */
    double bnode; /* B value at the current inverse-CDF table node            */
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
        moliere_build_invcdf_row(bnode, osh_moliere_invcdf[jb], &osh_moliere_m2[jb]);
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
    double gb; /* fractional B-grid coordinate (b − B_MIN)·inv_db      */
    double fb; /* interpolation weight within the B bracket            */
    double gu; /* fractional u-grid coordinate u·(NU−1)                */
    double fu; /* interpolation weight within the u bracket            */
    double r0; /* ϑ interpolated in u at the lower B node              */
    double r1; /* ϑ interpolated in u at the upper B node              */
    int jb;    /* lower B-grid index                                   */
    int ku;    /* lower u-grid index                                   */

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

/* ⟨ϑ²⟩ of the reduced Molière distribution at B (linear interp on the B grid). */
static double osh_moliere_mean_sq(double b) {
    double gb;
    double fb;
    int jb;

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
    return osh_moliere_m2[jb] * (1.0 - fb) + osh_moliere_m2[jb + 1] * fb;
}

int osh_physics_moliere_scatter(double const v[3],
                                double w[3],
                                double theta0,
                                double t_kin_mev,
                                double mass_mev,
                                double z_eff,
                                double d_gcm2,
                                double path_scale_gcm2,
                                double chic2_coeff,
                                double screen_z,
                                struct osh_rng *rng) {
    double p2;         /* momentum² p² = T(T+2m)                       [MeV²]   */
    double e_total;    /* total energy E = T + m                        [MeV]    */
    double pv;         /* βcp = p²/E                                     [MeV]    */
    double beta2;      /* β² = p²/E²                                     [—]      */
    double beta;       /* velocity β = v/c                              [—]      */
    double chic2_path; /* χ_c² over the macroscopic path scale          [rad²]   */
    double alpha;      /* Born parameter α = z·Z·α_fs/β (screening)      [—]      */
    double chi_a2;     /* Molière screening angle squared χ_a²          [rad²]   */
    double omega;      /* Ω = χ_c²/(1.167·χ_a²): effective # of scatters [—]      */
    double b;          /* Molière B: solves B − ln B = ln Ω             [—]      */
    double scale;      /* θ = scale·ϑ; anchors RMS to Highland θ₀        [rad]    */
    double theta;      /* sampled space polar deflection angle          [rad]    */
    double st;         /* sin(theta)                                              */
    double ct;         /* cos(theta)                                              */
    double tr;         /* sampled reduced angle ϑ (Molière shape)       [—]      */
    double th;         /* trial physical angle scale·ϑ before rejection [rad]    */
    double cos_phi;    /* azimuth cosine (uniform φ)                              */
    double sin_phi;    /* azimuth sine                                            */
    int tries;         /* sin θ/θ rejection-loop counter                          */

    /* Default: copy the incident direction (used on every early-out path). */
    w[0] = v[0];
    w[1] = v[1];
    w[2] = v[2];

    if (theta0 <= 0.0 || chic2_coeff <= 0.0 || z_eff <= 0.0 || screen_z <= 0.0 || d_gcm2 <= 0.0) {
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

    /* B (the distribution-shape parameter) from χ_c² over the macroscopic path
     * and the screening angle χ_a (Lynch & Dahl). */
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

    /* Anchor the *magnitude* to the validated Highland width: the Molière
     * distribution sets only the shape (Gaussian core + Rutherford tail), and
     * we rescale so the sampled space angle has variance 2·θ₀² — i.e. projected
     * RMS = θ₀ — using the precomputed ⟨ϑ²⟩(B).  This keeps mode 2's width equal
     * to mode 1 (and to data) while retaining the tail. */
    scale = theta0 * sqrt(2.0 / osh_moliere_mean_sq(b));

    theta = 0.0;
    for (tries = 0; tries < 16; ++tries) {
        tr = osh_physics_moliere_sample_reduced(b, rng);
        th = scale * tr;
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
