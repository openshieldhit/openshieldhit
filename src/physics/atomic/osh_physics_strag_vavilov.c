#include "physics/atomic/osh_physics_strag_vavilov.h"

#include <math.h>
#include <stddef.h>

#include "physics/atomic/osh_physics_strag_vavilov_coeffs.h"

/* Per-u-band transform x'(u); ids match emit_coeffs.py (0=-log u, 1=u, 2=-log(1-u)). */
static double _vav_xform(double u, int kind) {
    switch (kind) {
    case 0:
        return -log(u);
    case 1:
        return u;
    default:
        return -log1p(-u);
    }
}

double osh_physics_strag_vavilov_lambda(double kappa, double beta2, double u) {
    int b;                  /* selected κ-band index (0..NKB-1)                    */
    int ub;                 /* selected u-band index (0..NUB-1)                    */
    int k;                  /* polynomial power in the Horner loop                 */
    int i;                  /* Chebyshev-product / coefficient index (0..NT-1)     */
    int j;                  /* β Chebyshev order index (0..NBE-1)                  */
    double kn;              /* κ mapped to [-1,1]: normalised ln κ within the band */
    double bn;              /* β² mapped to [-1,1] over [BLO, BHI]                 */
    double xn;              /* transformed u mapped to ~[-1,1] within the u-band   */
    double lam;             /* result: reduced Vavilov variable λ                  */
    double tk[OSH_VAV_NKA]; /* Chebyshev polynomials T_i(kn), i = 0..NKA-1         */
    double tb[OSH_VAV_NBE]; /* Chebyshev polynomials T_j(bn), j = 0..NBE-1         */
    double t[OSH_VAV_NT];   /* Chebyshev product basis, t[i*NBE + j] = T_i·T_j     */
    double const *creg;     /* base of this (κ-band, u-band) region's coefficients */
    double const *ck;       /* coefficient row for the current polynomial power k  */
    double ak;              /* a_k(κ,β²) = Σ_i t[i]·c[k][i] — the k-th u-coefficient */

    if (u < OSH_VAV_UMIN) {
        u = OSH_VAV_UMIN;
    }
    if (u > OSH_VAV_UMAX) {
        u = OSH_VAV_UMAX;
    }

    /* κ-band dispatch (last band catches all remaining κ). */
    for (b = 0; b < OSH_VAV_NKB - 1; ++b) {
        if (kappa < osh_vav_kdisp_hi[b]) {
            break;
        }
    }
    kn = (2.0 * log(kappa) - log(osh_vav_klo[b]) - log(osh_vav_khi[b])) / (log(osh_vav_khi[b]) - log(osh_vav_klo[b]));

    /* u-band dispatch. */
    for (ub = 0; ub < OSH_VAV_NUB - 1; ++ub) {
        if (u <= osh_vav_uhi[ub]) {
            break;
        }
    }
    xn = (_vav_xform(u, osh_vav_ukind[ub]) - osh_vav_umid[ub]) / osh_vav_uhalf[ub];

    bn = (2.0 * beta2 - OSH_VAV_BLO - OSH_VAV_BHI) / (OSH_VAV_BHI - OSH_VAV_BLO);

    /* Chebyshev product basis, same order as the fit: t[i*NBE + j]. */
    tk[0] = 1.0;
    tk[1] = kn;
    tk[2] = 2.0 * kn * kn - 1.0;
    tk[3] = (4.0 * kn * kn - 3.0) * kn;
    tb[0] = 1.0;
    tb[1] = bn;
    tb[2] = 2.0 * bn * bn - 1.0;
    for (i = 0; i < OSH_VAV_NKA; ++i) {
        for (j = 0; j < OSH_VAV_NBE; ++j) {
            t[(i * OSH_VAV_NBE) + j] = tk[i] * tb[j];
        }
    }

    /* λ = Horner in xn over powers 0..udeg; each coefficient a_k = Σ_i t[i]·c[k][i]. */
    creg = &osh_vav_coef[(((size_t) b * OSH_VAV_NUB) + ub) * OSH_VAV_MAXP * OSH_VAV_NT];
    lam = 0.0;
    for (k = osh_vav_udeg[ub]; k >= 0; --k) {
        ak = 0.0;
        ck = creg + (size_t) k * OSH_VAV_NT;
        for (i = 0; i < OSH_VAV_NT; ++i) {
            ak += t[i] * ck[i];
        }
        lam = lam * xn + ak;
    }
    return lam;
}
