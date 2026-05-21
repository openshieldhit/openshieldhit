#include "physics/nuclear/osh_nuclear_tripathi.h"

#include <math.h>

#include "openshieldhit/const.h"

/*
 * Nuclear radius parameter [cm].  Tripathi 1999 uses r0 = 1.1 fm.
 */
#define OSH_TRIPATHI_R0_CM 1.1e-13

/*
 * Coulomb radius parameter [fm] used in the barrier calculation.
 * Tripathi 1999 eq. (4): r_p = r_t = 1.29 * A^(1/3) fm.
 */
#define OSH_TRIPATHI_RC_FM 1.29

/*
 * Coulomb constant [MeV·fm]: e²/(4πε₀) ≈ 1.44 MeV·fm.
 */
#define OSH_COULOMB_MEV_FM 1.44

double osh_nuclear_tripathi_sigma(unsigned int zp, unsigned int ap, double zt, double at, double e_lab_per_nucleon) {
    double ap3;
    double at3;
    double rp_fm;
    double rt_fm;
    double bc;
    double ecm;
    double s;
    double d;
    double ce;
    double delta_e;
    double arg;
    double sigma;

    if (ap == 0u || at <= 0.0 || e_lab_per_nucleon <= 0.0)
        return 0.0;

    ap3 = cbrt((double) ap);
    at3 = cbrt(at);

    /* Coulomb barrier [MeV] — eq. (4) in Tripathi 1999 */
    rp_fm = OSH_TRIPATHI_RC_FM * ap3;
    rt_fm = OSH_TRIPATHI_RC_FM * at3;
    bc = OSH_COULOMB_MEV_FM * (double) zp * zt / (rp_fm + rt_fm);

    /* Centre-of-mass kinetic energy per nucleon [MeV/nucleon] — non-relativistic */
    ecm = e_lab_per_nucleon * at / ((double) ap + at);

    if (ecm <= bc)
        return 0.0;

    /* Nuclear overlap parameter S — eq. (2) */
    s = ap3 * at3 / (ap3 + at3);

    /* Pauli blocking coefficient D — proton projectile uses 2.05, others 1.75 */
    d = (ap == 1u) ? 2.05 : 1.75;

    /* Energy-dependent correction CE — eq. (3) */
    ce = d * (1.0 - exp(-ecm / 40.0));

    /* Transparency and isospin correction delta_E — eq. (1) */
    delta_e = 1.85 * s + 0.16 * s / cbrt(ecm) - ce + 0.91 * (at - 2.0 * zt) * (double) zp / (at * (double) ap);

    arg = ap3 + at3 + delta_e;
    if (arg <= 0.0)
        return 0.0;

    sigma = OSH_M_PI * OSH_TRIPATHI_R0_CM * OSH_TRIPATHI_R0_CM * arg * arg * (1.0 - bc / ecm);

    return (sigma > 0.0) ? sigma : 0.0;
}

double osh_nuclear_lambda_gcm2(double at_g_per_mol, double sigma_cm2) {
    if (sigma_cm2 <= 0.0)
        return 1.0e30;
    return at_g_per_mol / (OSH_NAVOGADRO * sigma_cm2);
}

double osh_nuclear_survival_prob(double ds_gcm2, double lambda_gcm2) {
    if (ds_gcm2 <= 0.0)
        return 1.0;
    if (lambda_gcm2 <= 0.0)
        return 0.0;
    return exp(-ds_gcm2 / lambda_gcm2);
}
