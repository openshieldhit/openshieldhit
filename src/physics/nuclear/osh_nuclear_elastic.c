#include "physics/nuclear/osh_nuclear_elastic.h"

#include <math.h>

#include "physics/nuclear/osh_nuclear_tripathi.h"
#include "random/osh_rng.h"

/* sigma_el / sigma_reac prefactor.  Black-disk optical limit: sigma_tot =
 * 2*sigma_reac, so sigma_el = sigma_tot - sigma_reac = sigma_reac (factor 1).
 * Tunable at compile time for closer agreement with data/SH12A. */
#ifndef OSH_NUCLEAR_ELASTIC_SIGMA_FACTOR
#define OSH_NUCLEAR_ELASTIC_SIGMA_FACTOR 1.0
#endif

/* Nuclear radius parameter [fm]: R = r0 * A^(1/3), used for the diffraction slope. */
#ifndef OSH_NUCLEAR_ELASTIC_R0_FM
#define OSH_NUCLEAR_ELASTIC_R0_FM 1.3
#endif

/* hbar*c [MeV fm] — converts the slope from fm^2 to (MeV/c)^-2. */
#define OSH_HBARC_MEV_FM 197.327

double osh_nuclear_elastic_sigma(unsigned int zp, unsigned int ap, double zt, double at, double e_lab_per_nucleon) {
    double sigma_reac = osh_nuclear_tripathi_sigma(zp, ap, zt, at, e_lab_per_nucleon);
    if (!(sigma_reac > 0.0)) {
        return 0.0;
    }
    return OSH_NUCLEAR_ELASTIC_SIGMA_FACTOR * sigma_reac;
}

double osh_nuclear_elastic_lambda_gcm2(double at_g_per_mol, double sigma_cm2) {
    return osh_nuclear_lambda_gcm2(at_g_per_mol, sigma_cm2);
}

double osh_nuclear_elastic_slope(double at) {
    /* Diffraction slope of a disk of radius R: B ~ R^2/3.  R = r0*A^(1/3) [fm];
     * divide by (hbar c)^2 to express B in (MeV/c)^-2. */
    double r = OSH_NUCLEAR_ELASTIC_R0_FM * cbrt(at);
    double b_fm2 = (r * r) / 3.0;
    return b_fm2 / (OSH_HBARC_MEV_FM * OSH_HBARC_MEV_FM);
}

double osh_nuclear_elastic_sample_cos_cm(double p_cm, double at, struct osh_rng *rng) {
    double b;      /* diffraction slope [ (MeV/c)^-2 ] */
    double t_max;  /* |t| at 180 deg CM = 4 p_cm^2      */
    double norm;   /* CDF normalisation over [0, t_max] */
    double t;      /* sampled momentum transfer |t|     */
    double cos_cm;

    if (!(p_cm > 0.0)) {
        return 1.0;
    }
    b = osh_nuclear_elastic_slope(at);
    if (!(b > 0.0)) {
        return 1.0;
    }
    t_max = 4.0 * p_cm * p_cm;
    /* |t| from the truncated exponential dsigma/d|t| ~ exp(-B|t|) on [0, t_max]:
     * inverse-CDF sampling, t = -ln(1 - u*(1-exp(-B t_max)))/B. */
    norm = 1.0 - exp(-b * t_max);
    if (!(norm > 0.0)) {
        return 1.0;
    }
    t = -log(1.0 - (osh_rng_double(rng) * norm)) / b;
    cos_cm = 1.0 - (t / (2.0 * p_cm * p_cm));
    if (cos_cm < -1.0) {
        cos_cm = -1.0;
    }
    if (cos_cm > 1.0) {
        cos_cm = 1.0;
    }
    return cos_cm;
}
