#include "physics/nuclear/osh_nuclear_elastic.h"

#include <math.h>

#include "physics/nuclear/osh_nuclear_sigma_reac.h"
#include "physics/nuclear/osh_nuclear_tripathi.h" /* osh_nuclear_lambda_gcm2 */
#include "random/osh_rng.h"

/* Overall sigma_el scale on top of the energy-dependent ratio below.
 * Tunable at compile time for closer agreement with data/SH12A. */
#ifndef OSH_NUCLEAR_ELASTIC_SIGMA_FACTOR
#define OSH_NUCLEAR_ELASTIC_SIGMA_FACTOR 1.0
#endif

/*
 * Energy-dependent sigma_el / sigma_reac ratio (issue #277).
 *
 * The former constant black-disk factor (sigma_el = sigma_reac) overpredicts
 * the integrated nuclear elastic cross section by a factor 2-3 above
 * ~100 MeV: the only integrated (P,EL) measurement on a therapy target,
 * Garron 1962 (p+C-12 at 155 MeV, 75 +- 7 mb, EXFOR O0360005), gives
 * sigma_el / sigma_reac = 0.34 against the LA150 sigma_reac of 222 mb.
 * Toward low energies optical-model systematics bring the ratio up to the
 * black-disk scale (~1 at and below ~30 MeV).
 *
 * Model: ratio = RATIO_LOW for E <= E_LOW, RATIO_HIGH for E >= E_HIGH, and
 * log-E linear interpolation in between.  With the defaults the Garron point
 * evaluates to 0.38 * 222 mb = 84 mb — within 1.3 sigma of the measurement.
 * The transport-level acceptance (primary depth fluence vs SH12A, issue
 * #277) is the final calibration of these anchors.
 */
#ifndef OSH_NUCLEAR_ELASTIC_RATIO_LOW
#define OSH_NUCLEAR_ELASTIC_RATIO_LOW 1.0
#endif
#ifndef OSH_NUCLEAR_ELASTIC_RATIO_HIGH
#define OSH_NUCLEAR_ELASTIC_RATIO_HIGH 0.32
#endif
#ifndef OSH_NUCLEAR_ELASTIC_RATIO_E_LOW_MEV
#define OSH_NUCLEAR_ELASTIC_RATIO_E_LOW_MEV 30.0
#endif
#ifndef OSH_NUCLEAR_ELASTIC_RATIO_E_HIGH_MEV
#define OSH_NUCLEAR_ELASTIC_RATIO_E_HIGH_MEV 180.0
#endif

/* Nuclear radius parameter [fm]: R = r0 * A^(1/3), used for the diffraction slope. */
#ifndef OSH_NUCLEAR_ELASTIC_R0_FM
#define OSH_NUCLEAR_ELASTIC_R0_FM 1.3
#endif

/* hbar*c [MeV fm] — converts the slope from fm^2 to (MeV/c)^-2. */
#define OSH_HBARC_MEV_FM 197.327

/** sigma_el / sigma_reac ratio at lab energy @p e_mev (log-E interpolation
 * between the low- and high-energy anchors). */
static inline double _elastic_ratio(double e_mev) {
    double frac;

    if (e_mev <= OSH_NUCLEAR_ELASTIC_RATIO_E_LOW_MEV) {
        return OSH_NUCLEAR_ELASTIC_RATIO_LOW;
    }
    if (e_mev >= OSH_NUCLEAR_ELASTIC_RATIO_E_HIGH_MEV) {
        return OSH_NUCLEAR_ELASTIC_RATIO_HIGH;
    }
    frac = log(e_mev / OSH_NUCLEAR_ELASTIC_RATIO_E_LOW_MEV)
           / log(OSH_NUCLEAR_ELASTIC_RATIO_E_HIGH_MEV / OSH_NUCLEAR_ELASTIC_RATIO_E_LOW_MEV);
    return OSH_NUCLEAR_ELASTIC_RATIO_LOW + (frac * (OSH_NUCLEAR_ELASTIC_RATIO_HIGH - OSH_NUCLEAR_ELASTIC_RATIO_LOW));
}

double osh_nuclear_elastic_sigma_from_reac(double sigma_reac_cm2, double e_lab_per_nucleon) {
    if (!(sigma_reac_cm2 > 0.0)) {
        return 0.0;
    }
    return OSH_NUCLEAR_ELASTIC_SIGMA_FACTOR * _elastic_ratio(e_lab_per_nucleon) * sigma_reac_cm2;
}

double osh_nuclear_elastic_sigma(unsigned int zp, unsigned int ap, double zt, double at, double e_lab_per_nucleon) {
    return osh_nuclear_elastic_sigma_from_reac(osh_nuclear_sigma_reac(zp, ap, zt, at, e_lab_per_nucleon),
                                               e_lab_per_nucleon);
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
    double b;     /* diffraction slope [ (MeV/c)^-2 ] */
    double t_max; /* |t| at 180 deg CM = 4 p_cm^2      */
    double norm;  /* CDF normalisation over [0, t_max] */
    double t;     /* sampled momentum transfer |t|     */
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
