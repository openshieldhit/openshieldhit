#ifndef OSH_NUCLEAR_ELASTIC_H
#define OSH_NUCLEAR_ELASTIC_H

#ifdef __cplusplus
extern "C" {
#endif

struct osh_rng;

/**
 * @file osh_nuclear_elastic.h
 * @brief Proton-nucleus (p+A) elastic scattering: cross section, mean free path,
 *        and diffraction-model recoil-angle sampling.
 *
 * @details
 * Approximate model for issue #212 (projectile-agnostic API so nucleus-nucleus
 * elastic can reuse it later):
 *   - sigma_el = ratio(E) * sigma_reac, with an energy-dependent ratio
 *     calibrated against the measured integrated nuclear elastic scale
 *     (Garron 1962, p+C-12 at 155 MeV) and the SH12A-implied removal rate
 *     (issue #277): black-disk scale (~1) at and below ~30 MeV, falling to
 *     ~0.32 above ~180 MeV, log-E interpolated between.  A tunable overall
 *     prefactor scales it.
 *   - the angular distribution is nuclear diffraction, dsigma/dt ~ exp(B(A)*t),
 *     forward-peaked with slope B(A) set by the nuclear radius R = r0*A^(1/3);
 *     the recoil nucleus therefore carries a small (sub-MeV) kinetic energy.
 *
 * The scattered projectile stays a primary (elastic keeps the generation, FLUKA
 * convention); the recoil nucleus is emitted as a secondary.
 */

/** Elastic cross section [cm^2] for projectile (zp,ap) on target (zt,at). */
double osh_nuclear_elastic_sigma(unsigned int zp, unsigned int ap, double zt, double at, double e_lab_per_nucleon);

/**
 * @brief Elastic cross section [cm^2] from an already-evaluated reaction cross
 *        section @p sigma_reac_cm2, applying the energy-dependent sigma_el /
 *        sigma_reac ratio at @p e_lab_per_nucleon [MeV/u].
 *
 * Lets a caller that already holds σ_R (e.g. the shared rate loop that feeds both
 * the inelastic and the p+A-elastic hazard) obtain σ_el without re-evaluating the
 * reaction cross section. `osh_nuclear_elastic_sigma()` is exactly this applied
 * to a fresh osh_nuclear_sigma_reac() call.
 */
double osh_nuclear_elastic_sigma_from_reac(double sigma_reac_cm2, double e_lab_per_nucleon);

/** Mean free path [g/cm^2] from target molar mass [g/mol] and sigma [cm^2]. */
double osh_nuclear_elastic_lambda_gcm2(double at_g_per_mol, double sigma_cm2);

/** Diffraction slope B [ (MeV/c)^-2 ] for a target of mass number @p at. */
double osh_nuclear_elastic_slope(double at);

/**
 * @brief Sample the CM scattering-angle cosine from dsigma/dt ~ exp(B(A)*t).
 * @param p_cm  CM momentum magnitude [MeV/c].
 * @param at    target mass number (sets the slope).
 * @param rng   RNG state (mutated).
 * @returns cos(theta_CM) of the scattered projectile in [-1, 1].
 */
double osh_nuclear_elastic_sample_cos_cm(double p_cm, double at, struct osh_rng *rng);

#ifdef __cplusplus
}
#endif

#endif /* OSH_NUCLEAR_ELASTIC_H */
