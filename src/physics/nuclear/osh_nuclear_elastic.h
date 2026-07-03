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
 *   - sigma_el = sigma_tot - sigma_reac, with sigma_tot from a black-disk optical
 *     limit (sigma_tot = 2*sigma_reac), so sigma_el ~= sigma_reac (Tripathi).
 *     A single tunable prefactor scales it.
 *   - the angular distribution is nuclear diffraction, dsigma/dt ~ exp(B(A)*t),
 *     forward-peaked with slope B(A) set by the nuclear radius R = r0*A^(1/3);
 *     the recoil nucleus therefore carries a small (sub-MeV) kinetic energy.
 *
 * The scattered projectile stays a primary (elastic keeps the generation, FLUKA
 * convention); the recoil nucleus is emitted as a secondary.
 */

/** Elastic cross section [cm^2] for projectile (zp,ap) on target (zt,at). */
double osh_nuclear_elastic_sigma(unsigned int zp, unsigned int ap, double zt, double at, double e_lab_per_nucleon);

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
