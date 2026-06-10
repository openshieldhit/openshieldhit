#ifndef OSH_NUCLEAR_PP_H
#define OSH_NUCLEAR_PP_H

/**
 * @file osh_nuclear_pp.h
 *
 * @brief Proton-proton elastic scattering: cross section and angle sampling.
 *
 * @details
 * Data source: G4LEppData.hh (Arndt PSA 1998 + Jones 2010 extension).
 * Valid lab kinetic energy range: 10 MeV – 5 GeV.
 *
 * All hot-path functions are trig-free: cos(θ_CM) is returned directly from the
 * CDF table without any arccos call.
 */

struct osh_rng;

/**
 * @brief pp elastic cross section [cm²] at lab kinetic energy @p e_lab_mev.
 *
 * @details
 * Linearly interpolates σ_tot × (1 − σ_in_frac) from the tabulated data.
 * Returns 0 for energies below 10 MeV; clamps to the last bin above 5 GeV.
 *
 * @param e_lab_mev  Lab kinetic energy [MeV].
 * @returns σ_el [cm²].
 */
double osh_nuclear_pp_sigma_el(double e_lab_mev);

/**
 * @brief Sample the CM scattering angle cosine from the tabulated differential
 *        cross section, without any arccos call.
 *
 * @details
 * Bilinearly interpolates the CDF in (energy, angle); uses a forward scan for the
 * angle bracket to handle possible flat CDF segments (zero dσ/dΩ bins with equal
 * consecutive CDF entries), which would violate the strictly-increasing precondition
 * of osh_binary_search_f().  Consumes exactly one uniform RNG deviate.
 *
 * @param e_lab_mev  Lab kinetic energy [MeV].
 * @param rng        RNG state (consumes 1 deviate).
 * @returns cos(θ_CM) ∈ [cos(179.5°), cos(0.5°)].
 */
double osh_nuclear_pp_sample_cos_theta_cm(double e_lab_mev, struct osh_rng *rng);

/**
 * @brief pp elastic mean free path in a material with given hydrogen mass fraction.
 *
 * @details
 * λ = A_H / (f_H × N_A × σ_el)  [g/cm²]
 * where A_H ≈ 1.00794 g/mol and f_H is the hydrogen mass fraction.
 * Returns 1e30 when @p hydrogen_mass_fraction or @p sigma_el_cm2 is ≤ 0.
 *
 * @param hydrogen_mass_fraction  Sum of mass fractions for Z=1 elements [0, 1].
 * @param sigma_el_cm2            pp elastic cross section [cm²].
 * @returns Mean free path [g/cm²].
 */
double osh_nuclear_pp_lambda_gcm2(double hydrogen_mass_fraction, double sigma_el_cm2);

#endif /* OSH_NUCLEAR_PP_H */
