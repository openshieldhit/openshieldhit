#ifndef OSH_NUCLEAR_TRIPATHI_H
#define OSH_NUCLEAR_TRIPATHI_H

/**
 * @file osh_nuclear_tripathi.h
 *
 * @brief Tripathi total nuclear reaction cross section.
 *
 * @details
 * Implements the parametric total nuclear reaction cross section from:
 *
 *   Tripathi R.K., Cucinotta F.A., Wilson J.W.
 *   "Universal Parameterization of Absorption Cross Sections"
 *   NASA Technical Paper NASA/TP-1999-209726 (1999).
 *
 * The formula is energy-dependent and includes a natural Coulomb threshold,
 * making it accurate across the full clinical energy range.  Compound target
 * materials are handled via effective (non-integer) Z and A values.
 *
 * All functions are pure: no global state, no allocation.
 */

/**
 * @brief Total nuclear reaction cross section (Tripathi 1999).
 *
 * @param zp           Projectile atomic number.
 * @param ap           Projectile mass number.
 * @param zt           Target effective atomic number (may be non-integer for compounds).
 * @param at           Target effective mass number [g/mol] (may be non-integer).
 * @param e_lab_per_nucleon  Lab-frame kinetic energy per nucleon T/A [MeV/nucleon].
 *                          This is the total kinetic energy divided by the integer
 *                          mass number A — NOT energy per atomic mass unit (MeV/u),
 *                          which would differ by ~0.7 % for protons.
 *
 * @returns Cross section [cm²], or 0 if below the Coulomb threshold.
 */
double osh_nuclear_tripathi_sigma(unsigned int zp, unsigned int ap, double zt, double at, double e_lab_per_nucleon);

/**
 * @brief Nuclear interaction mean free path.
 *
 * @param at_g_per_mol  Target effective mass number [g/mol].
 * @param sigma_cm2     Total nuclear reaction cross section [cm²].
 *
 * @returns Mean free path [g/cm²].  Returns a large value when sigma_cm2 <= 0.
 */
double osh_nuclear_lambda_gcm2(double at_g_per_mol, double sigma_cm2);

/**
 * @brief Survival probability over an areal density step.
 *
 * @param ds_gcm2      Step areal density [g/cm²].
 * @param lambda_gcm2  Mean free path [g/cm²].
 *
 * @returns exp(-ds_gcm2 / lambda_gcm2), in [0, 1].
 */
double osh_nuclear_survival_prob(double ds_gcm2, double lambda_gcm2);

#endif /* OSH_NUCLEAR_TRIPATHI_H */
