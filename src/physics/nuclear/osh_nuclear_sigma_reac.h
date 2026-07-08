#ifndef OSH_NUCLEAR_SIGMA_REAC_H
#define OSH_NUCLEAR_SIGMA_REAC_H

/**
 * @file osh_nuclear_sigma_reac.h
 * @brief Proton-nucleus reaction (nonelastic) cross section with evaluated
 *        tables for selected targets (issue #277).
 *
 * The Tripathi parameterisation overpredicts the nonelastic cross section on
 * light therapy targets in the 50-150 MeV region (measured vs LA150/EXFOR:
 * p+O-16 up to +13% near 100 MeV, p+C-12 up to +24%).  For proton
 * projectiles on targets with a condensed evaluated table (O-16, C-12 —
 * ENDF/B-VIII.0 proton sublibrary, the LA150h lineage behind the ICRU 63
 * nonelastic tables), osh_nuclear_sigma_reac() interpolates the table;
 * everything else falls back to Tripathi unchanged.
 */

/**
 * @brief Reaction (nonelastic) cross section for a projectile on a target
 *        nucleus [cm^2].
 *
 * Drop-in replacement for osh_nuclear_tripathi_sigma(): evaluated-table
 * lookup (lin-lin in E, flat clamp above the table end) when @p zp / @p ap
 * is a proton and (@p zt, @p at) matches a tabulated target exactly;
 * Tripathi otherwise.  Returns 0 below the evaluation's threshold.
 *
 * @param zp                 Projectile charge number.
 * @param ap                 Projectile mass number.
 * @param zt                 Target charge number.
 * @param at                 Target mass number.
 * @param e_lab_per_nucleon  Lab kinetic energy per nucleon [MeV/u].
 * @return Cross section [cm^2]; 0 when the channel is closed.
 */
double osh_nuclear_sigma_reac(unsigned int zp, unsigned int ap, double zt, double at, double e_lab_per_nucleon);

#endif /* OSH_NUCLEAR_SIGMA_REAC_H */
