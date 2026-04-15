#ifndef OSH_PHYSICS_STRAGGLING_H
#define OSH_PHYSICS_STRAGGLING_H

/**
 * @file osh_physics_straggling.h
 *
 * @brief Gaussian energy-straggling width (Bohr variance).
 *
 * @details
 * Energy straggling is the statistical spread in energy loss around the mean
 * predicted by the Bethe stopping-power formula.  The exact distribution is
 * the Vavilov distribution, which interpolates between the thick-absorber
 * (Gaussian) and thin-absorber (Landau) limits.
 *
 * This module implements the Gaussian (thick-absorber) limit using Bohr's
 * formula for the variance [Boh15, PDG34.2.2]:
 *
 *   σ²_Bohr = C_bohr × z_eff² × (Z/A) × d   [MeV²]
 *
 * where:
 *   C_bohr = K_bethe × m_e c²
 *           = 0.307075 × 0.511 = 0.1569 MeV² cm²/g
 *   z_eff  = effective charge of the projectile (from osh_physics_bethe_z_eff())
 *   Z/A    = effective Z/A of the target [mol/g]
 *   d      = areal density of the step ρ·ds [g/cm²]
 *
 * The Gaussian approximation is most accurate for thick absorbers (κ >> 1
 * in Vavilov's notation).  For the step sizes used in CSDA transport
 * (DELTAE ~ 2% of range), κ is typically >> 1 for protons in water above
 * ~10 MeV, making this a reasonable first-order model.  A future Vavilov
 * sampler can replace this function without changing the call sites.
 *
 * @par References
 * [Boh15] N. Bohr, Phil. Mag. 30, 581 (1915).  Straggling formula.
 * [PDG34.2.2] Zyla et al. (PDG), Prog. Theor. Exp. Phys. 2020, 083C01.
 */

/**
 * @brief Compute the Gaussian (Bohr) standard deviation of energy loss.
 *
 * @details
 * Returns σ = sqrt(C_bohr × z_eff² × z_over_a × thickness_gcm2) [MeV].
 *
 * The caller should then draw ΔE_strag ~ N(0, σ) and add it to the CSDA
 * mean energy loss for the step.  The result must be clamped to keep the
 * exit energy in [E_cutoff, E_entry].
 *
 * @param[in] z_eff           Effective projectile charge (dimensionless).
 * @param[in] z_over_a        Effective Z/A of the target [mol/g].
 * @param[in] thickness_gcm2  Areal density of the step ρ·ds [g/cm²].
 *
 * @returns σ [MeV], always ≥ 0.
 */
double osh_physics_straggling_sigma(double z_eff, double z_over_a, double thickness_gcm2);

#endif /* OSH_PHYSICS_STRAGGLING_H */
