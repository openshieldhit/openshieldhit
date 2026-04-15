#ifndef OSH_PHYSICS_MOLIERE_H
#define OSH_PHYSICS_MOLIERE_H

/**
 * @file osh_physics_moliere.h
 *
 * @brief Highland / Molière multiple Coulomb scattering for ions in matter.
 *
 * @details
 * Implements the Highland approximation to Molière theory: the projected
 * scattering-angle distribution is treated as Gaussian with RMS width θ₀
 * given by the Highland formula [Hig75, PDG34.3]:
 *
 *   θ₀ = (13.6 MeV / βcp) × z_eff × √(d/X₀) × [1 + 0.038 ln(z_eff² d/X₀)]
 *
 * where:
 *   βcp   = β × c × p  [MeV]  — velocity × momentum of the projectile
 *   z_eff = effective charge of the ion (see osh_physics_bethe_z_eff())
 *   d     = areal density of the step [g/cm²] = ρ × ds
 *   X₀    = radiation length of the medium [g/cm²]
 *
 * The Gaussian approximation is valid for d/X₀ in roughly [0.001, 0.1] and
 * overestimates the probability of large-angle single scatters; for the
 * clinical-beam use case (d/X₀ typically < 0.01 per step) this is adequate.
 *
 * The angular deflection is sampled by drawing two independent projected
 * angles θ_x, θ_y ~ N(0, θ₀) in the local transverse plane and rotating
 * the particle direction accordingly.
 *
 * @par References
 * [Hig75] Highland, NIMB 129, 497 (1975).  Original proposal.  \n
 * [Ferm54] Fermi, Notes on the theory of multiple scattering. \n
 * [PDG34.3] Zyla et al. (PDG), Prog. Theor. Exp. Phys. 2020, 083C01.
 */

struct osh_rng;

/**
 * @brief Compute the Highland RMS projected scattering angle θ₀.
 *
 * @details
 * Returns the Gaussian sigma of the projected angular deflection per step.
 * The log correction factor [1 + 0.038 ln(z_eff² d/X₀)] is included; if it
 * drives θ₀ negative (ultra-thin slab), θ₀ is clamped to zero.
 *
 * @param[in] t_total_mev   Total kinetic energy of the projectile [MeV].
 * @param[in] mass_mev      Nuclear rest mass of the projectile [MeV/c²].
 * @param[in] z_eff         Effective charge z_eff = Z_proj × GAMMA (dimensionless).
 * @param[in] thickness_gcm2 Areal density of the step ρ·ds [g/cm²].
 * @param[in] x0_gcm2       Radiation length of the medium [g/cm²]; must be > 0.
 *
 * @returns θ₀ [rad], always ≥ 0.
 */
double
osh_physics_moliere_theta0(double t_total_mev, double mass_mev, double z_eff, double thickness_gcm2, double x0_gcm2);

/**
 * @brief Maximum step length such that θ₀ does not exceed a given threshold.
 *
 * @details
 * Inverts the Highland formula (without the log correction) to find the areal
 * density d [g/cm²] at which the uncorrected RMS projected angle equals
 * @p theta_max_rad, then converts to a geometric step length [cm]:
 *
 *   s_theta = (X₀ / ρ) × (θ_max × βcp / (13.6 × z_eff))²
 *
 * Ignoring the log correction is conservative: the corrected θ₀ is slightly
 * larger than the uncorrected one for thick slabs, so the true θ₀ at s_theta
 * is marginally below @p theta_max_rad — an acceptable overestimate that avoids
 * an iterative solve.
 *
 * This is the third substep criterion alongside DELTAE (energy loss) and the
 * geometric boundary.  It keeps θ₀ inside the regime where the Gaussian/
 * Highland approximation is valid and prevents runaway large-angle deflections
 * in single substeps.
 *
 * Returns 0 (no limit applied) if any input is invalid: rho ≤ 0, z_eff ≤ 0,
 * x0_gcm2 ≤ 0, or theta_max_rad ≤ 0.  The caller should treat 0 as "no
 * constraint from this criterion."
 *
 * @param[in] t_total_mev   Total kinetic energy of the projectile [MeV].
 * @param[in] mass_mev      Nuclear rest mass of the projectile [MeV/c²].
 * @param[in] z_eff         Effective charge (Hubert GAMMA × Z_proj).
 * @param[in] rho_gcm3      Material density [g/cm³].
 * @param[in] x0_gcm2       Radiation length of the medium [g/cm²].
 * @param[in] theta_max_rad Target maximum θ₀ [rad].
 *
 * @returns Maximum step length [cm], or 0 if inputs are invalid.
 */
double osh_physics_moliere_s_theta(
    double t_total_mev, double mass_mev, double z_eff, double rho_gcm3, double x0_gcm2, double theta_max_rad);

/**
 * @brief Sample a scattered exit direction from a Gaussian MCS deflection.
 *
 * @details
 * Draws two independent projected-plane scattering angles
 *   θ_x, θ_y ~ N(0, θ₀)
 * in the transverse frame of the incident direction @p v, then constructs
 * the scattered direction @p w and renormalises it to unit length.
 *
 * For θ₀ < 1e-9 rad the function copies @p v into @p w unchanged.
 *
 * A stable orthonormal transverse basis is built by choosing the axis
 * least aligned with @p v and taking the cross products; this avoids the
 * singularity that appears if one just uses a fixed global axis.
 *
 * @param[in]  v      Incident unit direction vector (length 3).
 * @param[out] w      Exit unit direction vector (length 3).
 * @param[in]  theta0 RMS projected scattering angle θ₀ [rad].
 * @param[in]  rng    RNG state; consumes 2 Gaussian deviates per call.
 */
void osh_physics_moliere_scatter(double const v[3], double w[3], double theta0, struct osh_rng *rng);

#endif /* OSH_PHYSICS_MOLIERE_H */
