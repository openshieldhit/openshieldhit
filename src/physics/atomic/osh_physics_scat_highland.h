#ifndef OSH_PHYSICS_SCAT_HIGHLAND_H
#define OSH_PHYSICS_SCAT_HIGHLAND_H

/**
 * @file osh_physics_scat_highland.h
 *
 * @brief Highland Gaussian multiple Coulomb scattering for ions in matter.
 *
 * @details
 * Implements the Highland approximation to Moliere theory: the projected
 * scattering-angle distribution is treated as Gaussian with RMS width theta0
 * given by the Highland formula [Hig75, PDG34.3]:
 *
 *   theta0 = (13.6 MeV / beta cp) * z_eff * sqrt(d/X0)
 *            * [1 + 0.038 ln(z_eff^2 s/X0)]
 *   (with s = path_scale_gcm2; s <= 0 falls back to d)
 *
 * The implementation can evaluate the logarithmic correction at a larger
 * macroscopic path scale while keeping the leading variance on the substep
 * thickness.  This lets independent substep variances add back to the same
 * full-path Highland value.
 */

struct osh_rng;

/**
 * @brief Compute the Highland RMS projected scattering angle theta0.
 *
 * @param[in] t_total_mev     Total kinetic energy of the projectile [MeV].
 * @param[in] mass_mev        Nuclear rest mass of the projectile [MeV/c^2].
 * @param[in] z_eff           Effective charge z_eff = Z_proj * GAMMA.
 * @param[in] thickness_gcm2  Areal density of the step rho*ds [g/cm^2].
 * @param[in] path_scale_gcm2 Macroscopic path scale for the log correction [g/cm^2];
 *                            <= 0 falls back to @p thickness_gcm2.
 * @param[in] x0_gcm2         Radiation length of the medium [g/cm^2].
 *
 * @returns theta0 [rad], always >= 0.
 */
double osh_physics_highland_theta0(
    double t_total_mev, double mass_mev, double z_eff, double thickness_gcm2, double path_scale_gcm2, double x0_gcm2);

/**
 * @brief Maximum step length such that theta0 does not exceed a threshold.
 *
 * @param[in] t_total_mev   Total kinetic energy of the projectile [MeV].
 * @param[in] mass_mev      Nuclear rest mass of the projectile [MeV/c^2].
 * @param[in] z_eff         Effective charge.
 * @param[in] rho_gcm3      Material density [g/cm^3].
 * @param[in] x0_gcm2       Radiation length of the medium [g/cm^2].
 * @param[in] theta_max_rad Target maximum theta0 [rad].
 *
 * @returns Maximum step length [cm], or 0 if inputs are invalid.
 */
double osh_physics_highland_s_theta(
    double t_total_mev, double mass_mev, double z_eff, double rho_gcm3, double x0_gcm2, double theta_max_rad);

/**
 * @brief Maximum step length such that the lateral displacement built up over
 *        one step stays below a tolerance.
 *
 * @details
 * The stepper models one substep as a straight leg plus a single deflection
 * (at a random hinge, or at the step end for boundary-limited steps), so the
 * modelled trajectory can miss the true random walk by up to theta0(s)*s in
 * the transverse direction anywhere inside the step.  Capping the step length
 * therefore caps that error, which is what makes scoring surfaces placed
 * strictly inside a homogeneous zone see multiple scattering at all (issue
 * #325).
 *
 * With the log correction evaluated at the macroscopic path scale, theta0
 * grows as sqrt(s), so theta0(s)*s = A*s^(3/2) with A = theta0(1 cm) and the
 * cap inverts in closed form:  s_max = (lateral_max / A)^(2/3).
 *
 * @param[in] t_total_mev     Total kinetic energy of the projectile [MeV].
 * @param[in] mass_mev        Nuclear rest mass of the projectile [MeV/c^2].
 * @param[in] z_eff           Effective charge.
 * @param[in] rho_gcm3        Material density [g/cm^3].
 * @param[in] path_scale_gcm2 Macroscopic path scale for the log correction [g/cm^2];
 *                            <= 0 falls back to the 1 cm reference thickness.
 * @param[in] x0_gcm2         Radiation length of the medium [g/cm^2].
 * @param[in] lateral_max_cm  Tolerated per-step lateral displacement [cm].
 *
 * @returns Maximum step length [cm], or 0 if inputs are invalid or the medium
 *          imposes no meaningful limit (treat 0 as "no limit").
 */
double osh_physics_highland_s_lateral(double t_total_mev,
                                      double mass_mev,
                                      double z_eff,
                                      double rho_gcm3,
                                      double path_scale_gcm2,
                                      double x0_gcm2,
                                      double lateral_max_cm);

/**
 * @brief Sample a scattered exit direction from a Gaussian MCS deflection.
 *
 * @param[in]  v      Incident unit direction vector (length 3).
 * @param[out] w      Exit unit direction vector (length 3).
 * @param[in]  theta0 RMS projected scattering angle [rad].
 * @param[in]  rng    RNG state; consumes 2 Gaussian deviates per call.
 */
void osh_physics_highland_scatter(double const v[3], double w[3], double theta0, struct osh_rng *rng);

#endif /* OSH_PHYSICS_SCAT_HIGHLAND_H */
