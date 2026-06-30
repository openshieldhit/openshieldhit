#ifndef OSH_PHYSICS_SCAT_MOLIERE_H
#define OSH_PHYSICS_SCAT_MOLIERE_H

/**
 * @file osh_physics_scat_moliere.h
 *
 * @brief Full Bethe-Moliere multiple Coulomb scattering.
 *
 * @details
 * Clean-room implementation from published theory only; no source code or
 * tabulated numbers from Geant3 / SHIELD-HIT12A / third-party ports.
 *
 * References:
 * - G. Moliere, Z. Naturforsch. 2a, 133 (1947); 3a, 78 (1948).
 * - H. A. Bethe, Phys. Rev. 89, 1256 (1953).
 * - W. T. Scott, Rev. Mod. Phys. 35, 231 (1963).
 * - G. R. Lynch, O. I. Dahl, NIM B58, 6 (1991).
 */

struct osh_rng;

/**
 * @brief Precompute the Moliere reduced-angle distribution tables.
 *
 * @details
 * Fills the static f0/f1/f2 grids by numerically evaluating their Bethe
 * integral definitions.  Idempotent; call once before the parallel transport
 * loop.  It is also invoked lazily for single-threaded test/helper use.
 */
void osh_physics_moliere_init(void);

/**
 * @brief Solve Moliere's B equation B - ln B = ln(omega).
 *
 * @param[in]  omega  Effective number of scatters.
 * @param[out] b_out  Solved B (the larger root).
 * @returns 1 on success; 0 if no useful Moliere root exists.
 */
int osh_physics_moliere_solve_b(double omega, double *b_out);

/**
 * @brief Reduced-angle distribution function f_n(theta_reduced), n = 0,1,2.
 *
 * @param[in] n              0, 1, or 2.
 * @param[in] theta_reduced  Reduced angle >= 0.
 * @returns f_n(theta_reduced), or 0 outside the tabulated range / invalid n.
 */
double osh_physics_moliere_reduced_f(int n, double theta_reduced);

/**
 * @brief Inverse of the reduced-angle CDF: ϑ for given B and uniform deviate u.
 *
 * @details
 * O(1) bilinear interpolation of the precomputed inverse-CDF support table
 * (built once in osh_physics_moliere_init from our own f⁽ⁿ⁾ formulas).  The
 * deterministic (u) form is exposed so tests can probe monotonicity/bounds.
 *
 * @param[in] b  Moliere B parameter (clamped to the tabulated range).
 * @param[in] u  Uniform cumulative probability in [0,1].
 * @returns reduced angle ϑ ≥ 0.
 */
double osh_physics_moliere_inv_cdf(double b, double u);

/**
 * @brief Sample a reduced angle from the Moliere distribution for given B.
 *
 * @param[in] b    Moliere B parameter.
 * @param[in] rng  RNG state.
 * @returns reduced angle >= 0.
 */
double osh_physics_moliere_sample_reduced(double b, struct osh_rng *rng);

/**
 * @brief Sample and apply a full Bethe-Moliere angular deflection.
 *
 * @param[in]  v                Incident unit direction (length 3).
 * @param[out] w                Exit unit direction (length 3); always set.
 * @param[in]  t_kin_mev        Kinetic energy at mid-step [MeV].
 * @param[in]  mass_mev         Projectile rest mass [MeV/c^2].
 * @param[in]  z_eff            Effective projectile charge.
 * @param[in]  d_gcm2           Substep areal density rho*ds [g/cm^2].
 * @param[in]  path_scale_gcm2  Macroscopic path scale [g/cm^2].
 * @param[in]  chic2_coeff      Medium chi_c^2 coefficient.
 * @param[in]  screen_z         Medium effective screening Z.
 * @param[in]  rng              RNG state.
 * @returns 1 if a deflection was applied, 0 if @p w was copied unchanged.
 */
int osh_physics_moliere_scatter(double const v[3],
                                double w[3],
                                double t_kin_mev,
                                double mass_mev,
                                double z_eff,
                                double d_gcm2,
                                double path_scale_gcm2,
                                double chic2_coeff,
                                double screen_z,
                                struct osh_rng *rng);

#endif /* OSH_PHYSICS_SCAT_MOLIERE_H */
