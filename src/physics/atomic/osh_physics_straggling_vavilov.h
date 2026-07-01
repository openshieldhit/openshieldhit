#ifndef OSH_PHYSICS_STRAGGLING_VAVILOV_H
#define OSH_PHYSICS_STRAGGLING_VAVILOV_H

/**
 * @file osh_physics_straggling_vavilov.h
 *
 * @brief Vavilov energy-straggling sampler (intermediate regime 0.01 ≤ κ < 10).
 *
 * @details
 * Returns the reduced Vavilov variable λ for a uniform deviate u, from our own
 * clean-room fit to the exact Vavilov distribution (Vavilov, Sov. Phys. JETP 5,
 * 749 (1957)).  The fit uses the vavinv *architecture* — region branches in
 * (u, κ); per region a polynomial in a transformed u (the pole-free Q=1 special
 * case of the rational form) whose coefficients are 2-D Chebyshev sums in
 * (ln κ, β²) — but none of Thomsen's / GEANT3's numbers.  Coefficients live in
 * the generated osh_physics_straggling_vavilov_coeffs.h.
 *
 * The caller forms the energy-loss fluctuation ΔE = ξ·(λ − λ̄) and adds it to the
 * CSDA mean loss (see osh_physics_straggling.h for ξ and λ̄).
 *
 * Pure function of (κ, β², u): no RNG, no state — the RNG draw stays in the
 * caller so a batched SoA form can be added without changing call sites.
 */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Reduced Vavilov variable λ(u; κ, β²).
 *
 * @param[in] kappa  Vavilov parameter κ ∈ [0.01, 10] (values outside are clamped
 *                   by the band dispatch to the nearest fitted band).
 * @param[in] beta2  Projectile β².
 * @param[in] u      Uniform deviate in (0, 1); clamped to the fitted range
 *                   [OSH_VAV_UMIN, 0.995].
 * @returns λ (dimensionless).
 */
double osh_physics_straggling_vavilov_lambda(double kappa, double beta2, double u);

#ifdef __cplusplus
}
#endif

#endif /* OSH_PHYSICS_STRAGGLING_VAVILOV_H */
