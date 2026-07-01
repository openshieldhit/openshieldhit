#ifndef OSH_PHYSICS_STRAG_LANDAU_H
#define OSH_PHYSICS_STRAG_LANDAU_H

/**
 * @file osh_physics_strag_landau.h
 *
 * @brief Landau energy-straggling sampler (thin-absorber regime κ < 0.01).
 *
 * @details
 * Returns the universal reduced Landau variable λ for a uniform deviate u, from
 * a reproducible fit of the standard (DENLAN/CERNLIB) Landau inverse CDF — the
 * κ→0 limit of the Vavilov distribution.  The Landau λ is universal (no κ or β²
 * dependence; those enter the energy-loss fluctuation only through ξ and λ̄), so
 * this is a 1-D piecewise polynomial in a transformed u, following an approach
 * by Bjarne Thomsen (Dept. of Physics and Astronomy, Aarhus University, Aarhus,
 * Denmark, 2012).  The coefficients (generated into
 * osh_physics_strag_landau_coeffs.h by tools/vavilov_fit) are
 * OpenShieldHIT-specific, regenerated from the distribution.
 *
 * The caller forms ΔE = ξ·(λ − λ̄) and adds it to the CSDA mean loss.
 *
 * Pure function of u: no RNG, no state (RNG draw stays in the caller).
 *
 * @par References
 * Landau, J. Phys. USSR 8, 201 (1944).  Kölbig & Schorr, Comput. Phys. Commun.
 * 31, 97 (1984) (the DENLAN convention this matches).
 */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Universal Landau variable λ(u).
 *
 * @param[in] u  Uniform deviate in (0, 1); clamped to the fitted range
 *               [OSH_LAN_UMIN, OSH_LAN_UMAX].
 * @returns λ (dimensionless).
 */
double osh_physics_strag_landau_lambda(double u);

#ifdef __cplusplus
}
#endif

#endif /* OSH_PHYSICS_STRAG_LANDAU_H */
