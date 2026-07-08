#ifndef OSH_PHYSICS_STRAG_H
#define OSH_PHYSICS_STRAG_H

/**
 * @file osh_physics_strag.h
 *
 * @brief Energy-straggling dispatcher and shared straggling kinematics.
 *
 * @details
 * Top-level module for the statistical spread in energy loss around the CSDA
 * mean.  It owns the shared kinematics (ξ, the maximum energy transfer E_max,
 * the Vavilov parameter κ = ξ/E_max, and the reduced-mean λ̄) and dispatches to
 * the per-model samplers, mirroring the multiple-scattering module split
 * (osh_physics_scat.h → highland/moliere):
 *
 *   - Gaussian (Bohr), κ ≳ 10      → osh_physics_strag_gauss.h
 *   - Vavilov, 0.01 ≤ κ < 10       → osh_physics_strag_vavilov.h
 *   - Landau, κ < 0.01             → osh_physics_strag_landau.h
 *
 * The per-model λ-samplers are pure functions of (κ, β², u) with no RNG or
 * mutable state, so the RNG draw lives in the caller/dispatcher and a batched
 * SoA form can be added later without touching call sites.
 *
 * @par References
 * PDG "Passage of Particles through Matter", Prog. Theor. Exp. Phys. 2020.
 * Vavilov, Sov. Phys. JETP 5, 749 (1957).  Landau, J. Phys. USSR 8, 201 (1944).
 */

#include "physics/atomic/osh_physics_strag_gauss.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Regime boundaries in the Vavilov parameter κ = ξ/E_max (STRAGG 2 auto-dispatch):
 * κ ≥ GAUSS → thick absorber, Bohr Gaussian; κ < LANDAU → thin absorber, Landau;
 * in between, Vavilov.  Matches the GEANT3 GFLUCT thresholds. */
#define OSH_STRAG_KAPPA_GAUSS 10.0
#define OSH_STRAG_KAPPA_LANDAU 0.01

/** Energy-loss straggling model; integer values mirror enum
 *  osh_transport_straggling_mode. */
enum osh_strag_model {
    OSH_STRAG_OFF = 0,
    OSH_STRAG_GAUSSIAN = 1, /**< Bohr Gaussian (thick absorber). */
    OSH_STRAG_VAVILOV = 2,  /**< Vavilov; auto Gaussian/Landau by κ. */
    OSH_STRAG_URBAN = 3     /**< Reserved (G4 Urbán); not implemented. */
};

/**
 * @brief Landau/Vavilov width parameter ξ for the step [MeV].
 *
 * @details ξ = (K_bethe/2) · (Z/A) · (z_eff²/β²) · d, with d = ρ·ds the areal
 * density [g/cm²] and Z/A in [mol/g].
 *
 * @param[in] z_eff           Effective projectile charge (dimensionless).
 * @param[in] z_over_a        Effective Z/A of the target [mol/g].
 * @param[in] thickness_gcm2  Areal density of the step ρ·ds [g/cm²].
 * @param[in] beta2           Projectile β² (v²/c²), in (0, 1).
 * @returns ξ [MeV], or 0 for non-physical inputs.
 */
double osh_physics_strag_xi(double z_eff, double z_over_a, double thickness_gcm2, double beta2);

/**
 * @brief Maximum energy transfer to a target electron E_max (W_max) [MeV].
 *
 * @details E_max = 2·mₑc²·β²γ² / (1 + 2γ·mₑ/M + (mₑ/M)²), with γ = 1 + T/M and
 * β² = 1 − 1/γ².
 *
 * @param[in] t_kin_mev  Projectile kinetic energy T [MeV].
 * @param[in] mass_mev   Projectile rest mass M [MeV/c²].
 * @returns E_max [MeV], or 0 for non-physical inputs.
 */
double osh_physics_strag_emax(double t_kin_mev, double mass_mev);

/**
 * @brief Maximum energy transfer when projectile β² is already known.
 *
 * @details Same formula and units as @ref osh_physics_strag_emax, but avoids
 * recomputing β² in callers that already need it for ξ and the samplers.
 *
 * @param[in] t_kin_mev  Projectile kinetic energy T [MeV].
 * @param[in] mass_mev   Projectile rest mass M [MeV/c²].
 * @param[in] beta2      Projectile β² (v²/c²), in (0, 1).
 * @returns E_max [MeV], or 0 for non-physical inputs.
 */
double osh_physics_strag_emax_beta2(double t_kin_mev, double mass_mev, double beta2);

/**
 * @brief Vavilov straggling parameter κ = ξ / E_max.
 *
 * @param[in] xi     Width parameter ξ [MeV].
 * @param[in] e_max  Maximum energy transfer E_max [MeV].
 * @returns κ (dimensionless), or 0 when @p e_max ≤ 0.
 */
double osh_physics_strag_kappa(double xi, double e_max);

/**
 * @brief Reduced-mean λ̄ of the Vavilov/Landau variable for this step.
 *
 * @details λ̄ = −β² − ln κ + (γ_Euler − 1).  The sampled fluctuation of the
 * energy loss is ΔE = ξ·(λ − λ̄), which is mean-preserving by construction.
 *
 * @param[in] kappa  Vavilov parameter κ (> 0).
 * @param[in] beta2  Projectile β².
 * @returns λ̄ (dimensionless).
 */
double osh_physics_strag_lambda_bar(double kappa, double beta2);

/**
 * @brief Reduced-mean λ̄ when ln κ is already known.
 *
 * @details Same formula as @ref osh_physics_strag_lambda_bar, but lets hot-path
 * callers share ln κ with the Vavilov sampler.
 *
 * @param[in] log_kappa  Natural logarithm of κ.
 * @param[in] beta2      Projectile β².
 * @returns λ̄ (dimensionless), or NAN for a non-finite @p log_kappa.
 */
double osh_physics_strag_lambda_bar_log(double log_kappa, double beta2);

#ifdef __cplusplus
}
#endif

#endif /* OSH_PHYSICS_STRAG_H */
