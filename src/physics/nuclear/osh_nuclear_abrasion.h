#ifndef OSH_NUCLEAR_ABRASION_H
#define OSH_NUCLEAR_ABRASION_H

/**
 * @file osh_nuclear_abrasion.h
 * @brief Nuclear abrasion step — fast-nucleon emission in the Abrasion-Ablation model.
 *
 * @details
 * Implements the abrasion (fast) stage of the Bowman-Swiatecki-Tsang (1973)
 * Abrasion-Ablation model.  Uses the wounded-nucleon picture (Glauber 1970):
 *
 *   ⟨ν⟩ = σ_pN · A / σ_pA
 *
 * where σ_pN = OSH_ABRASION_SIGMA_PN_MB is a compile-time constant tuned
 * empirically against secondary proton yield benchmarks (FLUKA / SH12A).
 * The residual nucleus (A−ν, Z−ν_p) is discarded; FermiBreakup de-excitation
 * is a planned follow-on.
 */

/** Inelastic p+nucleon cross section [mb] used for the mean participant count ⟨ν⟩.
 *  Tune against FLUKA / SH12A secondary proton yield data; override at compile time
 *  via -DOSH_ABRASION_SIGMA_PN_MB=<value>. */
#ifndef OSH_ABRASION_SIGMA_PN_MB
#define OSH_ABRASION_SIGMA_PN_MB 30.0
#endif

#ifdef __cplusplus
extern "C" {
#endif

struct osh_rng;
struct osh_nuclear_event;

/**
 * @brief Sample fast-nucleon emission for one nuclear abrasion event.
 *
 * @details
 * Samples ν ~ Poisson(⟨ν⟩), then for each participant emits a neutron (prob N/A)
 * or proton (prob Z/A) using quasi-elastic equal-mass kinematics with isotropic CM
 * scattering.  Sequential collisions degrade the proton energy for conservation.
 *
 * Always sets event_out->kind = OSH_NUCLEAR_EVENT_ABRASION and primary_energy = 0
 * (primary absorbed).
 *
 * @param T_lab_mev      Incident proton kinetic energy [MeV].
 * @param incident_dir   Incident unit direction (length 3).
 * @param a_eff          Number-density-weighted mean A of the material.
 * @param z_eff          Number-density-weighted mean Z of the material.
 * @param sigma_pa_cm2   Tripathi σ_R for this material [cm²].
 * @param rng            RNG state (mutated).
 * @param event_out      Written with ABRASION event and emitted secondaries.
 */
void osh_nuclear_abrasion_step(double T_lab_mev,
                               double const incident_dir[3],
                               double a_eff,
                               double z_eff,
                               double sigma_pa_cm2,
                               struct osh_rng *rng,
                               struct osh_nuclear_event *event_out);

#ifdef __cplusplus
}
#endif

#endif /* OSH_NUCLEAR_ABRASION_H */
