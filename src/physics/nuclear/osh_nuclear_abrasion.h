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
 * where σ_pN = OSH_ABRASION_SIGMA_PN_MB is a compile-time constant intended
 * for empirical tuning against public benchmark data.
 * The residual prefragment (A−ν, Z−ν_p) is handed to the Fermi break-up
 * stage with an excitation energy and a momentum from the event balance.
 *
 * This is deliberately a minimal development model: it gives the transport
 * layer realistic event topology (primary absorption plus fast nucleon
 * secondaries) while the full fragmentation/de-excitation machinery is being
 * built.  It is not intended to be a comprehensive nuclear reaction model.
 */

/** Inelastic p+nucleon cross section [mb] used for the mean participant count ⟨ν⟩.
 *  Tune against public secondary-nucleon yield benchmark data; override at compile time
 *  via -DOSH_ABRASION_SIGMA_PN_MB=<value>. */
#ifndef OSH_ABRASION_SIGMA_PN_MB
#define OSH_ABRASION_SIGMA_PN_MB 30.0
#endif

/** Mean prefragment excitation energy per abraded nucleon ("hole") [MeV],
 *  after Gaimard & Schmidt (1991).  The resulting E* is clamped to the
 *  leftover cascade-proton energy as an energy-conservation guard — that
 *  clamp is bookkeeping, not physics.  Override at compile time via
 *  -DOSH_ABRASION_EXCITATION_PER_HOLE_MEV=<value>. */
#ifndef OSH_ABRASION_EXCITATION_PER_HOLE_MEV
#define OSH_ABRASION_EXCITATION_PER_HOLE_MEV 13.3
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
 * scattering. The available proton energy is degraded after each sampled emission
 * so emitted kinetic energy does not exceed the incoming energy.
 *
 * Always sets event_out->kind = OSH_NUCLEAR_EVENT_ABRASION and primary_energy = 0
 * (primary absorbed).
 *
 * The surviving residual is written to event_out->fragments[0] with its mass
 * and atomic number reduced by the knocked-out nucleons, an excitation energy
 * E* = OSH_ABRASION_EXCITATION_PER_HOLE_MEV per abraded nucleon (clamped to
 * the leftover cascade-proton energy), and a lab momentum from the event
 * momentum balance — input for the Fermi break-up de-excitation stage.
 *
 * @param T_lab_mev      Incident proton kinetic energy [MeV].
 * @param incident_dir   Incident unit direction (length 3).
 * @param a_eff          Mass number of the struck target nucleus.
 * @param z_eff          Atomic number of the struck target nucleus.
 * @param sigma_pa_cm2   Tripathi σ_R for this target nucleus [cm²].
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
