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
 * for empirical tuning against public benchmark data.  ν is sampled
 * zero-truncated (the reaction cross-section already fired, so at least one
 * nucleon participates).
 *
 * Intranuclear-cascade picture: the proton undergoes ν quasi-elastic
 * collisions, knocking out one nucleon each and being deflected collision by
 * collision; the degraded cascade proton then escapes as a transportable
 * secondary (or is absorbed below ~1 MeV).  Each hole charges
 * OSH_ABRASION_EXCITATION_PER_HOLE_MEV to the cascade proton, booked as
 * prefragment excitation.  Low-energy knockouts are re-absorbed into the
 * prefragment (OSH_ABRASION_RETENTION_THRESHOLD_MEV, issue #263): their
 * kinetic energy funds E* and they become particle excitons.  Kinetic
 * energy is conserved exactly:
 *
 *   T_in = Σ KE_escaped + E*
 *
 * The prefragment exports its exciton configuration (particles, holes) in
 * struct osh_nuclear_fragment for the pre-equilibrium stage (issue #221).
 *
 * The residual prefragment (A−ν, Z−ν_p) is handed to the Fermi break-up
 * stage with this excitation energy and a momentum from the event balance.
 *
 * This is deliberately a minimal development model: it gives the transport
 * layer realistic event topology (primary termination plus fast nucleon
 * secondaries) while the full fragmentation/de-excitation machinery is being
 * built.  It is not intended to be a comprehensive nuclear reaction model.
 */

/** Inelastic p+nucleon cross section [mb] used for the mean participant count ⟨ν⟩.
 *  Tune against public secondary-nucleon yield benchmark data; override at compile time
 *  via -DOSH_ABRASION_SIGMA_PN_MB=<value>. */
#ifndef OSH_ABRASION_SIGMA_PN_MB
#define OSH_ABRASION_SIGMA_PN_MB 30.0
#endif

/** Prefragment excitation energy per abraded nucleon ("hole") [MeV], after
 *  Gaimard & Schmidt (1991).  Charged to the continuing cascade proton at
 *  each collision and booked as prefragment excitation, so the event energy
 *  balance is exact; if the proton cannot pay it is absorbed and its
 *  remaining energy funds E*.  Override at compile time via
 *  -DOSH_ABRASION_EXCITATION_PER_HOLE_MEV=<value>.
 *  The #260 knob sweep showed this constant must NOT be retuned: raising it
 *  fixes the alpha mean energy only by breaking alpha fluence and dose. */
#ifndef OSH_ABRASION_EXCITATION_PER_HOLE_MEV
#define OSH_ABRASION_EXCITATION_PER_HOLE_MEV 13.3
#endif

/** Knockout-nucleon retention threshold [MeV] (issue #263): a knocked-out
 *  nucleon with lab kinetic energy in the neighbourhood of this value or
 *  below is re-absorbed into the prefragment instead of escaping — the lean
 *  analog of INCL4.6's "back to spectator" recipe (Boudard et al., PRC 87
 *  (2013) 014606, Sec. II D1: neutrons return to spectator status below the
 *  neutron-emission threshold, xi ~= 7 MeV; protons below their emission
 *  threshold plus ~2/3 of the Coulomb barrier).  The same paper's earlier
 *  xi = 18 MeV variant improved low-energy cluster yields but distorted the
 *  nucleon spectra — do not raise this knob to chase E*.  The retained
 *  kinetic energy is booked into E* and the nucleon becomes a particle
 *  exciton, so the event energy balance stays exact.  Override via
 *  -DOSH_ABRASION_RETENTION_THRESHOLD_MEV=<value>; 0 disables retention. */
#ifndef OSH_ABRASION_RETENTION_THRESHOLD_MEV
#define OSH_ABRASION_RETENTION_THRESHOLD_MEV 7.0
#endif

/** Width [MeV] of the smooth (Fermi-function) turn-on of the retention
 *  probability around the threshold; a hard step would imprint a spurious
 *  edge on the escaping-nucleon spectrum.  Override via
 *  -DOSH_ABRASION_RETENTION_WIDTH_MEV=<value>. */
#ifndef OSH_ABRASION_RETENTION_WIDTH_MEV
#define OSH_ABRASION_RETENTION_WIDTH_MEV 2.0
#endif

/** Nuclear radius parameter [fm] for the retention Coulomb-barrier estimate
 *  V_C = e² Z / (r0 A^(1/3)); protons add ~2/3 of V_C to the retention
 *  threshold (INCL4.6 recipe). */
#define OSH_ABRASION_COULOMB_R0_FM 1.25

#ifdef __cplusplus
extern "C" {
#endif

struct osh_rng;
struct osh_nuclear_event;

/**
 * @brief Sample fast-nucleon emission for one nuclear abrasion event.
 *
 * @details
 * Samples ν ~ Poisson(⟨ν⟩ | ν ≥ 1), then for each participant emits a neutron
 * (prob N/A) or proton (prob Z/A) using quasi-elastic equal-mass kinematics
 * with isotropic CM scattering, deflecting the continuing cascade proton at
 * each collision.  A knockout near or below the retention threshold is
 * re-absorbed instead of emitted: its kinetic energy funds E*, it stays in
 * the residue (no A/Z decrement), and it counts as a particle exciton.
 * After the cascade, the degraded proton escapes as an additional secondary
 * (species proton, generation +1); below ~1 MeV it is absorbed into the
 * prefragment instead (also counted as a particle exciton).
 *
 * Always sets event_out->kind = OSH_NUCLEAR_EVENT_ABRASION and primary_energy = 0
 * (the primary slot is terminated; the cascade proton continues as a secondary).
 *
 * The surviving residual is written to event_out->fragments[0] with its mass
 * and atomic number reduced by the escaped knockouts — and increased by a
 * captured cascade proton (issue #265; skipped when the nuclide would leave
 * the de-excitation table domain, e.g. full capture on O-16) — the
 * accumulated excitation energy E* (hole charges + retained kinetic energy),
 * the exciton configuration (particles, holes), and a lab momentum from the
 * event momentum balance — input for the pre-equilibrium / Fermi break-up
 * de-excitation stages.  Kinetic energy is conserved exactly:
 * T_in = Σ KE_secondaries + E*.
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
