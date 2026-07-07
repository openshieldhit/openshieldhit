#ifndef OSH_NUCLEAR_PREEQ_H
#define OSH_NUCLEAR_PREEQ_H

/**
 * @file osh_nuclear_preeq.h
 * @brief Pre-equilibrium emission — single-component exciton model (issue #225).
 *
 * @details
 * De-excites the abrasion prefragment BEFORE the equilibrium Fermi break-up
 * stage: starting from the exciton configuration (particles p, holes h) and
 * excitation energy E* exported by the fast stage, the residue either emits a
 * fast ejectile {n, p, d, t, He-3, He-4} or thermalizes one step (a Delta-n
 * = +2 intranuclear collision), until it reaches the equilibrium exciton
 * number n_eq = sqrt(2 g E*), runs out of excitation, or exhausts its
 * particle excitons.  This is the mechanism that produces the fast light
 * clusters (in particular the hard alpha/deuteron components) that pure
 * statistical break-up cannot generate — see the #260 measurement and the
 * roadmap on #221.
 *
 * Model (approximate by design; the #221 tolerance note applies):
 *  - Exciton state density: Williams,
 *      omega(p, h, E) = g (gE)^(p+h-1) / (p! h! (p+h-1)!),
 *    with a single-particle level density g = OSH_PREEQ_LEVEL_DENSITY_PER_A
 *    * A (equidistant model; no pairing or finite-well corrections).
 *  - Emission width, detailed balance (Betak-Dobes form; equals the
 *    Geant4/CEM closed forms for nucleons when the same g is used):
 *      lambda_j(T) dT = (2 s_j + 1) mu_j T sigma_inv(T) / (pi^2 (hbar c)^3)
 *                       * gamma_j R_j [omega(p - a_j, h, U) / omega(p, h, E)] dT,
 *    U = E* - B_j - T, with B_j the exact separation energy from the isotope
 *    mass table.  R_j is a composition factor C(a_j, z_j) f_p^z_j f_n^(a_j-z_j)
 *    built from the residue proton/neutron fractions (single-component model;
 *    the exciton charge split is not tracked).
 *  - Cluster condensation probability gamma_j (Geant4 closed forms):
 *      gamma_d = 16/A, gamma_t = gamma_He3 = 243/A^2, gamma_alpha = 4096/A^3,
 *    each multiplied by one calibration scale knob (CEM fits the equivalent
 *    factor against data — expect to tune, see OSH_PREEQ_GAMMA_SCALE_*).
 *  - Inverse cross sections: Dostrovsky-style geometric form,
 *    sigma = pi R^2 alpha_n (1 + beta_n / T) for neutrons and
 *    sigma = pi R^2 (1 - V_C / T) (clamped at 0) for charged ejectiles,
 *    R = OSH_PREEQ_EMISSION_R0_FM (A_res^(1/3) + a_j^(1/3)).
 *  - Transition rate lambda_plus (Delta-n = +2 only, "never go back"):
 *    CEM form — averaged in-medium NN cross section at the mean exciton
 *    relative energy 1.6 E_F + E_star_per_exciton, Pauli-blocking factor,
 *    interaction volume from OSH_PREEQ_TRANSITIONS_R0_FM
 *    (Gudima-Mashnik-Toneev; the Geant4 G4PreCompoundTransitions CEM branch
 *    implements the same expressions).
 *  - Emission angles: isotropic (lab).  Kalbach continuum systematics are a
 *    planned refinement within #225 calibration.
 *
 * Contracts:
 *  - A fragment with (p, h) = (0, 0) is thermalized: the step is a no-op
 *    (this preserves the neutron-capture compound path and all break-up
 *    residues untouched).
 *  - Emission bookkeeping is exact: E* before equals E* after plus B_j
 *    plus T for every emission (mass-table B_j); fragment momentum is
 *    reduced by every emitted ejectile; A, Z and exciton counts are updated
 *    in place.
 *  - Only whitelist species {n, p, d, t, He-3, He-4} are emitted, appended
 *    to the caller-owned event; the step performs no allocation.
 */

#include <stddef.h>

#include "openshieldhit/status.h"
#include "particle/osh_particle.h"

/** Single-particle level density per nucleon [1/MeV]: g = value * A.
 *  Geant4 de-excitation default (G4DeexPrecoParameters). */
#ifndef OSH_PREEQ_LEVEL_DENSITY_PER_A
#define OSH_PREEQ_LEVEL_DENSITY_PER_A 0.075
#endif

/** Fermi energy [MeV] entering the CEM transition rate. */
#ifndef OSH_PREEQ_FERMI_ENERGY_MEV
#define OSH_PREEQ_FERMI_ENERGY_MEV 35.0
#endif

/** Radius parameter [fm] of the intranuclear interaction volume in the CEM
 *  transition rate (Geant4 default 0.6 fm). */
#ifndef OSH_PREEQ_TRANSITIONS_R0_FM
#define OSH_PREEQ_TRANSITIONS_R0_FM 0.6
#endif

/** Radius parameter [fm] of the emission inverse-cross-section estimate. */
#ifndef OSH_PREEQ_EMISSION_R0_FM
#define OSH_PREEQ_EMISSION_R0_FM 1.25
#endif

/** Per-species calibration scales on the cluster condensation probability
 *  gamma_j.  CEM03.03 multiplies its own gamma_j ("a rather crude estimate")
 *  by empirically fitted factors; these knobs are that fit surface.  1.0 =
 *  bare Geant4 closed forms. */
#ifndef OSH_PREEQ_GAMMA_SCALE_D
#define OSH_PREEQ_GAMMA_SCALE_D 1.0
#endif
#ifndef OSH_PREEQ_GAMMA_SCALE_T
#define OSH_PREEQ_GAMMA_SCALE_T 1.0
#endif
#ifndef OSH_PREEQ_GAMMA_SCALE_HE3
#define OSH_PREEQ_GAMMA_SCALE_HE3 1.0
#endif
#ifndef OSH_PREEQ_GAMMA_SCALE_ALPHA
#define OSH_PREEQ_GAMMA_SCALE_ALPHA 1.0
#endif

/** Dense mass-table domain (shared with the Fermi break-up stage). */
#define OSH_PREEQ_ZMAX 8
#define OSH_PREEQ_AMAX 16

/** Number of emittable ejectile species: n, p, d, t, He-3, He-4. */
#define OSH_PREEQ_NSPECIES 6

#ifdef __cplusplus
extern "C" {
#endif

struct osh_rng;
struct osh_nuclear_event;
struct osh_nuclear_fragment;

/** Compiled pre-equilibrium model: dense ground-state nuclear mass table
 *  [MeV/c^2] over z in [0, OSH_PREEQ_ZMAX], a in [0, OSH_PREEQ_AMAX], plus
 *  the emitted-species descriptors.  No heap allocation; compile fills the
 *  tables from the isotope database. */
struct osh_nuclear_preeq {
    double mass_mev[(OSH_PREEQ_ZMAX + 1) * (OSH_PREEQ_AMAX + 1)];
    struct particle species[OSH_PREEQ_NSPECIES];
    int compiled;
};

/**
 * @brief Build the mass tables from the isotope database.
 *
 * @param model  Model storage (caller-owned).
 * @returns OSH_OK on success.
 */
enum osh_status osh_nuclear_preeq_compile(struct osh_nuclear_preeq *model);

/**
 * @brief Run pre-equilibrium emission on one excited fragment, in place.
 *
 * @details
 * Emits fast ejectiles into event_out->secondaries and cools the fragment
 * (E*, A, Z, momentum, exciton counts updated) until equilibrium.  A
 * fragment outside the mass-table domain, without particle excitons, or
 * below every emission threshold is left untouched.  Sets
 * event_out->kind = OSH_NUCLEAR_EVENT_FRAGMENTATION when at least one
 * ejectile was emitted.  No allocation; only the RNG state and the given
 * structs are mutated.
 *
 * @param model     Compiled model.
 * @param fragment  Excited fragment (mutated in place).
 * @param rng       RNG state (mutated).
 * @param event_out Event to append emitted secondaries to.
 */
void osh_nuclear_preeq_step(struct osh_nuclear_preeq const *model,
                            struct osh_nuclear_fragment *fragment,
                            struct osh_rng *rng,
                            struct osh_nuclear_event *event_out);

#ifdef __cplusplus
}
#endif

#endif /* OSH_NUCLEAR_PREEQ_H */
