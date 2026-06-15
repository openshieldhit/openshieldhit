#ifndef OSH_NUCLEAR_FERMI_BREAKUP_H
#define OSH_NUCLEAR_FERMI_BREAKUP_H

/**
 * @file osh_nuclear_fermi_breakup.h
 * @brief Fermi break-up de-excitation of light excited prefragments.
 *
 * @details
 * De-excites the residual prefragment (Z, A, E*) produced by the abrasion
 * stage into nucleons and light ions (Fermi 1950 statistical break-up,
 * applicable to light nuclei A <= 16).
 *
 * This is a semiphysical placeholder for the full Fermi break-up model:
 * instead of sampling simultaneous n-body partitions of the parent, the
 * de-excitation is realised as sequential binary splits, each weighted by
 * the two-body phase-space factor
 *
 *   w = g1 * g2 * mu^(3/2) * sqrt(E_kin),   E_kin = E* + Q
 *
 * where g = 2J+1 are ground-state spin degeneracies and mu is the reduced
 * mass.  Repeated two-body decay can bias multiplicities and kinetic-energy
 * spectra relative to the canonical simultaneous break-up; this is accepted
 * for the placeholder model while the full n-body Fermi break-up / SMM work
 * remains pending.  Particle-unstable nuclides present in the isotope
 * database (Be-8, He-5, Li-5, ...) participate as intermediate states and
 * decay further through their own open channels (e.g. Be-8 -> 2 alpha).
 *
 * Final-product policy (strict): only the whitelist n, p, d, t, He-3, He-4
 * is ever emitted as transportable secondaries.  Any other nuclide without
 * an open channel (particle-stable residues such as Li-7 or Be-9, or a
 * numerically closed channel) is returned as an *unprocessed fragment* in
 * the event so the transport layer counts it in the fragment pool; it is
 * never pushed onto a particle stack.
 *
 * The channel table and product species descriptors are compiled once at
 * startup from the isotope database; the per-event step performs no
 * allocations and only mutates the RNG state, matching the pool-independent
 * wavefront design of the nuclear handler.
 */

#include <stddef.h>
#include <stdint.h>

#include "openshieldhit/status.h"

#ifdef __cplusplus
extern "C" {
#endif

struct osh_nuclear_event;
struct osh_nuclear_fragment;
struct osh_rng;
struct particle;

/** Maximum parent/product atomic number handled by the break-up table. */
#define OSH_FERMI_BREAKUP_ZMAX 8

/** Maximum parent/product mass number handled by the break-up table. */
#define OSH_FERMI_BREAKUP_AMAX 16

/** One binary decay channel parent -> (z1,a1) + (z2,a2). */
struct osh_fermi_channel {
    float q_mev;            /**< Q = M_parent - m1 - m2 (ground states) [MeV]. */
    float weight_prefactor; /**< g1*g2*mu^(3/2); halved for identical products. */
    uint8_t z1;             /**< First product atomic number.                   */
    uint8_t a1;             /**< First product mass number.                     */
    uint8_t z2;             /**< Second product atomic number.                  */
    uint8_t a2;             /**< Second product mass number.                    */
};

/**
 * @brief Compiled Fermi break-up model: channel table and product species.
 *
 * @details
 * All arrays are indexed by the dense (z, a) key z*(AMAX+1)+a for
 * z in [0, ZMAX], a in [0, AMAX].  Compiled once by
 * osh_nuclear_fermi_breakup_compile(); owned by the nuclear handler whose
 * lifetime must exceed the transport run (pool entries borrow species
 * pointers from this model).  Free with osh_nuclear_fermi_breakup_free().
 */
struct osh_nuclear_fermi_breakup {
    struct osh_fermi_channel *channel_pool; /**< Flat array of all parents' channels. */
    double *mass_mev;                       /**< Ground-state nuclear mass per dense (z,a); 0 if absent. */
    struct particle *species;               /**< Final-product descriptors (n, p, d, t, He-3, He-4).    */
    uint16_t *chan_offset;                  /**< chan_offset[idx]: start of parent idx in channel_pool.  */
    uint16_t *chan_count;                   /**< chan_count[idx]:  channel count of parent idx.          */
    size_t nchannels;                       /**< Total entries in channel_pool.                          */
};

/**
 * @brief Compile the break-up channel table from the isotope database.
 *
 * @details
 * Enumerates, for every parent (z = 1..ZMAX, a = 2..AMAX) present in the
 * isotope database, all binary partitions into two database nuclides and
 * stores their Q-value and phase-space prefactor.  Also builds the
 * final-product species descriptors via osh_particle_from_pdg().
 *
 * Must be called once at setup, never on the hot path.
 *
 * @param[out] out  Model to populate (must point to zero-initialised storage).
 * @returns OSH_OK on success, OSH_ENOMEM on allocation failure,
 *          OSH_ESTATE if the particle registry lacks a required species.
 */
enum osh_status osh_nuclear_fermi_breakup_compile(struct osh_nuclear_fermi_breakup *out);

/**
 * @brief Free all memory owned by the model.
 *
 * @param m  Model to free.  May be NULL.  Fields are zeroed after free.
 */
void osh_nuclear_fermi_breakup_free(struct osh_nuclear_fermi_breakup *m);

/**
 * @brief De-excite one prefragment, appending break-up products to the event.
 *
 * @details
 * If the fragment is outside the model domain (A < 2, A > AMAX or Z > ZMAX)
 * or has no open channel at its excitation energy, the event is left
 * untouched: the kind stays as set by the abrasion stage and the fragment
 * remains in event_out->fragments[] for the fragment-pool counter (its E* is
 * dropped as un-modelled gamma de-excitation — a documented limitation).
 *
 * Otherwise the fragment slot is consumed and the de-excitation chain is
 * resolved with a fixed-size work stack:
 *   - whitelist nuclides (n, p, d, t, He-3, He-4) without open channels are
 *     appended to event_out->secondaries[],
 *   - nuclides with an open channel are split further (two-body kinematics,
 *     isotropic in the parent rest frame, Lorentz-boosted to the lab),
 *   - anything else is appended to event_out->fragments[] as unprocessed.
 *
 * Truncation policy: once secondaries[] is full, remaining products are
 * routed to fragments[] (counted as unprocessed rather than silently losing
 * energy); if fragments[] is also full, only energy accounting is lost and
 * the product is dropped.
 *
 * event_out->kind is set to OSH_NUCLEAR_EVENT_FRAGMENTATION only when at
 * least one transportable secondary was emitted.
 *
 * @param model      Compiled model (must not be NULL).
 * @param fragment   Prefragment to de-excite; may alias event_out->fragments[0]
 *                   (it is copied before the event is modified).
 * @param rng        RNG state (mutated; only side effect besides event_out).
 * @param event_out  Event to append products to (secondaries/fragments/kind).
 */
void osh_nuclear_fermi_breakup_step(struct osh_nuclear_fermi_breakup const *model,
                                    struct osh_nuclear_fragment const *fragment,
                                    struct osh_rng *rng,
                                    struct osh_nuclear_event *event_out);

#ifdef __cplusplus
}
#endif

#endif /* OSH_NUCLEAR_FERMI_BREAKUP_H */
