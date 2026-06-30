#ifndef OSH_NUCLEAR_FERMI_BREAKUP_H
#define OSH_NUCLEAR_FERMI_BREAKUP_H

/**
 * @file osh_nuclear_fermi_breakup.h
 * @brief Fermi break-up de-excitation of light excited prefragments.
 *
 * @details
 * De-excites the residual prefragment (Z, A, E*) produced by the abrasion
 * stage into nucleons and light ions via the simultaneous N-body Fermi
 * break-up (Fermi 1950), applicable to light nuclei A <= 16.
 *
 * The model enumerates at compile time all N-body partitions (N = 2 ..
 * OSH_FERMI_BREAKUP_NMAX) of the parent (Z, A) into ground-state nuclear
 * species drawn from the full isotope database.  Particle-unstable members of
 * a partition (He-5, Li-5, Be-8, ...) are decayed further at runtime by a work
 * stack; bound non-transportable residues (Li-6, Li-7, Be-9, ...) are returned
 * as unprocessed fragments.
 *
 * At runtime, partition i is weighted by the canonical microcanonical
 * Fermi break-up phase-space probability
 *
 *   W_i = [V / ((2*pi)^(3/2) * hbar^3)]^(N_i - 1)
 *         * (prod g_k / prod n_j!) * (prod m_k)^(3/2)
 *         / Gamma(3 N_i / 2 - 3/2)
 *         * (E* + Q_i)^(3 N_i / 2 - 5/2)
 *
 * where the leading factor is the free-volume / density-of-states term that
 * sets the multiplicity scale and makes the weights of different N
 * dimensionally commensurable.  V = (4/3) pi r0^3 A is the break-up volume
 * (single physical knob r0, calibrated against Geant4 FermiBreakUp; see
 * OSH_FERMI_BREAKUP_R0_FM).  The remaining factors encode spin degeneracies,
 * identical-particle combinatorics, fragment masses, and the Gamma
 * normalisation.  At high E* the higher-N channels grow faster, driving the
 * multiplicity to rise as observed in the Geant4 reference data.
 *
 * Kinematics: two-body decay (osh_kinematics) for N=2 partitions; Kopylov
 * N-body phase space for N>=3 partitions.  Every product is fed back through
 * the work stack, so a single statistical picture applies at all N.
 *
 * Final-product policy (strict): only the whitelist {n, p, d, t, He-3, He-4}
 * is ever emitted as transportable secondaries.  Any other nuclide is returned
 * as an unprocessed fragment in the event.
 *
 * TODO(#196): fragments are produced in their ground state only.  Systematic
 * particle-unstable excited-state level tables (which fill the residual
 * high-E* multiplicity tail seen in G4FermiBreakUp) are a planned follow-up.
 *
 * The partition table and product species descriptors are compiled once at
 * startup from the isotope database; the per-event step performs no
 * allocations and only mutates the RNG state.
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

/** Maximum fragment multiplicity N enumerated per partition. */
#define OSH_FERMI_BREAKUP_NMAX 6

/** Effective break-up volume radius r0 [fm]; the sole multiplicity calibration
 *  knob (V = (4/3) pi r0^3 A).  Calibrated to G4FermiBreakUp mean multiplicity
 *  vs excitation per nucleon over C-12/C-13/N-12/N-13 (RMS ~0.37).  It is
 *  smaller than the geometric nuclear radius (~1.2 fm) because the
 *  ground-state-only approximation routes all of E*+Q into kinetic energy and
 *  would otherwise over-favour high multiplicity; the effective free volume
 *  absorbs that bias.  A future excited-state treatment (TODO) would restore a
 *  more physical value. */
#define OSH_FERMI_BREAKUP_R0_FM 0.50

/** One N-body decay partition: parent -> N fragments listed in fspec_pool. */
struct osh_fermi_partition {
    double weight_prefactor; /**< C^(N-1) * prod(g_i)/prod(n_k!) * prod(m_i^(3/2)) / Gamma(3N/2-3/2).
                              where C = V / ((2*pi)^(3/2) * (hbar*c)^3). */
    uint32_t fspec_offset;   /**< First index in fspec_pool for this partition's (z,a) list. */
    float q_mev;             /**< Q = M_parent_gs - sum(m_i_gs) [MeV]. */
    uint8_t n_frags;         /**< N: number of fragments. */
    uint8_t _pad[3];
};

/** Fragment species entry: (z, a) of one ground-state fragment in a partition. */
struct osh_fermi_frag_spec {
    uint8_t z;
    uint8_t a;
};

/**
 * @brief Compiled Fermi break-up model: partition table and product species.
 *
 * @details
 * All per-parent arrays are indexed by the dense (z, a) key z*(AMAX+1)+a for
 * z in [0, ZMAX], a in [0, AMAX].  Compiled once by
 * osh_nuclear_fermi_breakup_compile(); owned by the nuclear handler.
 * Free with osh_nuclear_fermi_breakup_free().
 */
struct osh_nuclear_fermi_breakup {
    double *mass_mev;                       /**< Ground-state nuclear masses per dense (z,a). */
    struct particle *species;               /**< Final-product descriptors (n, p, d, t, He-3, He-4). */
    struct osh_fermi_partition *part_pool;  /**< Flat array of all parents' partitions. */
    struct osh_fermi_frag_spec *fspec_pool; /**< Flat array of (z,a) fragment specs. */
    uint32_t *part_offset;                  /**< part_offset[idx]: start of parent idx in part_pool. */
    uint16_t *part_count;                   /**< part_count[idx]: partition count of parent idx. */
    size_t npartitions;                     /**< Total entries in part_pool. */
    size_t nfspecs;                         /**< Total entries in fspec_pool. */
};

/**
 * @brief Compile the break-up partition table from the isotope database.
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
 * If the fragment is outside the model domain or has no open channel at its
 * excitation energy, the event is left untouched.
 *
 * Otherwise the fragment slot is consumed and all N products of the sampled
 * partition are resolved: whitelist nuclides go to event_out->secondaries[],
 * other nuclides go to event_out->fragments[] as unprocessed residues.
 * N=2 partitions into particle-unstable products (He-5, Li-5, Be-8) are
 * handled via a small work stack that decays them further.
 *
 * event_out->kind is set to OSH_NUCLEAR_EVENT_FRAGMENTATION only when at
 * least one transportable secondary was emitted.
 *
 * @param model      Compiled model (must not be NULL).
 * @param fragment   Prefragment to de-excite; may alias event_out->fragments[0].
 * @param rng        RNG state (mutated).
 * @param event_out  Event to append products to.
 */
void osh_nuclear_fermi_breakup_step(struct osh_nuclear_fermi_breakup const *model,
                                    struct osh_nuclear_fragment const *fragment,
                                    struct osh_rng *rng,
                                    struct osh_nuclear_event *event_out);

#ifdef __cplusplus
}
#endif

#endif /* OSH_NUCLEAR_FERMI_BREAKUP_H */
