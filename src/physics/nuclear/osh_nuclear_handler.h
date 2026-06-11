#ifndef OSH_NUCLEAR_HANDLER_H
#define OSH_NUCLEAR_HANDLER_H

/**
 * @file osh_nuclear_handler.h
 *
 * @brief Nuclear interaction handler: per-material composition cache and
 *        per-step event selection for all nuclear channels.
 *
 * @details
 * The handler compiles per-material element composition from the material
 * workspace at startup, then dispatches the appropriate nuclear physics at
 * each transport step.
 *
 * Currently supported channels:
 *   - pp elastic scattering (hydrogen target, projectile = proton)
 *   - Tripathi inelastic absorption (all ion projectiles)
 *   - abrasion + Fermi break-up de-excitation (proton projectile)
 *
 * The step entry point is pool-independent: it takes scalar inputs and writes
 * a result struct; only RNG state is mutated as a side effect.  This makes the
 * nuclear phase trivially separable from atomic physics for future wavefront
 * batching or GPU offload.
 */

#include <stddef.h>

#include "openshieldhit/status.h"
#include "physics/nuclear/osh_nuclear_fermi_breakup.h"

#ifdef __cplusplus
extern "C" {
#endif

struct osh_material_workspace;
struct osh_rng;
struct osh_transport_params;
struct particle;

/* ---- Event result -------------------------------------------------------- */

/** Maximum secondaries per nuclear event (sized for Fermi breakup / SMM). */
#define OSH_NUCLEAR_MAX_SECONDARIES 32

/** Maximum residual nuclear fragments per event: the abrasion prefragment
 *  plus unprocessed Fermi break-up residues (non-whitelist nuclides without
 *  open channels, or truncation overflow). */
#define OSH_NUCLEAR_MAX_FRAGMENTS 4

/** Classification of the nuclear event that fired on a given step. */
enum osh_nuclear_event_kind {
    OSH_NUCLEAR_EVENT_NONE = 0,     /**< No nuclear event this step.                          */
    OSH_NUCLEAR_EVENT_ABSORB,       /**< Inelastic kill, no secondaries (threshold / fallback). */
    OSH_NUCLEAR_EVENT_ELASTIC_PP,   /**< pp elastic scatter.                                  */
    OSH_NUCLEAR_EVENT_ABRASION,     /**< Inelastic: primary absorbed, fast nucleons emitted (AA model). */
    OSH_NUCLEAR_EVENT_FRAGMENTATION /**< Inelastic: abrasion + Fermi break-up products emitted. */
};

/** One secondary particle produced by a nuclear event. */
struct osh_nuclear_secondary {
    double dir[3];                  /**< Exit unit direction.                    */
    double energy;                  /**< Kinetic energy [MeV].                   */
    struct particle const *species; /**< Borrowed from particle registry.        */
};

/** One residual nuclear fragment produced by an inelastic event. */
struct osh_nuclear_fragment {
    double excitation_energy; /**< Excitation energy E* [MeV].             */
    double p[3];              /**< Lab momentum vector [MeV/c].            */
    unsigned int z;           /**< Residual atomic number.                 */
    unsigned int a;           /**< Residual mass number.                   */
};

/**
 * @brief Result of one nuclear step.
 *
 * @details
 * Written by osh_nuclear_handler_step(); read by the transport commit phase.
 * When kind == OSH_NUCLEAR_EVENT_NONE the remaining fields are undefined.
 */
struct osh_nuclear_event {
    enum osh_nuclear_event_kind kind;
    double primary_energy; /**< Primary exit KE [MeV]; 0 if ABSORB.       */
    double primary_dir[3]; /**< Primary exit direction; unchanged if ABSORB. */
    size_t n_secondaries;
    struct osh_nuclear_secondary secondaries[OSH_NUCLEAR_MAX_SECONDARIES];
    size_t n_fragments;
    struct osh_nuclear_fragment fragments[OSH_NUCLEAR_MAX_FRAGMENTS];
};

/* ---- Handler ------------------------------------------------------------- */

/** One element entry in the per-material composition cache. */
struct osh_nuclear_elem {
    unsigned int z;      /**< Atomic number.                                  */
    unsigned int a;      /**< Approximate mass number.                        */
    float mass_fraction; /**< Element mass fraction in this material [0, 1].  */
};

/**
 * @brief Per-material element composition cache.
 *
 * @details
 * Compiled once from the material workspace.  The flat elem_pool array stores
 * all materials' elements contiguously; elem_offset[i] and elem_count[i] index
 * into it for material i.
 *
 * Owned by the transport context; lifetime must exceed the transport run.
 * Free with osh_nuclear_handler_free().
 */
struct osh_nuclear_handler {
    struct osh_nuclear_fermi_breakup fbu; /**< Fermi break-up channel table.            */
    struct osh_nuclear_elem *elem_pool;   /**< Flat array: all materials, all elements. */
    size_t *elem_offset;                  /**< elem_offset[i]: start of material i.    */
    size_t *elem_count;                   /**< elem_count[i]:  element count.          */
    size_t nmaterials;
};

/**
 * @brief Compile the nuclear handler from the material workspace.
 *
 * @details
 * Allocates and populates the element composition cache.  Must be called after
 * the material workspace is fully populated and before any transport step.
 *
 * @param[in]  ws   Material workspace (read-only).
 * @param[out] out  Handler to populate (must point to zero-initialised storage).
 * @returns OSH_OK on success, OSH_ENOMEM on allocation failure.
 */
enum osh_status osh_nuclear_handler_compile(struct osh_material_workspace const *ws, struct osh_nuclear_handler *out);

/**
 * @brief Free all memory owned by the handler.
 *
 * @param h  Handler to free.  May be NULL.  Fields are zeroed after free.
 */
void osh_nuclear_handler_free(struct osh_nuclear_handler *h);

/**
 * @brief Sample nuclear interactions over one transport step.
 *
 * @details
 * Determines whether a nuclear event fires and, if so, which channel.
 * Uses competing independent Poisson processes:
 *
 *   rate_tot = rate_inel + rate_pp
 *   p_event  = 1 − exp(−ds × rate_tot)
 *   p_pp|event = rate_pp / rate_tot   (exact for exponential hazards)
 *
 * Two energies are distinguished to avoid double-counting energy loss:
 *   @p rate_energy_mev  (= ctx->e0, pre-step) is used for cross-section lookup
 *   @p final_energy_mev (= ctx->exit_energy, post-CSDA/straggling) is used for
 *                        final-state kinematics; ensures
 *                        deposited + e1_out + e2_out == e0 within rounding.
 *
 * This function is pool-independent (scalar in → result struct out).
 * Only RNG state is mutated as a side effect.
 *
 * @param handler           Compiled handler (must not be NULL).
 * @param rate_energy_mev   Pre-step KE used for cross-section [MeV].
 * @param final_energy_mev  Post-CSDA/straggling KE used for kinematics [MeV].
 * @param incident_dir      Pre-step unit direction (length 3).
 * @param material_idx      Index into handler->elem_offset/count.
 * @param ds_gcm2           Step areal density [g/cm²].
 * @param projectile        Projectile species (borrowed).
 * @param params            Transport control parameters.
 * @param rng               RNG state (mutated).
 * @param event_out         Event result (always written; kind=NONE if no event).
 */
void osh_nuclear_handler_step(struct osh_nuclear_handler const *handler,
                              double rate_energy_mev,
                              double final_energy_mev,
                              double const incident_dir[3],
                              size_t material_idx,
                              double ds_gcm2,
                              struct particle const *projectile,
                              struct osh_transport_params const *params,
                              struct osh_rng *rng,
                              struct osh_nuclear_event *event_out);

#ifdef __cplusplus
}
#endif

#endif /* OSH_NUCLEAR_HANDLER_H */
