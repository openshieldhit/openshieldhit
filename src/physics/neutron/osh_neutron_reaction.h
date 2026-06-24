#ifndef OSH_NEUTRON_REACTION_H
#define OSH_NEUTRON_REACTION_H

/**
 * @file osh_neutron_reaction.h
 * @brief Neutron reaction sampling — target selection, channel sampling, final state.
 *
 * @details
 * Middle layer between the transport loop and the nuclear back-ends (FBU,
 * compound adapter, kinematics).  The transport loop owns path-length sampling
 * and pool management; this layer owns all reaction physics.
 *
 * Channel routing:
 *   σ_el        → ELASTIC: update neutron dir/energy; H-1 recoil proton secondary
 *   σ(n,γ)      → CAPTURE: kill neutron, deposit T_n locally
 *   σ(n,p/α)    → CHARGE_EXCHANGE: 2-body relativistic kinematics, ion secondary
 *   σ(n,n')+σ(n,2n) → COMPOUND: build (Z,A+1,E*) compound, route through
 *                       osh_nuclear_compound_step() → FBU or heavy-A sink
 *   remainder   → LOCAL_DEPOSIT (Tier-2 generic non-elastic, heavy-A sink output)
 *
 * Ion secondaries from CHARGE_EXCHANGE and COMPOUND are always generated;
 * the transport layer applies the ion-feedback flag when deciding whether to
 * push them to the ion pool.
 */

#include <stddef.h>

#include "physics/nuclear/osh_nuclear_handler.h"

#ifdef __cplusplus
extern "C" {
#endif

struct osh_neutron_xsec;
struct osh_nuclear_fermi_breakup;
struct osh_rng;

/** Outcome kind of a neutron reaction event. */
enum osh_neutron_reaction_kind {
    OSH_NEUTRON_REACTION_NONE = 0,
    OSH_NEUTRON_REACTION_ELASTIC,        /**< Elastic scatter; neutron continues with new dir/e */
    OSH_NEUTRON_REACTION_CAPTURE,        /**< (n,γ): kill neutron, deposit T_n locally          */
    OSH_NEUTRON_REACTION_CHARGE_EXCHANGE,/**< (n,p) or (n,α): ion secondary, neutron killed     */
    OSH_NEUTRON_REACTION_COMPOUND,       /**< (n,n')/(n,2n)/generic: FBU or heavy-A sink        */
    OSH_NEUTRON_REACTION_LOCAL_DEPOSIT   /**< Generic non-elastic with no tracked secondaries    */
};

/**
 * @brief Result of one neutron reaction.
 *
 * For ELASTIC the transport loop updates the neutron in-place with
 * neutron_dir/neutron_e_mev.  For all other kinds the neutron is killed.
 * Secondaries (from FBU or CHARGE_EXCHANGE) are in secondaries[].
 */
struct osh_neutron_reaction_event {
    enum osh_neutron_reaction_kind kind;
    double local_deposit_mev;                                      /**< Energy to deposit locally    */
    double neutron_dir[3];                                         /**< Scattered neutron direction   */
    double neutron_e_mev;                                          /**< Scattered neutron energy [MeV]*/
    size_t n_secondaries;
    struct osh_nuclear_secondary secondaries[OSH_NUCLEAR_MAX_SECONDARIES];
};

/**
 * @brief Sample one neutron reaction in the current material cell.
 *
 * Selects target nuclide proportional to nᵢ σ_tot,i, samples a channel,
 * and generates the final state.  Both the xsec model (warning tracker)
 * and the RNG are mutated.
 *
 * @param[in,out] xsec         Cross-section model (Tier-2 warning state mutated).
 * @param[in]     handler      Compiled nuclear handler (element cache + FBU).
 * @param[in]     material_idx Dense material index into handler element tables.
 * @param[in]     rho_g_cm3   Material density [g/cm³].
 * @param[in]     e_mev        Neutron kinetic energy [MeV].
 * @param[in]     dir          Neutron unit direction (length 3).
 * @param[in,out] rng          RNG state.
 * @param[out]    event_out    Filled on return; kind=NONE only on empty material.
 */
void osh_neutron_reaction_sample(struct osh_neutron_xsec *xsec,
                                 struct osh_nuclear_handler const *handler,
                                 size_t material_idx,
                                 double rho_g_cm3,
                                 double e_mev,
                                 double const dir[3],
                                 struct osh_rng *rng,
                                 struct osh_neutron_reaction_event *event_out);

#ifdef __cplusplus
}
#endif

#endif /* OSH_NEUTRON_REACTION_H */
