#ifndef OSH_NUCLEAR_COMPOUND_H
#define OSH_NUCLEAR_COMPOUND_H

/**
 * @file osh_nuclear_compound.h
 * @brief Compound nucleus adapter — routes to Fermi break-up or heavy-A sink.
 *
 * @details
 * Called by the neutron reaction layer after constructing the compound nucleus
 * (Z, A+1, E*) from a neutron capture event (n,n'), (n,2n), or generic
 * non-elastic.  Bridges the reaction layer to the Fermi break-up back-end
 * without exposing FBU types to the transport layer.
 *
 * Routing logic:
 *   A ≤ OSH_FERMI_BREAKUP_AMAX  →  osh_nuclear_fermi_breakup_step()
 *   A >  OSH_FERMI_BREAKUP_AMAX  →  heavy-A sink: OSH_NUCLEAR_EVENT_ABSORB,
 *                                    energy deposited locally.  One-time
 *                                    OSH_DIAG_INFOF per (Z,A) on first hit.
 *                                    Future entry point for SMM.
 *
 * Output: struct osh_nuclear_event (defined in osh_nuclear_handler.h).
 * The neutron reaction layer translates this into its own event struct.
 */

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

struct osh_diag_sink;
struct osh_nuclear_event;
struct osh_nuclear_fermi_breakup;
struct osh_rng;

/**
 * @brief De-excite a compound nucleus into final-state secondaries.
 *
 * @param[in]  z          Compound nucleus atomic number.
 * @param[in]  a          Compound nucleus mass number (= target A + 1).
 * @param[in]  e_star_mev Excitation energy E* [MeV].
 * @param[in]  p_lab_mev  Lab-frame momentum vector [MeV/c] (3 elements).
 * @param[in]  fbu        Compiled Fermi break-up model (borrowed).
 * @param[in]  diag       Diagnostic sink; may be NULL.
 * @param[in]  rng        RNG state (mutated in place).
 * @param[out] event_out  Filled with secondaries or ABSORB for heavy-A sink.
 */
void osh_nuclear_compound_step(unsigned int z,
                               unsigned int a,
                               double e_star_mev,
                               double const p_lab_mev[3],
                               struct osh_nuclear_fermi_breakup const *fbu,
                               struct osh_diag_sink const *diag,
                               struct osh_rng *rng,
                               struct osh_nuclear_event *event_out);

#ifdef __cplusplus
}
#endif

#endif /* OSH_NUCLEAR_COMPOUND_H */
