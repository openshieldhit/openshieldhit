#ifndef OSH_TRANSPORT_NEUTRON_H
#define OSH_TRANSPORT_NEUTRON_H

#include "openshieldhit/status.h"

#ifdef __cplusplus
extern "C" {
#endif

struct osh_transport_context;
struct osh_beam_runtime;
struct osh_gemca_runtime;
struct osh_material_runtime;
struct osh_scoring_runtime;
struct osh_score_target;

/**
 * @brief Minimal fast-neutron transport loop.
 *
 * @details
 * Drains the neutron pool supplied via transport_ctx->neutron_pool using a
 * straight-line geometric transport (no CSDA, no MCS):
 *
 *   1. Wavefront geometry batch: zone-ref + boundary-distance for all live slots.
 *   2. For each live neutron: energy cutoff, escape/blackhole kill, free-path
 *      sampling (exponential in 1/Σ_tot), boundary crossing or interaction.
 *   3. Interaction via osh_neutron_reaction_sample(): elastic, capture, (n,p),
 *      (n,α), compound (FBU or heavy-A sink).
 *   4. Neutron secondaries from compound events are pushed back to the pool;
 *      ion secondaries are deposited locally (ion-feedback not yet wired).
 *   5. Dead slots are compacted after each pass; new secondaries are processed
 *      in the next pass.
 *
 * Energy deposits from captures and recoils are currently not scored
 * (osh_scoring_score_point() is not yet implemented).
 *
 * @param transport_ctx  Must have neutron_pool and nuclear_handler set.
 * @param beam_rt        Unused (no primary beam refill for neutrons).
 * @param score_rt       Scoring runtime the step deposits into.
 * @param target         Caller-owned deposit target (accumulator set + scratch);
 *                       NULL (or NULL fields) falls back to @p score_rt's master
 *                       views, so the single-worker path is unchanged (issue #230).
 */
enum osh_status osh_transport_neutron_run(struct osh_transport_context *transport_ctx,
                                          struct osh_beam_runtime *beam_rt,
                                          struct osh_gemca_runtime const *geom_rt,
                                          struct osh_material_runtime const *material_rt,
                                          struct osh_scoring_runtime *score_rt,
                                          struct osh_score_target const *target);

#ifdef __cplusplus
}
#endif

#endif /* OSH_TRANSPORT_NEUTRON_H */
