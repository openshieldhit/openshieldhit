#ifndef OSH_TRANSPORT_PHOTON_H
#define OSH_TRANSPORT_PHOTON_H

#include "openshieldhit/status.h"

#ifdef __cplusplus
extern "C" {
#endif

struct osh_transport_context;
struct osh_beam_runtime;
struct osh_gemca_runtime;
struct osh_material_runtime;
struct osh_scoring_runtime;

/**
 * @brief Photon transport loop (stub — not yet implemented).
 *
 * @details
 * Placeholder family kernel used by the scheduler so unsupported photon work
 * reports a normal OSH_ENOTSUP status through the same dispatch path as other
 * particle families.
 *
 * @param[in,out] transport_ctx  Per-run transport context; diagnostics are used
 *                               for the unsupported-family message.
 * @param[in]     beam_rt        Unused until photon transport is implemented.
 * @param[in]     geom_rt        Unused until photon transport is implemented.
 * @param[in]     material_rt    Unused until photon transport is implemented.
 * @param[in,out] score_rt       Unused until photon transport is implemented.
 *
 * @returns OSH_ENOTSUP always.
 */
enum osh_status osh_transport_photon_run(struct osh_transport_context *transport_ctx,
                                         struct osh_beam_runtime *beam_rt,
                                         struct osh_gemca_runtime const *geom_rt,
                                         struct osh_material_runtime const *material_rt,
                                         struct osh_scoring_runtime *score_rt);

#ifdef __cplusplus
}
#endif

#endif /* OSH_TRANSPORT_PHOTON_H */
