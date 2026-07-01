#ifndef OSH_TRANSPORT_ION_H
#define OSH_TRANSPORT_ION_H

#include <stddef.h> /* size_t */

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
 * @brief Transport the ion primaries of one history range [@p hist_lo, @p hist_hi).
 *
 * @details
 * The range-aware core of the ion pass: it validates the run parameters, borrows
 * the simulation's pre-allocated ion pool and geometry scratch, and transports
 * exactly the primaries in [@p hist_lo, @p hist_hi) to termination.  Physics: CSDA
 * energy loss, Highland/Molière MCS (random-hinge method), Bohr Gaussian energy
 * straggling, plus the nuclear secondaries banked for the family scheduler.  Each
 * primary is seeded from its *global* history index, so splitting [0, nstat) into
 * disjoint sub-ranges and transporting them in any order reproduces the canonical
 * per-history streams (scored output invariant up to floating-point reduction
 * order).  This is the seam the batch-aware scheduler drives once per checkpoint
 * batch (issue #195); a whole-run pass is simply the range [0, nstat).
 *
 * The random-hinge treatment follows the fast proton-transport approach described
 * by Fippel and Soukup (Med Phys. 2004;31(8):2263-2273. doi:10.1118/1.1769631).
 * The caller must invoke osh_gemca_compile() before and osh_gemca_runtime_free()
 * after this function.
 *
 * The borrowed ion pool is reset to empty on entry, so a range call assumes no
 * live carry-over from a previous batch (the wavefront loop always drains the
 * pool to empty before returning).
 *
 * @param[in,out] transport_ctx  Per-run transport context (pools, params, diag).
 * @param[in]     beam_rt        Hot beam runtime for primary generation.
 * @param[in]     geom_rt        Compiled geometry runtime.
 * @param[in]     material_rt    Hot material runtime tables.
 * @param[in,out] score_rt       Scoring runtime.
 * @param[in]     hist_lo        Inclusive lower bound of the range.
 * @param[in]     hist_hi        Exclusive upper bound; must be > @p hist_lo.
 * @param[out]    completed_out  Receives the number of primaries whose histories
 *                               finished in this range (== hist_hi - hist_lo
 *                               unless a clean stop drained it early).  May be NULL.
 *
 * @returns OSH_OK on success, OSH_EINVAL on a bad argument/range, or an OSH_E*
 *          from validation or the wavefront loop.
 */
enum osh_status osh_transport_ion_run_range(struct osh_transport_context *transport_ctx,
                                            struct osh_beam_runtime *beam_rt,
                                            struct osh_gemca_runtime const *geom_rt,
                                            struct osh_material_runtime const *material_rt,
                                            struct osh_scoring_runtime *score_rt,
                                            size_t hist_lo,
                                            size_t hist_hi,
                                            size_t *completed_out);

#ifdef __cplusplus
}
#endif

#endif /* OSH_TRANSPORT_ION_H */
