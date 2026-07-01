#ifndef OSH_CHECKPOINT_POLICY_H
#define OSH_CHECKPOINT_POLICY_H

#include <stddef.h> /* size_t */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @file osh_checkpoint_policy.h
 * @brief The checkpoint / partial-results policy — issue #195.
 *
 * @details
 * There is exactly one dial in this design: *how often the whole simulation is
 * brought to a family-complete, quiescent point*.  Call it a **checkpoint**.  A
 * checkpoint is where (optionally) every transport family is drained, parallel
 * workers merge into the master (#161), one variance batch folds (#169), and an
 * (optional) dump fires (#193).  A **batch** `[b, b+K)` is the primary range
 * transported between two checkpoints, so the checkpoint barrier and the
 * parallel-partition barrier are the same mechanism.
 *
 * @b Why quiescence matters. A scoring snapshot is only *physically complete*
 * once every transport family (ions, then the neutrons/fragments they banked, …)
 * spawned by the completed primaries has been drained into scoring.  `DLET` is
 * `Σ(LET·dose)/Σdose`; dropping the neutron/fragment numerator *and* denominator
 * shifts the ratio, so an ion-only partial is *biased*, not merely noisier.  A
 * checkpoint is therefore the correct place to observe a partial result.
 *
 * @b The one constraint on cadence. A **wall-time** cadence self-bounds the
 * per-run overhead independent of core count ("every 10 min" costs the drain +
 * barrier + merge once per 10 min whether on 4 or 200 cores).  A **count**
 * cadence does not (more cores → faster → more checkpoints per wall-second).  So
 * time cadence is for production/parallel and count cadence is for deterministic
 * tests (#168): count-cadence batches over fixed sub-ranges are order-independent
 * because each history's RNG stream is a pure function of its global index.
 *
 * @b This module is data + arithmetic only. It carries the intent (the struct)
 * and computes the next batch size (@ref osh_checkpoint_next_batch_size); it
 * never touches files, clocks, threads, or scoring.  The transport scheduler
 * consumes it to size its outer loop; the file/variance/merge machinery that
 * hangs off a checkpoint lands in #193 / #169 / #161.
 *
 * @b Default = today. A zero-initialised policy (or a NULL policy pointer) is
 * @ref OSH_PARTIAL_NONE — FINAL-ONLY: one batch of `K = nstat`, the fastest
 * path, byte-for-byte identical to a run configured without a policy.
 */

/** Partial-results mode: how often the run reaches a quiescent checkpoint. */
enum osh_partial_mode {
    OSH_PARTIAL_NONE = 0, /**< FINAL-ONLY: one batch (K = nstat).  Fastest; = today. */
    OSH_PARTIAL_LIVE = 1  /**< LIVE previews: family-complete batches at a cadence. */
};

/** Completeness of a partial output — the honesty label for a checkpoint. */
enum osh_partial_completeness {
    OSH_PARTIAL_EXACT = 0, /**< Checkpoint drained all families → physically complete. */
    OSH_PARTIAL_APPROX = 1 /**< Dumped at the inner ion safe point, families pending → labelled. */
};

/**
 * @brief The single knob: batch cadence + completeness for one run.
 *
 * @details
 * The presence of a cadence selects the mode; an all-zero struct is FINAL-ONLY.
 * @ref every_s (time) and @ref every_primaries (count) are the two cadences; at
 * most one is expected to be set.  @ref batch decouples the batch size from the
 * cadence (e.g. batch for variance, write rarely); 0 means "auto — derive the
 * batch from the cadence × measured throughput".  @ref write_files says whether
 * a checkpoint also dumps or checkpoints silently (both consumed by #193).
 *
 * The struct shape is stable so downstream issues can grow into it without an
 * API churn — mirroring how @ref osh_run_control declared its dump-cadence
 * scalars ahead of their consumer.
 */
struct osh_checkpoint_policy {
    enum osh_partial_mode mode;                 /**< NONE = final-only (default).                         */
    enum osh_partial_completeness completeness; /**< EXACT drains families; APPROX labels (consumed #193).*/
    double every_s;                             /**< Wall-time cadence [s]; 0 = off.                      */
    size_t every_primaries;                     /**< Count cadence (primaries); 0 = off.                  */
    size_t batch;                               /**< Explicit batch size; 0 = auto (cadence × rate).      */
    int write_files;                            /**< Dump at each checkpoint vs. checkpoint silently (#193).*/
};

/**
 * @brief Zero-initialise a policy to FINAL-ONLY / EXACT (today's behaviour).
 *
 * @param[out] policy  Policy to initialise; no-op when NULL.
 */
void osh_checkpoint_policy_init(struct osh_checkpoint_policy *policy);

/**
 * @brief Is this run final-only (one batch of K = nstat, zero checkpoint overhead)?
 *
 * @param[in] policy  Policy, or NULL.
 * @returns 1 when @p policy is NULL or its mode is @ref OSH_PARTIAL_NONE, else 0.
 */
int osh_checkpoint_policy_is_final_only(struct osh_checkpoint_policy const *policy);

/**
 * @brief Size the next batch `[done, done + K)` for the scheduler's outer loop.
 *
 * @details
 * FINAL-ONLY (or NULL) returns @p remaining — the `K = nstat` fast path that
 * keeps the run a single pass.  In LIVE mode the batch size is, in priority
 * order: the explicit @ref osh_checkpoint_policy::batch; else the count cadence
 * @ref osh_checkpoint_policy::every_primaries; else the adaptive time cadence
 * `round(measured_rate_pps × every_s)` so "preview every 10 min" resolves to the
 * same wall interval on a laptop or a 200-core node.  The result is always
 * clamped to `[1, remaining]`.  When LIVE is requested but no cadence is usable
 * yet (e.g. the first batch of a time cadence, before any throughput has been
 * measured), the whole @p remaining is returned so the run still makes progress
 * and stays exact.
 *
 * @param[in] policy             Policy, or NULL (⇒ final-only).
 * @param[in] measured_rate_pps  Recent throughput [primaries/s] for the adaptive
 *                               time cadence; ignored unless a time cadence is
 *                               active and neither @c batch nor @c every_primaries
 *                               is set.  Pass 0 when unknown.
 * @param[in] remaining          Primaries left to transport (> 0 expected).
 * @returns The batch size K in `[1, remaining]`, or 0 when @p remaining is 0.
 */
size_t
osh_checkpoint_next_batch_size(struct osh_checkpoint_policy const *policy, double measured_rate_pps, size_t remaining);

/**
 * @brief Stable string label for a completeness value.
 *
 * @details
 * The honesty stamp for a partial output (issue #195 "Completeness labelling"):
 * BDO tag @c OSHBDO_RT_COMPLETENESS and the ASCII header @c "# COMPLETENESS:".
 * EXACT → "exact"; APPROX → "families_pending".  The file writers (#193) append
 * the pending family list; this is the leading token.
 *
 * @param[in] completeness  Completeness value.
 * @returns A static, never-NULL string.
 */
char const *osh_checkpoint_completeness_label(enum osh_partial_completeness completeness);

#ifdef __cplusplus
}
#endif

#endif /* OSH_CHECKPOINT_POLICY_H */
