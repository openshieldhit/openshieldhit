#ifndef OSH_WORKER_CONTEXT_H
#define OSH_WORKER_CONTEXT_H

#include <stddef.h>

#include "openshieldhit/status.h"

#ifdef __cplusplus
extern "C" {
#endif

struct osh_particle_pool;
struct osh_zone_ref;
struct osh_scoring_accumulator;
struct osh_scoring_scratch;
struct osh_transport_profile;

/**
 * @brief Everything one transport worker needs to run an assigned slice of histories.
 *
 * @details
 * Bundles the *transport-local* mutable state — the particle pool, per-step
 * batch scratch, and the run-profile counters — for one slice of the run, so
 * this state stops being a barrier to partitioning work across workers (threads,
 * MPI ranks, or sequential profiling replicas) over disjoint history ranges.  It
 * is one building block, not the whole story: the scoring accumulators in
 * score_rt are owned elsewhere and remain shared, so concurrent execution will
 * additionally need per-worker scoring state (the @ref accumulators field below
 * is where per-worker scoring memory will attach).  The profile counters are no
 * longer in that list — each worker accumulates into its own @ref profile and
 * the driver merges them (see osh_transport_profile_merge).  Beam primary
 * generation is no longer in that list either — a worker
 * fills its pool through osh_beam_runtime_fill_pool_at() with an explicit global
 * base from its own range, so the beam runtime is shared read-only.
 *
 * The assigned work is the half-open history range [@ref hist_lo, @ref hist_hi).
 * Per-history RNG seeding is derived from the global history index
 * (rndoffset + history), so disjoint ranges produce disjoint, non-overlapping
 * random streams and the scored result is independent of how the range is split.
 *
 * In the current single-worker build exactly one context covers the whole run,
 * [0, nstat), and deposits go straight into the shared master accumulators.
 */
struct osh_worker_context {
    /* Assigned half-open history range. */
    size_t hist_lo; /* First global history index (inclusive). */
    size_t hist_hi; /* One past the last global history index (exclusive). */

    /* Transport-local pool and geometry-query scratch for this worker.  In the
     * single-worker baseline these are borrowed from the simulation's
     * pre-allocated buffers (see osh_worker_context_attach); a future parallel
     * worker would instead own a private set.  Never resized on the hot path. */
    size_t capacity;                /* Pool capacity (a pure performance knob). */
    struct osh_particle_pool *pool; /* Live primaries/secondaries for this worker (borrowed). */
    struct osh_zone_ref *zone_refs; /* Batch zone-reference scratch, length >= @ref capacity (borrowed). */
    double *dist_batch;             /* Batch boundary-distance scratch, length >= @ref capacity (borrowed). */

    /* Per-worker private scoring memory (the deposit target).  When @ref
     * accumulators is non-NULL the worker deposits into this private set and the
     * driver folds it into the shared master with osh_scoring_accumulator_merge()
     * once the range completes; @ref scratch is the matching per-worker traversal
     * scratch the deposit path writes crossings into.  Both NULL means "use the
     * shared master views" — the single-worker baseline, bit-for-bit unchanged.
     * The sequential-replica driver (issue #230) is the first caller to populate
     * these; a threaded/MPI/WASM backend attaches its own set here the same way. */
    struct osh_scoring_accumulator *accumulators;
    size_t naccumulators;
    struct osh_scoring_scratch *scratch;

    /* Per-worker transport profile (phase timers + event counters).  The worker
     * accumulates here on the hot path, never into a shared profile, so workers
     * do not race on these counters; the driver folds them into the run's master
     * profile with osh_transport_profile_merge() once the range completes.  NULL
     * disables profiling.  In the single-worker baseline this points straight at
     * the simulation's master profile (transport_ctx->profile), so the lone
     * worker writes the master directly and no merge is needed. */
    struct osh_transport_profile *profile;
};

/**
 * @brief Populate a worker context for the history range [@p hist_lo, @p hist_hi),
 *        borrowing the caller's pre-allocated pool and geometry scratch.
 *
 * @details
 * The single-worker baseline does not allocate: it points the context at the
 * simulation's pre-allocated ion pool and scratch (transport_ctx->ion_pool /
 * zone_refs / dist_batch), so nothing here is owned and there is nothing to free.
 * @ref capacity is taken from @p pool->capacity; the scratch arrays must be at
 * least that long.  The private scoring accumulators and the profile are left
 * NULL; the caller wires @ref profile to the run's master profile (serial) or to
 * a private one (parallel), and the accumulators stay NULL for shared-master
 * scoring (the current baseline).  A future parallel worker that owns a private
 * pool would use a separate owning constructor.
 *
 * @param[out] wctx       Context to populate (overwritten; must be non-NULL).
 * @param[in]  hist_lo    Inclusive lower bound of the assigned range.
 * @param[in]  hist_hi    Exclusive upper bound; must be > @p hist_lo.
 * @param[in]  pool       Borrowed live-history pool (must be non-NULL).
 * @param[in]  zone_refs  Borrowed zone-reference scratch (length >= pool->capacity).
 * @param[in]  dist_batch Borrowed boundary-distance scratch (length >= pool->capacity).
 */
void osh_worker_context_attach(struct osh_worker_context *wctx,
                               size_t hist_lo,
                               size_t hist_hi,
                               struct osh_particle_pool *pool,
                               struct osh_zone_ref *zone_refs,
                               double *dist_batch);

#ifdef __cplusplus
}
#endif

#endif /* OSH_WORKER_CONTEXT_H */
