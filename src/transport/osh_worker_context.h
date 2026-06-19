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

/**
 * @brief Everything one transport worker needs to run an assigned slice of histories.
 *
 * @details
 * Bundles the per-worker *mutable* state so that a run can be partitioned across
 * workers (threads, or MPI ranks, or sequential profiling replicas) simply by
 * constructing several of these over disjoint history ranges and giving each its
 * own pool and scratch.  Nothing here is shared between workers, which is the
 * precondition for lock-free parallel transport.
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

    /* Per-worker live-history pool and geometry-query scratch, all sized to
     * @ref capacity and never resized on the hot path. */
    size_t capacity;                /* Pool capacity (a pure performance knob). */
    struct osh_particle_pool *pool; /* Live primaries/secondaries for this worker. */
    struct osh_zone_ref *zone_refs; /* Batch zone-reference scratch, length @ref capacity. */
    double *dist_batch;             /* Batch boundary-distance scratch, length @ref capacity. */

    /* Future per-worker private scoring memory.  When non-NULL the worker will
     * deposit into this private set and the driver will fold it into the shared
     * master with osh_scoring_accumulator_merge() once the range completes.  NULL
     * today: the single worker deposits straight into the shared master, so there
     * is nothing to merge.  Profiling per-worker vs shared scoring memory plugs in
     * here without touching the transport loop. */
    struct osh_scoring_accumulator *accumulators;
    size_t naccumulators;
};

/**
 * @brief Allocate a worker context for the history range [@p hist_lo, @p hist_hi).
 *
 * @details
 * Allocates the pool and batch scratch.  The pool capacity is
 * min(@p requested_capacity, range size) — never more slots than histories —
 * with @p requested_capacity of 0 meaning "use the range size".  The private
 * scoring accumulators are left NULL (shared-master scoring, the current
 * baseline).
 *
 * @param[out] wctx               Context to populate (overwritten; must be non-NULL).
 * @param[in]  hist_lo            Inclusive lower bound of the assigned range.
 * @param[in]  hist_hi            Exclusive upper bound; must be > @p hist_lo.
 * @param[in]  requested_capacity Desired pool capacity, or 0 for the range size.
 * @returns OSH_OK on success, OSH_EINVAL on a bad range/NULL, OSH_ENOMEM on
 *          allocation failure (nothing leaks; @p wctx is left freed).
 */
enum osh_status
osh_worker_context_init(struct osh_worker_context *wctx, size_t hist_lo, size_t hist_hi, size_t requested_capacity);

/**
 * @brief Release everything owned by a worker context.  Safe on NULL/zeroed.
 */
void osh_worker_context_free(struct osh_worker_context *wctx);

#ifdef __cplusplus
}
#endif

#endif /* OSH_WORKER_CONTEXT_H */
