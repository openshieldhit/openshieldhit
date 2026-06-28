#include "transport/osh_worker_context.h"

#include "common/osh_particle_pool.h"

void osh_worker_context_attach(struct osh_worker_context *wctx,
                               size_t hist_lo,
                               size_t hist_hi,
                               struct osh_particle_pool *pool,
                               struct osh_zone_ref *zone_refs,
                               double *dist_batch) {
    if (!wctx) {
        return;
    }
    wctx->hist_lo = hist_lo;
    wctx->hist_hi = hist_hi;
    wctx->capacity = pool ? pool->capacity : 0u;
    wctx->pool = pool;             /* borrowed: owned by the simulation */
    wctx->zone_refs = zone_refs;   /* borrowed */
    wctx->dist_batch = dist_batch; /* borrowed */
    wctx->accumulators = NULL;
    wctx->naccumulators = 0u;
    wctx->profile = NULL; /* caller wires this to the master (serial) or a private profile (parallel). */
}
