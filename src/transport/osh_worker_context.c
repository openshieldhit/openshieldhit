#include "transport/osh_worker_context.h"

#include <stdlib.h>

#include "common/osh_particle_pool.h"
#include "openshieldhit/geometry.h"

enum osh_status
osh_worker_context_init(struct osh_worker_context *wctx, size_t hist_lo, size_t hist_hi, size_t requested_capacity) {
    size_t range;
    size_t capacity;
    enum osh_status rc;

    if (!wctx || hist_hi <= hist_lo) {
        return OSH_EINVAL;
    }

    wctx->hist_lo = hist_lo;
    wctx->hist_hi = hist_hi;
    wctx->capacity = 0u;
    wctx->pool = NULL;
    wctx->zone_refs = NULL;
    wctx->dist_batch = NULL;
    wctx->accumulators = NULL;
    wctx->naccumulators = 0u;

    range = hist_hi - hist_lo;
    capacity = (requested_capacity != 0u) ? requested_capacity : range;
    if (capacity > range) {
        capacity = range; /* never allocate more slots than histories */
    }

    /* Guard the size product against overflow now that capacity is caller-driven. */
    if (capacity > ((size_t) -1) / sizeof(*wctx->zone_refs)) {
        return OSH_ENOMEM;
    }

    rc = osh_particle_pool_alloc(capacity, &wctx->pool);
    if (rc != OSH_OK) {
        return rc;
    }
    wctx->zone_refs = (struct osh_zone_ref *) malloc(capacity * sizeof(*wctx->zone_refs));
    wctx->dist_batch = (double *) malloc(capacity * sizeof(*wctx->dist_batch));
    if (!wctx->zone_refs || !wctx->dist_batch) {
        osh_worker_context_free(wctx);
        return OSH_ENOMEM;
    }

    wctx->capacity = capacity;
    return OSH_OK;
}

void osh_worker_context_free(struct osh_worker_context *wctx) {
    if (!wctx) {
        return;
    }
    free(wctx->zone_refs);
    free(wctx->dist_batch);
    osh_particle_pool_free(wctx->pool);
    wctx->zone_refs = NULL;
    wctx->dist_batch = NULL;
    wctx->pool = NULL;
    wctx->capacity = 0u;
    /* accumulators are not owned here yet (always NULL in the single-worker
     * baseline); when per-worker scoring memory lands its owner frees it. */
}
