#include "scoring/runtime/osh_scoring_shadow.h"

#include <stdlib.h>

#include "scoring/runtime/osh_scoring_postprocess.h"

enum osh_status osh_scoring_shadow_init(struct osh_scoring_shadow *shadow, struct osh_scoring_runtime const *live) {
    size_t i;

    if (!shadow || !live) {
        return OSH_EINVAL;
    }

    shadow->view = *live; /* alias outputs/geometries/settings/etc. */
    shadow->live = live;
    shadow->npages = live->npages;
    shadow->pages = NULL;
    shadow->scratch = NULL;
    shadow->scratch_ready = 0;

    if (live->npages == 0u) {
        shadow->view.pages = NULL;
        return OSH_OK;
    }

    shadow->pages = (struct osh_scoring_page_runtime *) calloc(live->npages, sizeof(*shadow->pages));
    if (!shadow->pages) {
        return OSH_ENOMEM;
    }
    shadow->scratch = (double **) calloc(live->npages, sizeof(*shadow->scratch));
    if (!shadow->scratch) {
        free(shadow->pages);
        shadow->pages = NULL;
        return OSH_ENOMEM;
    }

    for (i = 0; i < live->npages; ++i) {
        shadow->pages[i] = live->pages[i]; /* struct copy: aliases data/data2/var + scalars */
    }
    shadow->view.pages = shadow->pages;

    return OSH_OK;
}

enum osh_status osh_scoring_shadow_refresh(struct osh_scoring_shadow *shadow) {
    if (!shadow || !shadow->live) {
        return OSH_EINVAL;
    }

    /* Lazy, allocate-once: the first refresh grabs a private data buffer for each
     * page whose postprocess writes data and redirects the view page's data
     * pointer at it.  Non-transformed pages keep aliasing the live data. */
    if (!shadow->scratch_ready && shadow->npages > 0u) {
        size_t i;
        for (i = 0; i < shadow->npages; ++i) {
            struct osh_scoring_page_runtime const *lp = &shadow->live->pages[i];
            if (osh_scoring_postprocess_writes_data(lp->score_kind)) {
                size_t const n = lp->acc.len ? lp->acc.len : 1u;
                double *buf = (double *) malloc(n * sizeof(*buf));
                if (!buf) {
                    return OSH_ENOMEM; /* partially-allocated scratch is released by _free */
                }
                shadow->scratch[i] = buf;
                shadow->pages[i].acc.data = buf; /* redirect the write target off the live array */
            }
        }
        shadow->scratch_ready = 1;
    }

    return osh_scoring_postprocess_into(&shadow->view, shadow->live);
}

uint64_t osh_scoring_shadow_bytes(struct osh_scoring_shadow const *shadow) {
    uint64_t total = 0u;
    size_t i;

    if (!shadow || !shadow->live) {
        return 0u;
    }
    for (i = 0; i < shadow->npages; ++i) {
        struct osh_scoring_page_runtime const *lp = &shadow->live->pages[i];
        if (osh_scoring_postprocess_writes_data(lp->score_kind)) {
            uint64_t const n = (uint64_t) lp->acc.len;
            uint64_t const b =
                (n <= UINT64_MAX / (uint64_t) sizeof(double)) ? (n * (uint64_t) sizeof(double)) : UINT64_MAX;
            total = (total <= UINT64_MAX - b) ? (total + b) : UINT64_MAX;
        }
    }
    return total;
}

struct osh_scoring_runtime const *osh_scoring_shadow_view(struct osh_scoring_shadow const *shadow) {
    return shadow ? &shadow->view : NULL;
}

void osh_scoring_shadow_free(struct osh_scoring_shadow *shadow) {
    size_t i;

    if (!shadow) {
        return;
    }
    if (shadow->scratch) {
        for (i = 0; i < shadow->npages; ++i) {
            free(shadow->scratch[i]); /* only buffers we own; NULL for aliased pages */
        }
        free((void *) shadow->scratch); /* double** -> void*: explicit cast keeps clang-tidy happy */
    }
    free(shadow->pages);

    shadow->pages = NULL;
    shadow->scratch = NULL;
    shadow->npages = 0u;
    shadow->live = NULL;
    shadow->scratch_ready = 0;
}
