#include "scoring/runtime/osh_scoring_shadow.h"

#include <stdlib.h>

#include "scoring/runtime/osh_scoring_postprocess.h"

enum osh_status osh_scoring_shadow_init(struct osh_scoring_shadow *shadow, struct osh_scoring_runtime const *live) {
    size_t i;

    if (!shadow || !live) {
        return OSH_EINVAL;
    }

    shadow->view = *live; /* alias outputs/geometries/settings/etc. */
    /* The view is a read-only presentation snapshot consumed only through
     * view.pages (osh_scoring_postprocess_into + the sink).  The struct-copy
     * above also carried master_acc (a shallow alias of the *live* per-page
     * accumulators) and master_scratch.crossing_buf (the *live* traversal
     * scratch), so a consumer reaching the view via
     * osh_scoring_runtime_master_accumulators() / _master_scratch() would
     * silently read live state, not the snapshot.  Drop those aliases so the
     * view is self-consistent; redirect happens only on the page array below. */
    shadow->view.master_acc = NULL;
    shadow->view.master_scratch.crossing_buf = NULL;
    shadow->view.master_scratch.crossing_cap = 0u;
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
        /* A mid-run dump is a data-only preview: the shadow copies only the `data`
         * array (osh_scoring_shadow_refresh), not the Welford M2 state, and the
         * live run keeps folding batches after the dump.  Drop the variance aliases
         * so a dump never (a) emits a standard-error column it cannot fill, nor
         * (b) lets a finalize on the view corrupt the still-live M2 arrays.  Final
         * (end-of-run) error bars come from the live runtime, not this shadow. */
        shadow->pages[i].variance = 0;
        shadow->pages[i].acc.data_var = NULL;
        shadow->pages[i].acc.data2_var = NULL;
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
                size_t n = lp->acc.len; /* never allocate a zero-length buffer */
                double *buf;
                if (n == 0u) {
                    n = 1u;
                }
                buf = (double *) malloc(n * sizeof(*buf));
                if (!buf) {
                    /* Roll back this attempt: free the scratch allocated so far and
                     * restore the aliased data pointers, leaving the shadow in its
                     * pristine all-aliased state.  scratch_ready stays 0, so a later
                     * refresh (or free) is well-defined and leaks nothing. */
                    size_t j;
                    for (j = 0; j < shadow->npages; ++j) {
                        if (shadow->scratch[j]) {
                            free(shadow->scratch[j]);
                            shadow->scratch[j] = NULL;
                            shadow->pages[j].acc.data = shadow->live->pages[j].acc.data;
                        }
                    }
                    return OSH_ENOMEM;
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
            uint64_t n = (uint64_t) lp->acc.len;
            /* osh_scoring_shadow_refresh() allocates one double even for a len==0
             * page (it never makes a zero-length buffer), so count one here too;
             * plain len would advertise 0 bytes and desync from the real scratch. */
            if (n == 0u) {
                n = 1u;
            }
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

    /* view.pages was set to shadow->pages (just freed); leaving it dangling
     * makes a post-free osh_scoring_shadow_view() hand back a runtime whose
     * page array points at freed memory.  Reset the view to a benign empty
     * snapshot so inspecting a freed shadow is well-defined, not a UAF. */
    shadow->view.pages = NULL;
    shadow->view.npages = 0u;
    shadow->view.master_acc = NULL;
    shadow->view.master_scratch.crossing_buf = NULL;
    shadow->view.master_scratch.crossing_cap = 0u;
}
