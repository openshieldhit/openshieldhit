#include "scoring/runtime/osh_scoring_postprocess.h"

#include "openshieldhit/const.h"

/* Guard for near-zero denominators in LET averaging. */
#define OSH_LET_DENOM_EPS 1.0e-300

static enum osh_status page_postprocess(struct osh_scoring_page_runtime *page);

enum osh_status osh_scoring_postprocess(struct osh_scoring_runtime *rt) {
    size_t i;
    enum osh_status rc;

    if (!rt) {
        return OSH_EINVAL;
    }

    for (i = 0; i < rt->npages; ++i) {
        rc = page_postprocess(&rt->pages[i]);
        if (rc != OSH_OK) {
            return rc;
        }
    }

    return OSH_OK;
}

static enum osh_status page_postprocess(struct osh_scoring_page_runtime *page) {
    size_t i;

    if (!page) {
        return OSH_EINVAL;
    }

    switch (page->score_kind) {

    case OSH_SCORING_SCORE_DOSE:
        /* Convert accumulated MeV/g to Gy once per bin, not per transport step. */
        for (i = 0; i < page->len; ++i) {
            page->data[i] *= OSH_MEVG2GY;
        }
        return OSH_OK;

    case OSH_SCORING_SCORE_DLET:
    case OSH_SCORING_SCORE_TLET:
    case OSH_SCORING_SCORE_DQEFF:
    case OSH_SCORING_SCORE_TQEFF:
        /* Finalise two-pass average: data = weighted_sum, data2 = weight_sum. */
        for (i = 0; i < page->len; ++i) {
            if (page->data2[i] > OSH_LET_DENOM_EPS) {
                page->data[i] /= page->data2[i];
            } else {
                page->data[i] = 0.0;
            }
        }
        /* Clear flags so save-layer validation passes. */
        page->has_data2 = 0;
        page->divide = 0;
        return OSH_OK;

    default:
        /* Guard against unhandled two-pass accumulators. */
        if (page->divide || page->has_data2) {
            return OSH_ENOTSUP;
        }

        /* Simple accumulators (ENERGY, FLUENCE, COUNT, …) are already in
         * their final per-step form.  The postproc mode (NORM, SUM, APPEND, …)
         * describes how to combine results across multiple simulation runs —
         * that logic belongs in the merge/save layer, not here.
         *
         * TODO: when multi-run merging is implemented (native multithreaded or
         * embarrassingly-parallel user-managed runs), this switch will perform
         * the weighted combination.  Each partial BDO file carries the total
         * primary count (nstat) so the merge layer knows the correct weight for
         * NORM vs. AVER vs. SUM, etc. */
        switch (page->postproc) {
        case OSH_SCORING_POSTPROC_NONE:
        case OSH_SCORING_POSTPROC_SUM:
        case OSH_SCORING_POSTPROC_NORM:
        case OSH_SCORING_POSTPROC_APPEND:
            return OSH_OK;
        case OSH_SCORING_POSTPROC_AVER:
        default:
            return OSH_ENOTSUP;
        }
    }
}
