#include "scoring/runtime/osh_scoring_postprocess.h"

static enum osh_status page_postprocess_status(struct osh_scoring_page_runtime const *page);

enum osh_status osh_scoring_postprocess(struct osh_scoring_runtime *rt) {
    size_t i;
    enum osh_status rc;

    if (!rt) {
        return OSH_EINVAL;
    }

    for (i = 0; i < rt->npages; ++i) {
        rc = page_postprocess_status(&rt->pages[i]);
        if (rc != OSH_OK) {
            return rc;
        }
    }

    return OSH_OK;
}

static enum osh_status page_postprocess_status(struct osh_scoring_page_runtime const *page) {
    if (!page) {
        return OSH_EINVAL;
    }

    /* TODO: implement actual postprocessing once save/runtime semantics are
     * settled. For now, simple accumulators are already in their final form.
     * Guard the cases that would be incorrect to silently ignore. */
    if (page->divide || page->has_data2) {
        return OSH_ENOTSUP;
    }

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
