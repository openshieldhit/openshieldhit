#include "scoring/runtime/osh_scoring_postprocess.h"

#include <string.h>

#include "openshieldhit/const.h"
#include "scoring/runtime/osh_scoring_accumulator.h"

/* Guard for near-zero denominators in LET averaging. */
#define OSH_LET_DENOM_EPS 1.0e-300

int osh_scoring_postprocess_writes_data(enum osh_scoring_score_kind kind) {
    switch (kind) {
    case OSH_SCORING_SCORE_DOSEGY:
    case OSH_SCORING_SCORE_DLET:
    case OSH_SCORING_SCORE_TLET:
    case OSH_SCORING_SCORE_DQEFF:
    case OSH_SCORING_SCORE_TQEFF:
        return 1;
    default:
        return 0;
    }
}

/* Compute the presentation form of one page: read src, write dst->acc.data.
 * dst and src may be the same page (in-place) or dst->acc.data may alias / own a
 * separate buffer (out-of-place shadow). */
static enum osh_status page_postprocess_into(struct osh_scoring_page_runtime *dst,
                                             struct osh_scoring_page_runtime const *src) {
    size_t i;
    size_t len;

    if (!dst || !src) {
        return OSH_EINVAL;
    }
    len = src->acc.len;
    if (dst->acc.len != len) {
        return OSH_EINVAL;
    }

    switch (src->score_kind) {

    case OSH_SCORING_SCORE_DOSEGY:
        /* Convert accumulated MeV/g to Gy once per bin, not per transport step. */
        if (!dst->acc.data || !src->acc.data) {
            return OSH_EINVAL;
        }
        for (i = 0; i < len; ++i) {
            dst->acc.data[i] = src->acc.data[i] * OSH_MEVG2GY;
        }
        return OSH_OK;

    case OSH_SCORING_SCORE_DLET:
    case OSH_SCORING_SCORE_TLET:
    case OSH_SCORING_SCORE_DQEFF:
    case OSH_SCORING_SCORE_TQEFF:
        /* Finalise two-pass average: data = weighted sum, data2 = weight sum.
         * data2 is read as the divisor and never written. */
        if (!dst->acc.data || !src->acc.data || !src->acc.data2) {
            return OSH_EINVAL;
        }
        for (i = 0; i < len; ++i) {
            dst->acc.data[i] = (src->acc.data2[i] > OSH_LET_DENOM_EPS) ? (src->acc.data[i] / src->acc.data2[i]) : 0.0;
        }
        /* Clear flags on the dst page so save-layer validation passes; the src
         * page is left untouched, so live accumulation continues correctly. */
        dst->has_data2 = 0;
        dst->divide = 0;
        return OSH_OK;

    default:
        /* Guard against unhandled two-pass accumulators. */
        if (src->divide || src->has_data2) {
            return OSH_ENOTSUP;
        }

        /* Simple accumulators (ENERGY, FLUENCE, COUNT, …) are already in their
         * final per-step form.  When dst aliases src->acc.data (the in-place
         * wrapper, or a shadow's non-transformed page) there is nothing to do;
         * when dst owns a distinct buffer, mirror the raw values so the saved
         * view is complete.  The postproc mode (NORM/SUM/APPEND) describes how to
         * combine across runs — that belongs in the merge/save layer, not here. */
        switch (src->postproc) {
        case OSH_SCORING_POSTPROC_NONE:
        case OSH_SCORING_POSTPROC_SUM:
        case OSH_SCORING_POSTPROC_NORM:
        case OSH_SCORING_POSTPROC_APPEND:
            if (dst->acc.data != src->acc.data && dst->acc.data && src->acc.data) {
                memcpy(dst->acc.data, src->acc.data, len * sizeof(*dst->acc.data));
            }
            return OSH_OK;
        case OSH_SCORING_POSTPROC_AVER:
        default:
            return OSH_ENOTSUP;
        }
    }
}

enum osh_status osh_scoring_postprocess_into(struct osh_scoring_runtime *dst, struct osh_scoring_runtime const *src) {
    size_t i;
    enum osh_status rc;

    if (!dst || !src) {
        return OSH_EINVAL;
    }
    if (dst->npages != src->npages) {
        return OSH_EINVAL;
    }

    for (i = 0; i < src->npages; ++i) {
        rc = page_postprocess_into(&dst->pages[i], &src->pages[i]);
        if (rc != OSH_OK) {
            return rc;
        }
    }

    return OSH_OK;
}

enum osh_status osh_scoring_postprocess(struct osh_scoring_runtime *rt) {
    /* In place is the dst == src case: every page's dst.data aliases src.data,
     * so DOSEGY/LET write back over the live arrays exactly as before, and simple
     * scorers skip the (self-)copy. */
    return osh_scoring_postprocess_into(rt, rt);
}
