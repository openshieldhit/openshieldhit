#include "scoring/runtime/osh_scoring_postprocess.h"

#include <string.h>

#include "openshieldhit/const.h"
#include "scoring/runtime/osh_scoring_accumulator.h"

/* Guard for near-zero denominators in LET averaging. */
#define OSH_LET_DENOM_EPS 1.0e-300

int osh_scoring_postprocess_writes_data(enum osh_scoring_score_kind kind) {
    switch (kind) {
    case OSH_SCORING_SCORE_DOSE:
    case OSH_SCORING_SCORE_DOSEGY:
    case OSH_SCORING_SCORE_FLUENCE:
    case OSH_SCORING_SCORE_NKERMA:
    case OSH_SCORING_SCORE_DLET:
    case OSH_SCORING_SCORE_TLET:
    case OSH_SCORING_SCORE_DQEFF:
    case OSH_SCORING_SCORE_TQEFF:
        return 1;
    default:
        return 0;
    }
}

/* Volume-normalise an extensive tally in place or out-of-place:
 * dst[i] = src[i] * bin_vol_inv[i % diff_stride] * extra_factor.  The spatial bin
 * of flat index i is i % diff_stride (holds for diff1/diff2 layouts).  Used by the
 * DOSE/FLUENCE/NKERMA/DOSEGY estimators, which deposit the extensive quantity at
 * score time and divide by the bin volume exactly once here. */
static enum osh_status page_scale_by_bin_volume(struct osh_scoring_page_runtime *dst,
                                                struct osh_scoring_page_runtime const *src,
                                                struct osh_scoring_geometry_runtime const *geo,
                                                double extra_factor) {
    size_t i;
    size_t len;
    size_t stride;
    double const *vol;

    if (!dst->acc.data || !src->acc.data) {
        return OSH_EINVAL;
    }
    /* A compiled scored geometry always has bin_vol_inv; a hand-built runtime with
     * no geometry (e.g. a shadow/postprocess unit test) does not — then there is no
     * volume to divide by and only the extra factor applies. */
    vol = NULL;
    if (geo) {
        vol = geo->bin_vol_inv;
    }
    stride = src->diff_stride;
    if (vol && stride == 0u) {
        return OSH_EINVAL; /* real per-bin volume but no spatial stride — malformed */
    }
    len = src->acc.len;
    for (i = 0; i < len; ++i) {
        double vinv;
        if (vol) {
            vinv = vol[i % stride];
        } else {
            vinv = 1.0;
        }
        dst->acc.data[i] = src->acc.data[i] * vinv * extra_factor;
    }
    return OSH_OK;
}

/* Compute the presentation form of one page: read src, write dst->acc.data.
 * dst and src may be the same page (in-place) or dst->acc.data may alias / own a
 * separate buffer (out-of-place shadow). */
static enum osh_status page_postprocess_into(struct osh_scoring_page_runtime *dst,
                                             struct osh_scoring_page_runtime const *src,
                                             struct osh_scoring_geometry_runtime const *geo) {
    size_t i;
    size_t len;

    if (!dst || !src) {
        return OSH_EINVAL;
    }
    len = src->acc.len;
    if (dst->acc.len != len) {
        return OSH_EINVAL;
    }
    /* dst and src must describe the same page.  Only src->score_kind / postproc
     * are read for the transform, but a kind mismatch means the caller paired the
     * wrong pages, so fail fast rather than silently transform under the wrong
     * rule. */
    if (dst->score_kind != src->score_kind) {
        return OSH_EINVAL;
    }

    switch (src->score_kind) {

    case OSH_SCORING_SCORE_DOSE:
    case OSH_SCORING_SCORE_FLUENCE:
    case OSH_SCORING_SCORE_NKERMA:
        /* Divide each bin by its volume once: the scorer deposited the extensive
         * quantity (energy/rho, track length) and postprocess makes it intensive —
         * dose [MeV/g], fluence [1/cm^2], kerma [MeV/g]. */
        return page_scale_by_bin_volume(dst, src, geo, 1.0);

    case OSH_SCORING_SCORE_DOSEGY:
        /* As DOSE, then convert MeV/g -> Gy, once per bin. */
        return page_scale_by_bin_volume(dst, src, geo, OSH_MEVG2GY);

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
    /* A non-zero page count with a NULL page array is a malformed runtime; guard
     * it here so the per-page loop never dereferences NULL. The geometries array
     * may legitimately be absent in a hand-built runtime (unit tests). */
    if (src->npages > 0u && (!dst->pages || !src->pages)) {
        return OSH_EINVAL;
    }

    for (i = 0; i < src->npages; ++i) {
        struct osh_scoring_geometry_runtime const *geo;

        geo = NULL;
        if (src->geometries) {
            geo = &src->geometries[src->pages[i].geometry_idx];
        }
        rc = page_postprocess_into(&dst->pages[i], &src->pages[i], geo);
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
