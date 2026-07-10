#include "scoring/runtime/osh_scoring_postprocess.h"

#include <math.h>
#include <string.h>

#include "openshieldhit/const.h"
#include "scoring/runtime/osh_scoring_accumulator.h"
#include "scoring/runtime/osh_scoring_estimator.h"

/* Guard for near-zero denominators in LET averaging. */
#define OSH_LET_DENOM_EPS 1.0e-300

int osh_scoring_postprocess_writes_data(enum osh_scoring_score_kind kind) {
    switch (kind) {
    case OSH_SCORING_SCORE_DOSE:
    case OSH_SCORING_SCORE_DOSEGY:
    case OSH_SCORING_SCORE_DIRTYDOSE:
    case OSH_SCORING_SCORE_DIRTYDOSEGY:
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

static int page_needs_diff_normalization(struct osh_scoring_page_runtime const *page) {
    if (!page) {
        return 0;
    }
    if (page->diff_nbins == 0u && page->diff2_nbins == 0u) {
        return 0;
    }
    /* AVER pages already report ratios such as LET; APPEND pages are raw payloads.
     * Dividing either by a differential-axis width would change the physics/API. */
    return page->postproc != OSH_SCORING_POSTPROC_AVER && page->postproc != OSH_SCORING_POSTPROC_APPEND;
}

int osh_scoring_postprocess_page_writes_data(struct osh_scoring_page_runtime const *page) {
    if (!page) {
        return 0;
    }
    return osh_scoring_postprocess_writes_data(page->score_kind) || page_needs_diff_normalization(page);
}

static double diff_axis_width(double lo, double hi, size_t nbins, int is_log, size_t bin) {
    double t0;
    double t1;
    double ratio;

    if (nbins == 0u || bin >= nbins || !(hi > lo)) {
        return 0.0;
    }
    t0 = (double) bin / (double) nbins;
    t1 = (double) (bin + 1u) / (double) nbins;
    if (is_log) {
        if (!(lo > 0.0)) {
            return 0.0;
        }
        ratio = hi / lo;
        return lo * (pow(ratio, t1) - pow(ratio, t0));
    }
    return (hi - lo) / (double) nbins;
}

static double page_diff_normalization_factor(struct osh_scoring_page_runtime const *page, size_t data_idx) {
    double factor;

    if (!page_needs_diff_normalization(page)) {
        return 1.0;
    }
    factor = 1.0;
    if (page->diff_nbins > 0u) {
        size_t db;
        double width;
        if (page->diff_stride == 0u) {
            return 0.0;
        }
        db = (data_idx / page->diff_stride) % page->diff_nbins;
        width = diff_axis_width(page->diff_lo, page->diff_hi, page->diff_nbins, page->diff_log, db);
        if (!(width > 0.0)) {
            return 0.0;
        }
        factor /= width;
    }
    if (page->diff2_nbins > 0u) {
        size_t db;
        double width;
        if (page->diff2_stride == 0u) {
            return 0.0;
        }
        db = (data_idx / page->diff2_stride) % page->diff2_nbins;
        width = diff_axis_width(page->diff2_lo, page->diff2_hi, page->diff2_nbins, page->diff2_log, db);
        if (!(width > 0.0)) {
            return 0.0;
        }
        factor /= width;
    }
    return factor;
}

/* Volume-normalise an extensive tally in place or out-of-place:
 * dst[i] = src[i] * bin_vol_inv[i % diff_stride] * extra_factor * diff_factor.
 * The spatial bin of flat index i is i % diff_stride (holds for diff1/diff2
 * layouts).  Used by the DOSE/FLUENCE/NKERMA/DOSEGY estimators, which deposit
 * the extensive quantity at score time and divide by the bin volume exactly once
 * here.  Differential additive pages are additionally divided by differential
 * bin width(s), turning per-bin totals into differential quantities. */
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
        dst->acc.data[i] = src->acc.data[i] * vinv * extra_factor * page_diff_normalization_factor(src, i);
    }
    return OSH_OK;
}

static enum osh_status page_copy_with_diff_normalization(struct osh_scoring_page_runtime *dst,
                                                         struct osh_scoring_page_runtime const *src) {
    size_t i;
    size_t len;

    if (!dst->acc.data || !src->acc.data) {
        return OSH_EINVAL;
    }
    len = src->acc.len;
    for (i = 0u; i < len; ++i) {
        dst->acc.data[i] = src->acc.data[i] * page_diff_normalization_factor(src, i);
    }
    return OSH_OK;
}

/* ---- Per-estimator postprocess handlers (registered in osh_scoring_estimator.c) ---- */

enum osh_status postprocess_volume(struct osh_scoring_page_runtime *dst,
                                   struct osh_scoring_page_runtime const *src,
                                   struct osh_scoring_geometry_runtime const *geo) {
    /* DOSE/FLUENCE/NKERMA: divide each bin by its volume once, turning the extensive
     * deposit (energy/rho, track length) into the intensive quantity. */
    return page_scale_by_bin_volume(dst, src, geo, 1.0);
}

enum osh_status postprocess_dosegy(struct osh_scoring_page_runtime *dst,
                                   struct osh_scoring_page_runtime const *src,
                                   struct osh_scoring_geometry_runtime const *geo) {
    /* As DOSE, then convert MeV/g -> Gy once per bin. */
    return page_scale_by_bin_volume(dst, src, geo, OSH_MEVG2GY);
}

enum osh_status postprocess_ratio(struct osh_scoring_page_runtime *dst,
                                  struct osh_scoring_page_runtime const *src,
                                  struct osh_scoring_geometry_runtime const *geo) {
    /* DLET/TLET/DQEFF/TQEFF: finalise the two-pass average data/data2 (weighted sum
     * over weight sum).  Volume cancels in the ratio, so geo is unused. */
    size_t i;
    size_t len;

    (void) geo;
    if (!dst->acc.data || !src->acc.data || !src->acc.data2) {
        return OSH_EINVAL;
    }
    len = src->acc.len;
    for (i = 0; i < len; ++i) {
        if (src->acc.data2[i] > OSH_LET_DENOM_EPS) {
            dst->acc.data[i] = src->acc.data[i] / src->acc.data2[i];
        } else {
            dst->acc.data[i] = 0.0;
        }
    }
    /* Clear flags on dst so save-layer validation passes; src is left untouched so
     * live accumulation continues correctly. */
    dst->has_data2 = 0;
    dst->divide = 0;
    return OSH_OK;
}

/* Compute the presentation form of one page: read src, write dst->acc.data.
 * dst and src may be the same page (in-place) or dst->acc.data may alias / own a
 * separate buffer (out-of-place shadow). */
static enum osh_status page_postprocess_into(struct osh_scoring_page_runtime *dst,
                                             struct osh_scoring_page_runtime const *src,
                                             struct osh_scoring_geometry_runtime const *geo) {
    struct osh_scoring_estimator const *est;
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

    /* Dispatch to the estimator's postprocess handler.  A NULL handler means the
     * accumulator is already in its final per-step form (ENERGY, COUNT, …). */
    est = osh_scoring_estimator_for(src->score_kind);
    if (est && est->postprocess) {
        return est->postprocess(dst, src, geo);
    }

    /* No transform.  Guard against an unhandled two-pass accumulator, then mirror
     * the raw values when dst owns a distinct buffer.  The postproc mode
     * (NORM/SUM/APPEND) describes cross-run combination and belongs in the
     * merge/save layer, not here. */
    if (src->divide || src->has_data2) {
        return OSH_ENOTSUP;
    }
    switch (src->postproc) {
    case OSH_SCORING_POSTPROC_NONE:
    case OSH_SCORING_POSTPROC_SUM:
    case OSH_SCORING_POSTPROC_NORM:
    case OSH_SCORING_POSTPROC_APPEND:
        if (page_needs_diff_normalization(src)) {
            return page_copy_with_diff_normalization(dst, src);
        }
        if (dst->acc.data != src->acc.data && dst->acc.data && src->acc.data) {
            memcpy(dst->acc.data, src->acc.data, len * sizeof(*dst->acc.data));
        }
        return OSH_OK;
    case OSH_SCORING_POSTPROC_AVER:
    default:
        return OSH_ENOTSUP;
    }
}

enum osh_status osh_scoring_postprocess_into(struct osh_scoring_runtime *dst, struct osh_scoring_runtime const *src) {
    size_t i;
    enum osh_status rc;
    int in_place;

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

    /* dst == src is the in-place, destructive path: every page's dst.data aliases
     * src.data, so ÷volume / MeV·g→Gy / data÷data2 rewrite the live arrays and a
     * second pass would double-apply. The single-shot guard and its flag live here,
     * not just in the osh_scoring_postprocess() wrapper, so a direct in-place
     * osh_scoring_postprocess_into(rt, rt) is protected identically — the two
     * in-place entry points can no longer diverge. An out-of-place call (a distinct
     * dst shadow) never mutates src and stays repeatable for dump/checkpoint (#170)
     * and future merge. */
    in_place = ((void const *) dst == (void const *) src);
    if (in_place && dst->postprocessed) {
        return OSH_ESTATE;
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

    if (in_place) {
        dst->postprocessed = 1;
    }
    return OSH_OK;
}

enum osh_status osh_scoring_postprocess(struct osh_scoring_runtime *rt) {
    /* Thin wrapper over the dst == src case: osh_scoring_postprocess_into() owns
     * the single-shot guard and the postprocessed flag, so in-place finalisation
     * behaves identically whether reached through here or called directly. */
    return osh_scoring_postprocess_into(rt, rt);
}

/* Relative standard error of one additive estimator from its raw sum and M2.
 *   sum   = Σ x           (the raw accumulated sum, e.g. acc.data[i])
 *   m2    = Welford M2 across batches (acc.data_var[i])
 *   w     = Σ history count over all batches (acc.weight)
 *   dof   = B - 1 (> 0; the caller guards nbatch >= 2)
 * SE(mean) = sqrt(m2 / (dof·w)); relative = SE(mean) / |mean| with mean = sum/w,
 * which simplifies to sqrt(m2·w/dof)/|sum| — invariant under any linear rescale of
 * the reported value.  A zero sum (empty bin) has no defined relative error → 0. */
static double rel_error(double sum, double m2, double w, double dof) {
    double se;
    if (sum == 0.0 || m2 <= 0.0) {
        return 0.0;
    }
    se = sqrt(m2 * w / dof); /* = |mean| · relative, i.e. the SE of the raw sum */
    return se / fabs(sum);
}

enum osh_status osh_scoring_finalize_errors(struct osh_scoring_runtime *rt) {
    size_t p;
    size_t i;

    if (!rt) {
        return OSH_EINVAL;
    }
    /* Ordering guard: the relative error is derived from the raw sums (data, data2),
     * so this must run before an in-place postprocess rescales data / collapses the
     * two-pass ratio.  The postprocessed flag is set only by that destructive path
     * (dst==src); once set, the raw sums are gone and finalising would report error
     * against the rescaled value — reject rather than silently mis-report.  An
     * out-of-place postprocess into a distinct dst never sets the flag and leaves
     * src's raw sums intact, so finalising src stays valid.  Mirrors the single-shot
     * guard on the in-place postprocess path. */
    if (rt->postprocessed) {
        return OSH_ESTATE;
    }

    for (p = 0; p < rt->npages; ++p) {
        struct osh_scoring_accumulator *acc = &rt->pages[p].acc;
        double dof;

        /* No M2 array ⇒ variance off for this page (or a dump shadow, whose var
         * arrays are cleared): nothing to finalise. */
        if (!acc->data_var) {
            continue;
        }
        /* Fewer than two batches is zero degrees of freedom — no error is
         * defined.  Overwrite the raw M2 with zeros so the writer reports a clean
         * 0 rather than leaking sum-of-squares into the error column. */
        if (acc->nbatch < 2u || acc->weight <= 0.0) {
            memset(acc->data_var, 0, acc->len * sizeof(*acc->data_var));
            continue;
        }
        dof = (double) (acc->nbatch - 1u);

        if (acc->data2 && acc->data2_var) {
            /* Two-pass AVER quantity (DLET/TLET/Qeff): reported value is the ratio
             * data/data2.  Combine the numerator and denominator relative errors in
             * quadrature (covariance ignored ⇒ conservative for correlated num/den). */
            for (i = 0; i < acc->len; ++i) {
                double const rel_num = rel_error(acc->data[i], acc->data_var[i], acc->weight, dof);
                double const rel_den = rel_error(acc->data2[i], acc->data2_var[i], acc->weight, dof);
                acc->data_var[i] = sqrt(rel_num * rel_num + rel_den * rel_den);
            }
        } else {
            /* Additive quantity (DOSE/FLUENCE/ENERGY/COUNT): relative error of the
             * single accumulated sum. */
            for (i = 0; i < acc->len; ++i) {
                acc->data_var[i] = rel_error(acc->data[i], acc->data_var[i], acc->weight, dof);
            }
        }
    }

    return OSH_OK;
}
