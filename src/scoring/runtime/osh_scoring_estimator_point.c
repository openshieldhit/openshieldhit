#include "scoring/runtime/osh_scoring_estimator_common.h"
#include "scoring/runtime/osh_scoring_estimator_internal.h"
#include "scoring/runtime/osh_scoring_kernels.h"

/**
 * @brief Point-deposit counterpart of score_step_energy.
 *
 * Books the whole released energy (st->de) at the single located bin.  A point
 * has no track, so there is no path fraction.
 */
enum osh_status score_point_energy(OSH_SCORING_DEPOSIT_PARAMS) {
    size_t i;
    size_t j;
    size_t db;
    size_t db2;
    struct osh_scoring_page_runtime const *page;
    struct osh_scoring_accumulator *acc;
    (void) score_len; /* a point deposit has no track length */

    for (i = 0; i < group->npages; ++i) {
        page = &rt->pages[group->first_page + i];
        acc = &acc_set[group->first_page + i];
        if (!osh_scoring_page_passes_filters(page, part, st)) {
            continue;
        }
        if (!osh_scoring_estimator_resolve_diff_bins(page, rt, part, st, &db, &db2)) {
            continue;
        }
        for (j = 0; j < ncross; ++j) {
            if (crossings[j].idx >= page->diff_stride) {
                return OSH_ESTATE;
            }
            osh_score_deposit(acc->data,
                              osh_scoring_estimator_flat_bin(page, crossings[j].idx, db, db2),
                              osh_kernel_point_energy(st->de));
        }
    }
    return OSH_OK;
}

/**
 * @brief Point-deposit counterpart of score_step_dose.
 *
 * Books de/rho (times a dose-to-medium SP ratio) at the single located bin.
 * Neutral particles release energy here but do not deposit charged-particle dose
 * locally, so they book energy but not dose.
 */
enum osh_status score_point_dose(OSH_SCORING_DEPOSIT_PARAMS) {
    size_t i;
    size_t j;
    size_t db;
    size_t db2;
    double sp_ratio;
    struct osh_scoring_dose_sp_ctx sp;
    struct osh_scoring_page_runtime const *page;
    struct osh_scoring_accumulator *acc;
    (void) score_len; /* a point deposit has no track length */

    if (part->charge == 0) {
        return OSH_OK;
    }
    if (!(st->rho > 0.0)) {
        return OSH_OK;
    }
    sp = osh_scoring_estimator_dose_sp_gather(rt, part, st);

    /* group is a run of DOSE or DOSEGY pages for the same geometry, not a run of
     * bins.  Each page may have different filters, Settings overrides, and
     * differential axes, so the page loop remains even though point scoring
     * usually has only one located crossing. */
    for (i = 0; i < group->npages; ++i) {
        page = &rt->pages[group->first_page + i];
        acc = &acc_set[group->first_page + i];
        if (!osh_scoring_page_passes_filters(page, part, st)) {
            continue;
        }
        sp_ratio = osh_scoring_estimator_dose_sp_ratio(&sp, rt, page);
        if (!osh_scoring_estimator_resolve_diff_bins(page, rt, part, st, &db, &db2)) {
            continue;
        }
        for (j = 0; j < ncross; ++j) {
            if (crossings[j].idx >= page->diff_stride) {
                return OSH_ESTATE;
            }
            osh_score_deposit(acc->data,
                              osh_scoring_estimator_flat_bin(page, crossings[j].idx, db, db2),
                              osh_kernel_point_dose(st->de, st->rho, sp_ratio));
        }
    }
    return OSH_OK;
}
