#include "scoring/runtime/osh_scoring_estimator_common.h"
#include "scoring/runtime/osh_scoring_estimator_internal.h"

/**
 * @brief Point-deposit counterpart of osh_scoring_estimator_step_energy.
 *
 * Books the whole released energy (st->de) at the single located bin.  A point
 * has no track, so there is no path fraction.
 */
enum osh_status osh_scoring_estimator_point_energy(struct osh_scoring_runtime const *rt,
                                                   struct osh_scoring_accumulator *acc_set,
                                                   struct osh_scoring_geometry_score_group const *group,
                                                   size_t spatial_idx,
                                                   struct particle const *part,
                                                   struct step const *st) {
    size_t i;
    size_t db;           /* Index into the differential bins for the current spatial bin */
    size_t db2;          /* Index into the differential bins for the current spatial bin */
    size_t score_idx;    /* Index into the scoring array for the current spatial bin */
    double energy_score; /* Score for the current spatial bin */
    struct osh_scoring_page_runtime const *page;
    struct osh_scoring_accumulator *acc;

    for (i = 0; i < group->npages; ++i) {
        page = &rt->pages[group->first_page + i];
        acc = &acc_set[group->first_page + i];
        if (!osh_scoring_page_passes_filters(page, part, st)) {
            continue;
        }
        if (!osh_scoring_estimator_resolve_diff_bins(page, rt, part, st, &db, &db2)) {
            continue;
        }
        if (spatial_idx >= page->diff_stride) {
            return OSH_ESTATE;
        }
        /* Point ENERGY books the whole local energy release at one spatial bin.
         * There is no path fraction on the point-deposit path. */
        score_idx = osh_scoring_estimator_flat_bin(page, spatial_idx, db, db2);
        energy_score = st->de;
        osh_score_deposit(acc->data, score_idx, energy_score);
    }
    return OSH_OK;
}

/**
 * @brief Point-deposit counterpart of osh_scoring_estimator_step_dose.
 *
 * Books de/rho (times a dose-to-medium SP ratio) at the single located bin.
 * Neutral particles release energy here but do not deposit charged-particle dose
 * locally, so they book energy but not dose.
 */
enum osh_status osh_scoring_estimator_point_dose(struct osh_scoring_runtime const *rt,
                                                 struct osh_scoring_accumulator *acc_set,
                                                 struct osh_scoring_geometry_score_group const *group,
                                                 size_t spatial_idx,
                                                 struct particle const *part,
                                                 struct step const *st) {
    size_t i;
    size_t db;
    size_t db2;
    size_t score_idx;  /* Index into the scoring array for the current spatial bin */
    double sp_ratio;   /* Stopping-power ratio for converting dose-to-transport into dose-to-medium */
    double dose_score; /* Score for the current spatial bin */
    struct osh_scoring_dose_sp_ctx sp;
    struct osh_scoring_page_runtime const *page;
    struct osh_scoring_accumulator *acc;

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

        /* The stopping-power ratio is computed once per page, not per bin, because the
         * page's Settings override (if any) determines the medium for the SP
         * calculation.  The located bin may be in a different medium than the
         * transport medium, so the SP ratio is needed to convert dose-to-transport
         * into dose-to-medium. */
        sp_ratio = osh_scoring_estimator_dose_sp_ratio(&sp, rt, page);
        if (!osh_scoring_estimator_resolve_diff_bins(page, rt, part, st, &db, &db2)) {
            continue;
        }
        if (spatial_idx >= page->diff_stride) {
            return OSH_ESTATE;
        }
        /* Point DOSE books de/rho [MeV*cm3/g] at one spatial bin.  The optional
         * stopping-power ratio converts dose-to-transport into dose-to-medium;
         * postprocess divides by bin volume. */
        score_idx = osh_scoring_estimator_flat_bin(page, spatial_idx, db, db2);
        dose_score = (st->de / st->rho) * sp_ratio;
        osh_score_deposit(acc->data, score_idx, dose_score);
    }
    return OSH_OK;
}
