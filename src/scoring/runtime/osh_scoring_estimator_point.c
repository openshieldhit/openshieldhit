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

/**
 * @brief Point-deposit counterpart of osh_scoring_estimator_step_dirtydose.
 *
 * Identical to osh_scoring_estimator_point_dose(), except a page books dose only
 * when the projectile's mass stopping power in that page's scoring medium exceeds
 * OSH_DIRTYDOSE_MASS_SP_THRESHOLD [MeV*cm^2/g].  Neutrals are excluded up front
 * (charge == 0); the mass-SP gate additionally drops charged particles below the
 * threshold and any species without an SP-table entry.
 */
enum osh_status osh_scoring_estimator_point_dirtydose(struct osh_scoring_runtime const *rt,
                                                      struct osh_scoring_accumulator *acc_set,
                                                      struct osh_scoring_geometry_score_group const *group,
                                                      size_t spatial_idx,
                                                      struct particle const *part,
                                                      struct step const *st) {
    size_t i;
    size_t db;
    size_t db2;
    size_t score_idx;
    double sp_ratio;
    double mass_sp;
    double dose_score;
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

    for (i = 0; i < group->npages; ++i) {
        page = &rt->pages[group->first_page + i];
        acc = &acc_set[group->first_page + i];
        if (!osh_scoring_page_passes_filters(page, part, st)) {
            continue;
        }
        /* Mass-LET in this page's scoring medium [MeV*cm^2/g]; the override (if any)
         * is folded in exactly as it is for the dose scale. */
        sp_ratio = osh_scoring_estimator_dose_sp_ratio(&sp, rt, page);
        mass_sp = sp.sp_tr * sp_ratio;
        if (!(mass_sp > OSH_DIRTYDOSE_MASS_SP_THRESHOLD)) {
            continue; /* below the dirty-dose LET threshold (or mass SP unavailable) */
        }
        if (!osh_scoring_estimator_resolve_diff_bins(page, rt, part, st, &db, &db2)) {
            continue;
        }
        if (spatial_idx >= page->diff_stride) {
            return OSH_ESTATE;
        }
        score_idx = osh_scoring_estimator_flat_bin(page, spatial_idx, db, db2);
        dose_score = (st->de / st->rho) * sp_ratio;
        osh_score_deposit(acc->data, score_idx, dose_score);
    }
    return OSH_OK;
}

/**
 * @brief Point-deposit counterpart of osh_scoring_estimator_step_dlet.
 *
 * A sub-threshold recoil or fragment deposits its whole birth energy locally with
 * no track length (issue #227).  Its representative LET is the stopping power at
 * the birth energy S(medium, E_birth) * rho, read from the same SP table the step
 * scorer uses; because @c st->p[3] == st->q[3] on the point path, the shared LET
 * gather resolves that value at the birth energy (clamped to the table's
 * low-energy end for very slow recoils).  The dose weight is the whole local
 * energy release @c st->de, exactly like a track step's per-bin dose weight, so
 * DLET books (LET * de, de) into (data, data2) and the ratio postprocess yields
 * LET for an isolated deposit.
 *
 * A recoil species with no SP-table column has no representable dE/dx and is
 * skipped here; ENERGY and DOSE still score it on their own point paths.
 */
enum osh_status osh_scoring_estimator_point_dlet(struct osh_scoring_runtime const *rt,
                                                 struct osh_scoring_accumulator *acc_set,
                                                 struct osh_scoring_geometry_score_group const *group,
                                                 size_t spatial_idx,
                                                 struct particle const *part,
                                                 struct step const *st) {
    size_t i;
    double let_point;        /* LET this page scores for the deposit [MeV/cm] */
    double dose_weight;      /* whole local energy release [MeV] */
    double dlet_numerator;   /* LET * dose_weight, booked into acc->data */
    double dlet_denominator; /* dose_weight, booked into acc->data2 */
    struct osh_scoring_step_let_ctx let_ctx;
    struct osh_scoring_page_runtime const *page;
    struct osh_scoring_accumulator *acc;

    if (part->z == 0 || part->a == 0 || !(st->rho > 0.0)) {
        return OSH_OK;
    }
    /* Mean step energy equals the birth energy here (p[3] == q[3]).  The de/score_len
     * fallback baked into the gather is unused: we require an SP-table column
     * (have_proj) below, since a point deposit has no track to derive LET from. */
    let_ctx = osh_scoring_estimator_step_let_gather(rt, part, st, 1.0);
    if (!let_ctx.have_proj) {
        return OSH_OK;
    }

    for (i = 0; i < group->npages; ++i) {
        page = &rt->pages[group->first_page + i];
        acc = &acc_set[group->first_page + i];
        if (!osh_scoring_page_passes_filters(page, part, st)) {
            continue;
        }
        /* A per-page Settings override changes the LET value, honoured by _apply. */
        let_point = osh_scoring_estimator_step_let_apply(&let_ctx, rt, page);
        if (!(let_point > 0.0)) {
            continue;
        }
        if (spatial_idx >= page->diff_stride) {
            return OSH_ESTATE;
        }
        /* Dose-weighted LET sample.  Differential axes are not valid for DLET pages,
         * so spatial_idx is the complete bin index. */
        dose_weight = st->de;
        dlet_numerator = let_point * dose_weight;
        dlet_denominator = dose_weight;
        osh_score_deposit(acc->data, spatial_idx, dlet_numerator);
        osh_score_deposit(acc->data2, spatial_idx, dlet_denominator);
    }
    return OSH_OK;
}

/**
 * @brief Point-deposit counterpart of osh_scoring_estimator_step_tlet.
 *
 * Same representative birth-energy LET as osh_scoring_estimator_point_dlet, but
 * track-averaged.  A stopping recoil traverses an effective track length equal to
 * its residual range, approximated by de/LET (constant-LET slowing down).  This is
 * exactly "like a track step" with LET = de/ds, so TLET books
 * (LET * (de/LET), de/LET) == (de, de/LET) into (data, data2): the ratio
 * postprocess collapses to LET for an isolated deposit, while the de/LET
 * denominator gives the correct track-length weight against real track steps and
 * other recoils sharing the bin (a high-LET, short-range recoil counts less than a
 * low-LET, long-range one).
 */
enum osh_status osh_scoring_estimator_point_tlet(struct osh_scoring_runtime const *rt,
                                                 struct osh_scoring_accumulator *acc_set,
                                                 struct osh_scoring_geometry_score_group const *group,
                                                 size_t spatial_idx,
                                                 struct particle const *part,
                                                 struct step const *st) {
    size_t i;
    double let_point;        /* LET this page scores for the deposit [MeV/cm] */
    double track_weight;     /* effective track length de/LET [cm] */
    double tlet_numerator;   /* LET * track_weight, booked into acc->data */
    double tlet_denominator; /* track_weight, booked into acc->data2 */
    struct osh_scoring_step_let_ctx let_ctx;
    struct osh_scoring_page_runtime const *page;
    struct osh_scoring_accumulator *acc;

    if (part->z == 0 || part->a == 0 || !(st->rho > 0.0)) {
        return OSH_OK;
    }
    /* See osh_scoring_estimator_point_dlet: birth-energy LET, have_proj required. */
    let_ctx = osh_scoring_estimator_step_let_gather(rt, part, st, 1.0);
    if (!let_ctx.have_proj) {
        return OSH_OK;
    }

    for (i = 0; i < group->npages; ++i) {
        page = &rt->pages[group->first_page + i];
        acc = &acc_set[group->first_page + i];
        if (!osh_scoring_page_passes_filters(page, part, st)) {
            continue;
        }
        let_point = osh_scoring_estimator_step_let_apply(&let_ctx, rt, page);
        if (!(let_point > 0.0)) {
            continue; /* no representable LET -> also avoids the de/LET divide */
        }
        if (spatial_idx >= page->diff_stride) {
            return OSH_ESTATE;
        }
        /* Effective track length of a stopping recoil: de/LET (constant-LET range).
         * Differential axes are not valid for TLET pages, so spatial_idx is the
         * complete bin index. */
        track_weight = st->de / let_point;
        tlet_numerator = let_point * track_weight;
        tlet_denominator = track_weight;
        osh_score_deposit(acc->data, spatial_idx, tlet_numerator);
        osh_score_deposit(acc->data2, spatial_idx, tlet_denominator);
    }
    return OSH_OK;
}
