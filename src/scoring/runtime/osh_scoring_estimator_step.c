#include "material/runtime/osh_material_runtime.h"
#include "scoring/runtime/osh_scoring_estimator_common.h"
#include "scoring/runtime/osh_scoring_estimator_internal.h"

/**
 * @brief Accumulate energy deposition [MeV] into the ENERGY scorer pages.
 *
 * Distributes st->de proportionally to path length in each crossed voxel.
 */
enum osh_status osh_scoring_estimator_step_energy(struct osh_scoring_runtime const *rt,
                                                  struct osh_scoring_accumulator *acc_set,
                                                  struct osh_scoring_geometry_score_group const *group,
                                                  struct osh_voxel_crossing const *crossings,
                                                  size_t ncross,
                                                  struct particle const *part,
                                                  struct step const *st,
                                                  double score_len) {
    size_t i;
    size_t j;
    size_t db;           /* Index into the differential bins for the current spatial bin */
    size_t db2;          /* Index into the differential bins for the current spatial bin */
    size_t score_idx;    /* Index into the scoring array for the current spatial bin */
    double energy_score; /* Score for the current spatial bin */
    struct osh_scoring_page_runtime const *page;
    struct osh_scoring_accumulator *acc;

    /* A group is a run of ENERGY pages for one geometry.  Filters and
     * differential axes remain page-local, so resolve them before crossing I/O. */
    for (i = 0; i < group->npages; ++i) {
        page = &rt->pages[group->first_page + i];
        acc = &acc_set[group->first_page + i];
        if (!osh_scoring_page_passes_filters(page, part, st)) {
            continue;
        }
        if (!osh_scoring_estimator_resolve_diff_bins(page, rt, part, st, &db, &db2)) {
            continue;
        }
        /* The geometry driver already produced path lengths per spatial bin.
         * ENERGY distributes the step energy by each crossing's path fraction. */
        for (j = 0; j < ncross; ++j) {
            if (crossings[j].idx >= page->diff_stride) {
                return OSH_ESTATE;
            }
            /* ENERGY is extensive.  A transport step releases st->de [MeV];
             * each crossed bin receives the path fraction path_len/score_len. */
            score_idx = osh_scoring_estimator_flat_bin(page, crossings[j].idx, db, db2);
            energy_score = st->de * (crossings[j].path_len / score_len);
            osh_score_deposit(acc->data, score_idx, energy_score);
        }
    }
    return OSH_OK;
}

/**
 * @brief Accumulate fluence [1/cm2] into the FLUENCE scorer pages.
 *
 * Deposits the track length per crossing; the estimator postprocess divides by
 * the per-bin volume (geo->bin_vol_inv) once, yielding fluence [1/cm2].
 */
enum osh_status osh_scoring_estimator_step_fluence(struct osh_scoring_runtime const *rt,
                                                   struct osh_scoring_accumulator *acc_set,
                                                   struct osh_scoring_geometry_score_group const *group,
                                                   struct osh_voxel_crossing const *crossings,
                                                   size_t ncross,
                                                   struct particle const *part,
                                                   struct step const *st,
                                                   double score_len) {
    size_t i;
    size_t j;
    size_t db;            /* Index into the differential bins for the current spatial bin */
    size_t db2;           /* Index into the differential bins for the current spatial bin */
    size_t score_idx;     /* Index into the scoring array for the current spatial bin */
    double fluence_score; /* Score for the current spatial bin */
    struct osh_scoring_page_runtime const *page;
    struct osh_scoring_accumulator *acc;
    (void) score_len; /* fluence books raw track length; volume division is done in postprocess */

    /* FLUENCE pages share the same geometry traversal.  Only filters and
     * differential axes can change from page to page inside the group. */
    for (i = 0; i < group->npages; ++i) {
        page = &rt->pages[group->first_page + i];
        acc = &acc_set[group->first_page + i];
        if (!osh_scoring_page_passes_filters(page, part, st)) {
            continue;
        }
        if (!osh_scoring_estimator_resolve_diff_bins(page, rt, part, st, &db, &db2)) {
            continue;
        }
        /* Raw track length is accumulated here; postprocess applies the
         * geometry volume normalization. */
        for (j = 0; j < ncross; ++j) {
            if (crossings[j].idx >= page->diff_stride) {
                return OSH_ESTATE;
            }
            /* FLUENCE raw score is track length [cm].  Postprocess divides the
             * accumulated length by bin volume to obtain fluence [1/cm2]. */
            score_idx = osh_scoring_estimator_flat_bin(page, crossings[j].idx, db, db2);
            fluence_score = crossings[j].path_len;
            osh_score_deposit(acc->data, score_idx, fluence_score);
        }
    }
    return OSH_OK;
}

/**
 * @brief Accumulate dose [MeV/g] into the DOSE scorer pages.
 *
 * Deposits the extensive quantity path_len * de/(score_len*rho); the per-bin
 * volume division is applied once in postprocess (geo->bin_vol_inv), not here.
 *
 * When @c mat_tables is available and a Settings block specifies a medium,
 * applies a stopping-power ratio correction S(ovr,E)/S(tr,E) for
 * dose-to-medium scoring.  Pure density overrides do not change the dose
 * (Fano theorem).
 */
enum osh_status osh_scoring_estimator_step_dose(struct osh_scoring_runtime const *rt,
                                                struct osh_scoring_accumulator *acc_set,
                                                struct osh_scoring_geometry_score_group const *group,
                                                struct osh_voxel_crossing const *crossings,
                                                size_t ncross,
                                                struct particle const *part,
                                                struct step const *st,
                                                double score_len) {
    size_t i;
    size_t j;
    size_t db;        /* Index into the differential bins for the current spatial bin */
    size_t db2;       /* Index into the differential bins for the current spatial bin */
    size_t score_idx; /* Index into the scoring array for the current spatial bin */
    double base_scale;
    double dose_scale;
    double dose_score;
    struct osh_scoring_dose_sp_ctx sp;
    struct osh_scoring_page_runtime const *page;
    struct osh_scoring_accumulator *acc;

    if (!(st->rho > 0.0)) {
        return OSH_OK;
    }
    /* base_scale = de / (rho * score_len).  Volume division is deferred to
     * postprocess (accumulates in [MeV*cm3/g]; postprocess yields [MeV/g]). */
    base_scale = st->de / (score_len * st->rho);
    sp = osh_scoring_estimator_dose_sp_gather(rt, part, st);

    /* group is a run of DOSE or DOSEGY pages for this geometry, not a run of
     * bins.  Each page may have different filters, Settings overrides, and
     * differential axes, so dose_scale and diff bins are page-local. */
    for (i = 0; i < group->npages; ++i) {
        page = &rt->pages[group->first_page + i];
        acc = &acc_set[group->first_page + i];
        if (!osh_scoring_page_passes_filters(page, part, st)) {
            continue;
        }
        /* Dose-to-medium is page-local because a Settings block can request a
         * different scoring medium than the transport medium. */
        dose_scale = base_scale * osh_scoring_estimator_dose_sp_ratio(&sp, rt, page);
        if (!osh_scoring_estimator_resolve_diff_bins(page, rt, part, st, &db, &db2)) {
            continue;
        }
        /* Accumulate dose times volume; postprocess performs the bin-volume
         * division for Mesh/Cyl/Zone consistently. */
        for (j = 0; j < ncross; ++j) {
            if (crossings[j].idx >= page->diff_stride) {
                return OSH_ESTATE;
            }
            /* DOSE raw score is path_len * de/(score_len*rho) [MeV*cm3/g],
             * with any dose-to-medium stopping-power ratio folded into
             * dose_scale.  Postprocess divides by bin volume. */
            score_idx = osh_scoring_estimator_flat_bin(page, crossings[j].idx, db, db2);
            dose_score = crossings[j].path_len * dose_scale;
            osh_score_deposit(acc->data, score_idx, dose_score);
        }
    }
    return OSH_OK;
}

/**
 * @brief Accumulate LET-gated ("dirty") dose [MeV/g] into the DIRTYDOSE/DIRTYDOSEGY pages.
 *
 * Identical to osh_scoring_estimator_step_dose(), except a page books dose only
 * when the projectile's mass stopping power in that page's scoring medium exceeds
 * OSH_DIRTYDOSE_MASS_SP_THRESHOLD [MeV*cm^2/g].  The gate uses the same medium as
 * the dose: sp_tr for the transport medium, or sp_tr*ratio == S(override) under a
 * Settings medium override (dose-to-water then also gates on mass-LET in water).
 * A density-only override is Fano-invariant (ratio == 1), so the gate uses the
 * transport-medium mass SP.  When the mass SP is unavailable (neutral, no SP-table
 * entry, or medium < 0) sp_tr is 0 and nothing is scored — dirty dose is a
 * charged-particle, table-defined quantity.
 */
enum osh_status osh_scoring_estimator_step_dirtydose(struct osh_scoring_runtime const *rt,
                                                     struct osh_scoring_accumulator *acc_set,
                                                     struct osh_scoring_geometry_score_group const *group,
                                                     struct osh_voxel_crossing const *crossings,
                                                     size_t ncross,
                                                     struct particle const *part,
                                                     struct step const *st,
                                                     double score_len) {
    size_t i;
    size_t j;
    size_t db;        /* Index into the differential bins for the current spatial bin */
    size_t db2;       /* Index into the differential bins for the current spatial bin */
    size_t score_idx; /* Index into the scoring array for the current spatial bin */
    double base_scale;
    double sp_ratio;
    double mass_sp;
    double dose_scale;
    double dose_score;
    struct osh_scoring_dose_sp_ctx sp;
    struct osh_scoring_page_runtime const *page;
    struct osh_scoring_accumulator *acc;

    if (!(st->rho > 0.0)) {
        return OSH_OK;
    }
    base_scale = st->de / (score_len * st->rho);
    sp = osh_scoring_estimator_dose_sp_gather(rt, part, st);

    for (i = 0; i < group->npages; ++i) {
        page = &rt->pages[group->first_page + i];
        acc = &acc_set[group->first_page + i];
        if (!osh_scoring_page_passes_filters(page, part, st)) {
            continue;
        }
        /* Mass-LET in this page's scoring medium [MeV*cm^2/g]; the override (if any)
         * is folded in exactly as it is for the dose scale below. */
        sp_ratio = osh_scoring_estimator_dose_sp_ratio(&sp, rt, page);
        mass_sp = sp.sp_tr * sp_ratio;
        if (!(mass_sp > OSH_DIRTYDOSE_MASS_SP_THRESHOLD)) {
            continue; /* below the dirty-dose LET threshold (or mass SP unavailable) */
        }
        dose_scale = base_scale * sp_ratio;
        if (!osh_scoring_estimator_resolve_diff_bins(page, rt, part, st, &db, &db2)) {
            continue;
        }
        for (j = 0; j < ncross; ++j) {
            if (crossings[j].idx >= page->diff_stride) {
                return OSH_ESTATE;
            }
            score_idx = osh_scoring_estimator_flat_bin(page, crossings[j].idx, db, db2);
            dose_score = crossings[j].path_len * dose_scale;
            osh_score_deposit(acc->data, score_idx, dose_score);
        }
    }
    return OSH_OK;
}

/**
 * @brief Accumulate dose-averaged LET [MeV/cm] via a two-pass accumulator.
 *
 * Uses S(medium,E)*rho from the SP tables when available; falls back to de/score_len.
 * Per-page Settings overrides apply S(ovr,E)*rho_ovr (medium) or S(tr,E)*rho_ovr (density-only).
 *
 * Raw accumulator meaning before postprocess:
 *   acc->data  = sum(LET * dose_weight) per spatial bin
 *   acc->data2 = sum(dose_weight) per spatial bin
 *
 * postprocess_ratio() later returns acc->data / acc->data2, i.e. LETd.
 */
enum osh_status osh_scoring_estimator_step_dlet(struct osh_scoring_runtime const *rt,
                                                struct osh_scoring_accumulator *acc_set,
                                                struct osh_scoring_geometry_score_group const *group,
                                                struct osh_voxel_crossing const *crossings,
                                                size_t ncross,
                                                struct particle const *part,
                                                struct step const *st,
                                                double score_len) {
    size_t i;
    size_t j;
    double let_step;         /* LET this page scores for the step [MeV/cm] */
    double dose_weight;      /* energy-deposition weight for the crossed bin [MeV] */
    double dlet_numerator;   /* LET * dose_weight, booked into acc->data */
    double dlet_denominator; /* dose_weight, booked into acc->data2 */
    struct osh_scoring_step_let_ctx let_ctx;
    struct osh_scoring_page_runtime const *page;
    struct osh_scoring_accumulator *acc;

    if (part->z == 0 || part->a == 0 || !(st->rho > 0.0)) {
        return OSH_OK;
    }
    /* LET default (transport medium, table-based; de/score_len fallback) is a
     * per-step quantity; the per-page Settings override is applied inside the loop. */
    let_ctx = osh_scoring_estimator_step_let_gather(rt, part, st, score_len);

    for (i = 0; i < group->npages; ++i) {
        page = &rt->pages[group->first_page + i];
        acc = &acc_set[group->first_page + i];
        if (!osh_scoring_page_passes_filters(page, part, st)) {
            continue;
        }
        /* A per-page Settings override changes the LET value in the numerator, not
         * the dose weight: DLET is still weighted by where energy was deposited. */
        let_step = osh_scoring_estimator_step_let_apply(&let_ctx, rt, page);
        /* DLET stores a dose-weighted average in two explicit accumulators.
         *
         *   acc->data[crossing.idx]  += LET [MeV/cm] * dose_weight [MeV]
         *   acc->data2[crossing.idx] += dose_weight [MeV]
         *
         * dose_weight is the energy-deposition share assigned to this crossed bin,
         * st->de * crossing.path_len / score_len.  Differential axes are not valid
         * for DLET pages, so crossing.idx is the complete bin index. */
        for (j = 0; j < ncross; ++j) {
            if (crossings[j].idx >= page->diff_stride) {
                return OSH_ESTATE;
            }
            dose_weight = st->de * crossings[j].path_len / score_len;
            dlet_numerator = let_step * dose_weight;
            dlet_denominator = dose_weight;
            osh_score_deposit(acc->data, crossings[j].idx, dlet_numerator);
            osh_score_deposit(acc->data2, crossings[j].idx, dlet_denominator);
        }
    }
    return OSH_OK;
}

/**
 * @brief Accumulate track-averaged LET [MeV/cm] via a two-pass accumulator.
 *
 * Same table lookup as osh_scoring_estimator_step_dlet(); uses track-length ds_vox as the weight
 * rather than dose weight.
 *
 * Raw accumulator meaning before postprocess:
 *   acc->data  = sum(LET * track_weight) per spatial bin
 *   acc->data2 = sum(track_weight) per spatial bin
 *
 * postprocess_ratio() later returns acc->data / acc->data2, i.e. LETt.
 */
enum osh_status osh_scoring_estimator_step_tlet(struct osh_scoring_runtime const *rt,
                                                struct osh_scoring_accumulator *acc_set,
                                                struct osh_scoring_geometry_score_group const *group,
                                                struct osh_voxel_crossing const *crossings,
                                                size_t ncross,
                                                struct particle const *part,
                                                struct step const *st,
                                                double score_len) {
    size_t i;
    size_t j;
    double let_step;         /* LET this page scores for the step [MeV/cm] */
    double track_weight;     /* track-length weight for the crossed bin [cm] */
    double tlet_numerator;   /* LET * track_weight, booked into acc->data */
    double tlet_denominator; /* track_weight, booked into acc->data2 */
    struct osh_scoring_step_let_ctx let_ctx;
    struct osh_scoring_page_runtime const *page;
    struct osh_scoring_accumulator *acc;

    if (part->z == 0 || part->a == 0 || !(st->rho > 0.0)) {
        return OSH_OK;
    }
    /* Same transport-medium LET default as DLET; TLET differs only in weighting by
     * track length rather than energy deposition. */
    let_ctx = osh_scoring_estimator_step_let_gather(rt, part, st, score_len);

    for (i = 0; i < group->npages; ++i) {
        page = &rt->pages[group->first_page + i];
        acc = &acc_set[group->first_page + i];
        if (!osh_scoring_page_passes_filters(page, part, st)) {
            continue;
        }
        /* Settings can change which medium/density defines the scored LET; the
         * track-weight denominator stays geometric and override-independent. */
        let_step = osh_scoring_estimator_step_let_apply(&let_ctx, rt, page);
        /* TLET stores a track-length-weighted average in two explicit accumulators.
         *
         *   acc->data[crossing.idx]  += LET [MeV/cm] * track_weight [cm]
         *   acc->data2[crossing.idx] += track_weight [cm]
         *
         * track_weight is the physical step length assigned to this crossed bin,
         * st->ds * crossing.path_len / score_len.  Differential axes are not valid
         * for TLET pages, so crossing.idx is the complete bin index. */
        for (j = 0; j < ncross; ++j) {
            if (crossings[j].idx >= page->diff_stride) {
                return OSH_ESTATE;
            }
            track_weight = st->ds * crossings[j].path_len / score_len;
            tlet_numerator = let_step * track_weight;
            tlet_denominator = track_weight;
            osh_score_deposit(acc->data, crossings[j].idx, tlet_numerator);
            osh_score_deposit(acc->data2, crossings[j].idx, tlet_denominator);
        }
    }
    return OSH_OK;
}

/**
 * @brief Accumulate dose-averaged (z_eff/beta)^2 [dimensionless] via a two-pass accumulator.
 *
 * Requires mat_tables for the rest mass needed to compute beta.  Neutrals and
 * particles without SP table entries contribute nothing.
 *
 * Raw accumulator meaning before postprocess:
 *   acc->data  = sum(qeff * dose_weight) per spatial bin
 *   acc->data2 = sum(dose_weight) per spatial bin
 *
 * postprocess_ratio() later returns acc->data / acc->data2, i.e. DQEFF.
 *
 * @ref Kalholm et al. Medical Physics. 2023 Jan;50(1):651-9.
 *      https://doi.org/10.1002/mp.16029
 */
enum osh_status osh_scoring_estimator_step_dqeff(struct osh_scoring_runtime const *rt,
                                                 struct osh_scoring_accumulator *acc_set,
                                                 struct osh_scoring_geometry_score_group const *group,
                                                 struct osh_voxel_crossing const *crossings,
                                                 size_t ncross,
                                                 struct particle const *part,
                                                 struct step const *st,
                                                 double score_len) {
    size_t i;
    size_t j;
    double mean_energy;
    double beta;
    double qeff;
    double dose_weight;
    double dqeff_numerator;
    double dqeff_denominator;
    struct osh_scoring_page_runtime const *page;
    struct osh_scoring_accumulator *acc;

    if (part->z == 0) {
        return OSH_OK;
    }
    if (!(part->mass > 0.0)) {
        return OSH_OK;
    }

    /* Qeff depends on projectile charge and beta at the step midpoint.  It is a
     * single scalar for this step; only the averaging weight varies per crossing. */
    mean_energy = 0.5 * (st->p[3] + st->q[3]);
    beta = osh_scoring_estimator_particle_beta(mean_energy, part->mass);
    if (!(beta > 0.0)) {
        return OSH_OK;
    }
    qeff = osh_scoring_estimator_particle_qeff((int) part->z, beta);

    for (i = 0; i < group->npages; ++i) {
        page = &rt->pages[group->first_page + i];
        acc = &acc_set[group->first_page + i];
        if (!osh_scoring_page_passes_filters(page, part, st)) {
            continue;
        }
        /* DQEFF stores a dose-weighted qeff average in two explicit
         * accumulators.
         *
         * dose_weight is the energy-deposition share assigned to this crossed
         * spatial bin: st->de * crossing.path_len / score_len [MeV].
         *
         * Numerator:
         *   acc->data[crossing.idx]  += qeff * dose_weight [MeV]
         *
         * Denominator:
         *   acc->data2[crossing.idx] += dose_weight [MeV]
         *
         * Differential axes are not valid for DQEFF pages, so crossing.idx is
         * the complete bin index here. */
        for (j = 0; j < ncross; ++j) {
            if (crossings[j].idx >= page->diff_stride) {
                return OSH_ESTATE;
            }
            dose_weight = st->de * crossings[j].path_len / score_len;
            dqeff_numerator = qeff * dose_weight;
            dqeff_denominator = dose_weight;
            osh_score_deposit(acc->data, crossings[j].idx, dqeff_numerator);
            osh_score_deposit(acc->data2, crossings[j].idx, dqeff_denominator);
        }
    }
    return OSH_OK;
}

/**
 * @brief Accumulate track-averaged (z_eff/beta)^2 [dimensionless] via a two-pass accumulator.
 *
 * Same as osh_scoring_estimator_step_dqeff() but uses track-length ds_vox as the weight.
 *
 * Raw accumulator meaning before postprocess:
 *   acc->data  = sum(qeff * track_weight) per spatial bin
 *   acc->data2 = sum(track_weight) per spatial bin
 *
 * postprocess_ratio() later returns acc->data / acc->data2, i.e. TQEFF.
 *
 * @ref Kalholm et al. Medical Physics. 2023 Jan;50(1):651-9.
 *      https://doi.org/10.1002/mp.16029
 */
enum osh_status osh_scoring_estimator_step_tqeff(struct osh_scoring_runtime const *rt,
                                                 struct osh_scoring_accumulator *acc_set,
                                                 struct osh_scoring_geometry_score_group const *group,
                                                 struct osh_voxel_crossing const *crossings,
                                                 size_t ncross,
                                                 struct particle const *part,
                                                 struct step const *st,
                                                 double score_len) {
    size_t i;
    size_t j;
    double mean_energy;
    double beta; /* Particle relative speed at the step midpoint.  Needed to compute qeff. */
    double qeff;
    double track_weight;      /* Track-length weight for the current spatial bin */
    double tqeff_numerator;   /* Scoring numerator for the TQEFF calculation */
    double tqeff_denominator; /* Scoring denominator for the TQEFF calculation */
    struct osh_scoring_page_runtime const *page;
    struct osh_scoring_accumulator *acc;

    if (part->z == 0) {
        return OSH_OK;
    }
    if (!(part->mass > 0.0)) {
        return OSH_OK;
    }

    /* Qeff is computed once for the step midpoint.  TQEFF then averages that
     * value with geometric track-length weights in each crossed bin. */
    mean_energy = 0.5 * (st->p[3] + st->q[3]);
    beta = osh_scoring_estimator_particle_beta(mean_energy, part->mass);
    if (!(beta > 0.0)) {
        return OSH_OK;
    }
    qeff = osh_scoring_estimator_particle_qeff((int) part->z, beta);

    for (i = 0; i < group->npages; ++i) {
        page = &rt->pages[group->first_page + i];
        acc = &acc_set[group->first_page + i];
        if (!osh_scoring_page_passes_filters(page, part, st)) {
            continue;
        }
        /* TQEFF stores a track-length-weighted qeff average in two explicit
         * accumulators.
         *
         * track_weight is the physical step length assigned to this crossed
         * spatial bin: st->ds * crossing.path_len / score_len [cm].
         *
         * Numerator:
         *   acc->data[crossing.idx]  += qeff * track_weight [cm]
         *
         * Denominator:
         *   acc->data2[crossing.idx] += track_weight [cm]
         *
         * Differential axes are not valid for TQEFF pages, so crossing.idx is
         * the complete bin index here. */
        for (j = 0; j < ncross; ++j) {
            if (crossings[j].idx >= page->diff_stride) {
                return OSH_ESTATE;
            }
            track_weight = st->ds * crossings[j].path_len / score_len;
            tqeff_numerator = qeff * track_weight;
            tqeff_denominator = track_weight;
            osh_score_deposit(acc->data, crossings[j].idx, tqeff_numerator);
            osh_score_deposit(acc->data2, crossings[j].idx, tqeff_denominator);
        }
    }
    return OSH_OK;
}

/**
 * @brief Accumulate dose-averaged kinetic energy [MeV] via a two-pass accumulator.
 *
 * Averages the step-midpoint kinetic energy 0.5*(p[3]+q[3]) [MeV] — the same
 * scalar already computed for DQEFF — weighted by the per-crossing
 * energy-deposition share.  Unlike LET/QEFF, kinetic energy is well-defined
 * for neutrals too, so this handler applies no charge/mass gate.
 *
 * Raw accumulator meaning before postprocess:
 *   acc->data  = sum(mean_energy * dose_weight) per spatial bin
 *   acc->data2 = sum(dose_weight) per spatial bin
 *
 * postprocess_ratio() later returns acc->data / acc->data2, i.e. DAVGE.
 */
enum osh_status osh_scoring_estimator_step_davge(struct osh_scoring_runtime const *rt,
                                                 struct osh_scoring_accumulator *acc_set,
                                                 struct osh_scoring_geometry_score_group const *group,
                                                 struct osh_voxel_crossing const *crossings,
                                                 size_t ncross,
                                                 struct particle const *part,
                                                 struct step const *st,
                                                 double score_len) {
    size_t i;
    size_t j;
    double mean_energy;       /* kinetic energy at the step midpoint [MeV] */
    double dose_weight;       /* energy-deposition weight for the crossed bin [MeV] */
    double davge_numerator;   /* mean_energy * dose_weight, booked into acc->data */
    double davge_denominator; /* dose_weight, booked into acc->data2 */
    struct osh_scoring_page_runtime const *page;
    struct osh_scoring_accumulator *acc;

    /* Kinetic energy is a single scalar for this step; only the averaging
     * weight varies per crossing. */
    mean_energy = 0.5 * (st->p[3] + st->q[3]);

    for (i = 0; i < group->npages; ++i) {
        page = &rt->pages[group->first_page + i];
        acc = &acc_set[group->first_page + i];
        if (!osh_scoring_page_passes_filters(page, part, st)) {
            continue;
        }
        /* DAVGE stores a dose-weighted average in two explicit accumulators.
         * Differential axes are not valid for DAVGE pages, so crossing.idx is
         * the complete bin index. */
        for (j = 0; j < ncross; ++j) {
            if (crossings[j].idx >= page->diff_stride) {
                return OSH_ESTATE;
            }
            dose_weight = st->de * crossings[j].path_len / score_len;
            davge_numerator = mean_energy * dose_weight;
            davge_denominator = dose_weight;
            osh_score_deposit(acc->data, crossings[j].idx, davge_numerator);
            osh_score_deposit(acc->data2, crossings[j].idx, davge_denominator);
        }
    }
    return OSH_OK;
}

/**
 * @brief Accumulate track-averaged kinetic energy [MeV] via a two-pass accumulator.
 *
 * Same as osh_scoring_estimator_step_davge() but uses track-length ds_vox as
 * the weight rather than dose weight.
 *
 * Raw accumulator meaning before postprocess:
 *   acc->data  = sum(mean_energy * track_weight) per spatial bin
 *   acc->data2 = sum(track_weight) per spatial bin
 *
 * postprocess_ratio() later returns acc->data / acc->data2, i.e. TAVGE.
 */
enum osh_status osh_scoring_estimator_step_tavge(struct osh_scoring_runtime const *rt,
                                                 struct osh_scoring_accumulator *acc_set,
                                                 struct osh_scoring_geometry_score_group const *group,
                                                 struct osh_voxel_crossing const *crossings,
                                                 size_t ncross,
                                                 struct particle const *part,
                                                 struct step const *st,
                                                 double score_len) {
    size_t i;
    size_t j;
    double mean_energy;       /* kinetic energy at the step midpoint [MeV] */
    double track_weight;      /* track-length weight for the crossed bin [cm] */
    double tavge_numerator;   /* mean_energy * track_weight, booked into acc->data */
    double tavge_denominator; /* track_weight, booked into acc->data2 */
    struct osh_scoring_page_runtime const *page;
    struct osh_scoring_accumulator *acc;

    mean_energy = 0.5 * (st->p[3] + st->q[3]);

    for (i = 0; i < group->npages; ++i) {
        page = &rt->pages[group->first_page + i];
        acc = &acc_set[group->first_page + i];
        if (!osh_scoring_page_passes_filters(page, part, st)) {
            continue;
        }
        for (j = 0; j < ncross; ++j) {
            if (crossings[j].idx >= page->diff_stride) {
                return OSH_ESTATE;
            }
            track_weight = st->ds * crossings[j].path_len / score_len;
            tavge_numerator = mean_energy * track_weight;
            tavge_denominator = track_weight;
            osh_score_deposit(acc->data, crossings[j].idx, tavge_numerator);
            osh_score_deposit(acc->data2, crossings[j].idx, tavge_denominator);
        }
    }
    return OSH_OK;
}

/**
 * @brief Accumulate dose-averaged relative speed beta = v/c [dimensionless]
 *        via a two-pass accumulator.
 *
 * beta is derived from the same step-midpoint kinetic energy as DAVGE via
 * osh_scoring_estimator_particle_beta(), which is exact for massless
 * particles too (part->mass == 0 for photons), so no charge/mass gate is
 * needed.  The isnan-style guard below only rejects the degenerate case where
 * both kinetic energy and rest mass are zero (0/0).
 *
 * Raw accumulator meaning before postprocess:
 *   acc->data  = sum(beta * dose_weight) per spatial bin
 *   acc->data2 = sum(dose_weight) per spatial bin
 *
 * postprocess_ratio() later returns acc->data / acc->data2, i.e. DBETA.
 */
enum osh_status osh_scoring_estimator_step_dbeta(struct osh_scoring_runtime const *rt,
                                                 struct osh_scoring_accumulator *acc_set,
                                                 struct osh_scoring_geometry_score_group const *group,
                                                 struct osh_voxel_crossing const *crossings,
                                                 size_t ncross,
                                                 struct particle const *part,
                                                 struct step const *st,
                                                 double score_len) {
    size_t i;
    size_t j;
    double mean_energy;
    double beta;
    double dose_weight;
    double dbeta_numerator;
    double dbeta_denominator;
    struct osh_scoring_page_runtime const *page;
    struct osh_scoring_accumulator *acc;

    mean_energy = 0.5 * (st->p[3] + st->q[3]);
    beta = osh_scoring_estimator_particle_beta(mean_energy, part->mass);
    if (!(beta >= 0.0)) {
        return OSH_OK;
    }

    for (i = 0; i < group->npages; ++i) {
        page = &rt->pages[group->first_page + i];
        acc = &acc_set[group->first_page + i];
        if (!osh_scoring_page_passes_filters(page, part, st)) {
            continue;
        }
        for (j = 0; j < ncross; ++j) {
            if (crossings[j].idx >= page->diff_stride) {
                return OSH_ESTATE;
            }
            dose_weight = st->de * crossings[j].path_len / score_len;
            dbeta_numerator = beta * dose_weight;
            dbeta_denominator = dose_weight;
            osh_score_deposit(acc->data, crossings[j].idx, dbeta_numerator);
            osh_score_deposit(acc->data2, crossings[j].idx, dbeta_denominator);
        }
    }
    return OSH_OK;
}

/**
 * @brief Accumulate track-averaged relative speed beta = v/c [dimensionless]
 *        via a two-pass accumulator.
 *
 * Same as osh_scoring_estimator_step_dbeta() but uses track-length ds_vox as
 * the weight rather than dose weight.
 *
 * Raw accumulator meaning before postprocess:
 *   acc->data  = sum(beta * track_weight) per spatial bin
 *   acc->data2 = sum(track_weight) per spatial bin
 *
 * postprocess_ratio() later returns acc->data / acc->data2, i.e. TBETA.
 */
enum osh_status osh_scoring_estimator_step_tbeta(struct osh_scoring_runtime const *rt,
                                                 struct osh_scoring_accumulator *acc_set,
                                                 struct osh_scoring_geometry_score_group const *group,
                                                 struct osh_voxel_crossing const *crossings,
                                                 size_t ncross,
                                                 struct particle const *part,
                                                 struct step const *st,
                                                 double score_len) {
    size_t i;
    size_t j;
    double mean_energy;
    double beta;
    double track_weight;
    double tbeta_numerator;
    double tbeta_denominator;
    struct osh_scoring_page_runtime const *page;
    struct osh_scoring_accumulator *acc;

    mean_energy = 0.5 * (st->p[3] + st->q[3]);
    beta = osh_scoring_estimator_particle_beta(mean_energy, part->mass);
    if (!(beta >= 0.0)) {
        return OSH_OK;
    }

    for (i = 0; i < group->npages; ++i) {
        page = &rt->pages[group->first_page + i];
        acc = &acc_set[group->first_page + i];
        if (!osh_scoring_page_passes_filters(page, part, st)) {
            continue;
        }
        for (j = 0; j < ncross; ++j) {
            if (crossings[j].idx >= page->diff_stride) {
                return OSH_ESTATE;
            }
            track_weight = st->ds * crossings[j].path_len / score_len;
            tbeta_numerator = beta * track_weight;
            tbeta_denominator = track_weight;
            osh_score_deposit(acc->data, crossings[j].idx, tbeta_numerator);
            osh_score_deposit(acc->data2, crossings[j].idx, tbeta_denominator);
        }
    }
    return OSH_OK;
}
