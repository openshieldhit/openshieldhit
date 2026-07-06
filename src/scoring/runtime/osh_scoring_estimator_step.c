#include "material/runtime/osh_material_runtime.h"
#include "scoring/runtime/osh_scoring_estimator_common.h"
#include "scoring/runtime/osh_scoring_estimator_internal.h"

/**
 * @brief Accumulate energy deposition [MeV] into the ENERGY scorer pages.
 *
 * Distributes st->de proportionally to path length in each crossed voxel.
 */
enum osh_status score_step_energy(struct osh_scoring_runtime const *rt,
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
enum osh_status score_step_fluence(struct osh_scoring_runtime const *rt,
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
 * @c vol_inv is read per-crossing from @c crossings[j].vol_inv (pre-filled by
 * the caller: uniform scalar for Mesh; per-R-bin LUT for Cyl).
 *
 * When @c mat_tables is available and a Settings block specifies a medium,
 * applies a stopping-power ratio correction S(ovr,E)/S(tr,E) for
 * dose-to-medium scoring.  Pure density overrides do not change the dose
 * (Fano theorem).
 */
enum osh_status score_step_dose(struct osh_scoring_runtime const *rt,
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
enum osh_status score_step_dlet(struct osh_scoring_runtime const *rt,
                                struct osh_scoring_accumulator *acc_set,
                                struct osh_scoring_geometry_score_group const *group,
                                struct osh_voxel_crossing const *crossings,
                                size_t ncross,
                                struct particle const *part,
                                struct step const *st,
                                double score_len) {
    size_t i;
    size_t j;
    double sp_transport;     /* Stopping power in the transport medium */
    double let_default;      /* Default LET value */
    double let_step;         /* LET value for the current step */
    double mean_energy;      /* Mean energy of the particle */
    double e_per_nuc;        /* Energy per nucleon */
    double rho_ovr;          /* density used for the per-page LET override */
    double dose_weight;      /* Dose weight for the current step */
    double dlet_numerator;   /* Numerator for the dose-averaged LET */
    double dlet_denominator; /* Denominator for the dose-averaged LET */
    size_t proj_idx;
    int have_proj;
    struct osh_scoring_page_runtime const *page;
    struct osh_scoring_accumulator *acc;
    struct osh_scoring_page_override const *sset; /* per-page settings override pointer */
    struct osh_material_runtime const *mat_tables;

    mat_tables = rt->mat_tables;
    e_per_nuc = 0.0; /* initialized to satisfy MSVC C4701; overwritten when table data is available */

    if (part->z == 0 || part->a == 0 || !(st->rho > 0.0)) {
        return OSH_OK;
    }

    /* Determine the LET value carried by this transport step.  The default is
     * table-based LET in the transport medium when tables are available.  If the
     * transport material data are missing, use de/score_len as a geometric
     * fallback so old no-table workflows still score a sensible LET. */
    have_proj = 0;
    proj_idx = 0;
    sp_transport = 0.0;
    let_default = st->de / score_len; /* geometric fallback */
    if (mat_tables && st->medium >= 0) {
        mean_energy = 0.5 * (st->p[3] + st->q[3]);
        e_per_nuc = mean_energy / (double) part->a;
        if (osh_scoring_estimator_find_proj_idx(mat_tables, (unsigned int) part->z, &proj_idx)) {
            have_proj = 1;
            sp_transport = osh_material_runtime_sp_lookup(mat_tables, (size_t) st->medium, proj_idx, e_per_nuc);
            let_default = sp_transport * st->rho;
        }
    }

    for (i = 0; i < group->npages; ++i) {
        page = &rt->pages[group->first_page + i];
        acc = &acc_set[group->first_page + i];
        if (!osh_scoring_page_passes_filters(page, part, st)) {
            continue;
        }
        /* Each page can override the scoring medium or density through Settings.
         * That changes the LET value in the numerator, but not the dose weight:
         * DLET is still weighted by where this step deposited energy. */
        let_step = let_default;
        if (have_proj && page->has_sset) {
            sset = &page->sset;
            if (sset->has_medium && sset->medium >= 0) {
                rho_ovr = sset->has_density_g_cm3 ? sset->density_g_cm3 : mat_tables->rho[sset->medium];
                let_step =
                    osh_material_runtime_sp_lookup(mat_tables, (size_t) sset->medium, proj_idx, e_per_nuc) * rho_ovr;
            } else if (sset->has_density_g_cm3) {
                let_step = sp_transport * sset->density_g_cm3;
            }
        }
        /* DLET stores a dose-weighted average in two explicit accumulators.
         *
         * dose_weight is the energy-deposition share assigned to this crossed
         * spatial bin: st->de * crossing.path_len / score_len [MeV].
         *
         * Numerator:
         *   acc->data[crossing.idx]  += LET [MeV/cm] * dose_weight [MeV]
         *
         * Denominator:
         *   acc->data2[crossing.idx] += dose_weight [MeV]
         *
         * Differential axes are not valid for DLET pages, so crossing.idx is the
         * complete bin index here; no flat diff-bin offset is needed. */
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
 * Same table lookup as score_step_dlet(); uses track-length ds_vox as the weight
 * rather than dose weight.
 *
 * Raw accumulator meaning before postprocess:
 *   acc->data  = sum(LET * track_weight) per spatial bin
 *   acc->data2 = sum(track_weight) per spatial bin
 *
 * postprocess_ratio() later returns acc->data / acc->data2, i.e. LETt.
 */
enum osh_status score_step_tlet(struct osh_scoring_runtime const *rt,
                                struct osh_scoring_accumulator *acc_set,
                                struct osh_scoring_geometry_score_group const *group,
                                struct osh_voxel_crossing const *crossings,
                                size_t ncross,
                                struct particle const *part,
                                struct step const *st,
                                double score_len) {
    size_t i;
    size_t j;
    double sp_transport;
    double let_default;
    double let_step;
    double mean_energy;
    double e_per_nuc;
    double rho_ovr; /* density used for the per-page LET override */
    double track_weight;
    double tlet_numerator;
    double tlet_denominator;
    size_t proj_idx;
    int have_proj;
    struct osh_scoring_page_runtime const *page;
    struct osh_scoring_accumulator *acc;
    struct osh_scoring_page_override const *sset; /* per-page settings override pointer */
    struct osh_material_runtime const *mat_tables;

    mat_tables = rt->mat_tables;
    e_per_nuc = 0.0; /* initialized to satisfy MSVC C4701; overwritten when table data is available */

    if (part->z == 0 || part->a == 0 || !(st->rho > 0.0)) {
        return OSH_OK;
    }

    /* Determine the LET value carried by this transport step, using the same
     * transport-medium default as DLET.  TLET differs only in the denominator:
     * it weights by track length rather than energy deposition. */
    have_proj = 0;
    proj_idx = 0;
    sp_transport = 0.0;
    let_default = st->de / score_len;
    if (mat_tables && st->medium >= 0) {
        mean_energy = 0.5 * (st->p[3] + st->q[3]);
        e_per_nuc = mean_energy / (double) part->a;
        if (osh_scoring_estimator_find_proj_idx(mat_tables, (unsigned int) part->z, &proj_idx)) {
            have_proj = 1;
            sp_transport = osh_material_runtime_sp_lookup(mat_tables, (size_t) st->medium, proj_idx, e_per_nuc);
            let_default = sp_transport * st->rho;
        }
    }

    for (i = 0; i < group->npages; ++i) {
        page = &rt->pages[group->first_page + i];
        acc = &acc_set[group->first_page + i];
        if (!osh_scoring_page_passes_filters(page, part, st)) {
            continue;
        }
        /* Settings can change which medium/density defines the scored LET.
         * The track-weight denominator remains geometric and independent of
         * those material overrides. */
        let_step = let_default;
        if (have_proj && page->has_sset) {
            sset = &page->sset;
            if (sset->has_medium && sset->medium >= 0) {
                rho_ovr = sset->has_density_g_cm3 ? sset->density_g_cm3 : mat_tables->rho[sset->medium];
                let_step =
                    osh_material_runtime_sp_lookup(mat_tables, (size_t) sset->medium, proj_idx, e_per_nuc) * rho_ovr;
            } else if (sset->has_density_g_cm3) {
                let_step = sp_transport * sset->density_g_cm3;
            }
        }
        /* TLET stores a track-length-weighted average in two explicit
         * accumulators.
         *
         * track_weight is the physical step length assigned to this crossed
         * spatial bin: st->ds * crossing.path_len / score_len [cm].
         *
         * Numerator:
         *   acc->data[crossing.idx]  += LET [MeV/cm] * track_weight [cm]
         *
         * Denominator:
         *   acc->data2[crossing.idx] += track_weight [cm]
         *
         * Differential axes are not valid for TLET pages, so crossing.idx is the
         * complete bin index here; no flat diff-bin offset is needed. */
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
enum osh_status score_step_dqeff(struct osh_scoring_runtime const *rt,
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
    size_t proj_idx;
    struct osh_scoring_page_runtime const *page;
    struct osh_scoring_accumulator *acc;
    struct osh_material_runtime const *mat_tables;

    mat_tables = rt->mat_tables;
    if (part->z == 0) {
        return OSH_OK;
    }
    if (!mat_tables) {
        return OSH_OK;
    }
    if (!osh_scoring_estimator_find_proj_idx(mat_tables, (unsigned int) part->z, &proj_idx)) {
        return OSH_OK;
    }

    /* Qeff depends on projectile charge and beta at the step midpoint.  It is a
     * single scalar for this step; only the averaging weight varies per crossing. */
    mean_energy = 0.5 * (st->p[3] + st->q[3]);
    beta = osh_scoring_estimator_particle_beta(mean_energy, mat_tables->projectile_mass_mev[proj_idx]);
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
 * Same as score_step_dqeff() but uses track-length ds_vox as the weight.
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
enum osh_status score_step_tqeff(struct osh_scoring_runtime const *rt,
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
    size_t proj_idx;          /* Column index for given projectile in the stopping-power table */
    struct osh_scoring_page_runtime const *page;
    struct osh_scoring_accumulator *acc;
    struct osh_material_runtime const *mat_tables;
    mat_tables = rt->mat_tables;
    if (part->z == 0) {
        return OSH_OK;
    }
    if (!mat_tables) {
        return OSH_OK;
    }
    if (!osh_scoring_estimator_find_proj_idx(mat_tables, (unsigned int) part->z, &proj_idx)) {
        return OSH_OK;
    }

    /* Qeff is computed once for the step midpoint.  TQEFF then averages that
     * value with geometric track-length weights in each crossed bin. */
    mean_energy = 0.5 * (st->p[3] + st->q[3]);
    beta = osh_scoring_estimator_particle_beta(mean_energy, mat_tables->projectile_mass_mev[proj_idx]);
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
