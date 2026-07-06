#include "scoring/runtime/osh_scoring_estimator_common.h"

#include <math.h>

/** @brief Relativistic speed beta = v/c from total kinetic energy and rest mass [MeV/c2]. */
double osh_scoring_estimator_particle_beta(double e_kin_mev, double rest_mass_mev) {
    double gamma_inv;

    gamma_inv = rest_mass_mev / (e_kin_mev + rest_mass_mev);
    return sqrt(1.0 - (gamma_inv * gamma_inv));
}

/**
 * @brief Barkas effective charge: z_eff = z * (1 - exp(-125 * beta * z^(-2/3))).
 * @pre   z > 0
 */
static double particle_zeff(int z, double beta) {
    double cbrt_z;

    cbrt_z = cbrt((double) z);
    return z * (1.0 - exp(-125.0 * beta / (cbrt_z * cbrt_z)));
}

/**
 * @brief Dose/track quality factor (z_eff/beta)^2 [dimensionless].
 * @pre   z > 0 (delegates to particle_zeff)
 */
double osh_scoring_estimator_particle_qeff(int z, double beta) {
    double zeff;

    zeff = particle_zeff(z, beta);
    return (zeff * zeff) / (beta * beta);
}

/**
 * @brief Find the column index for projectile atomic number z in the SP table.
 * @return 1 if found and proj_idx_out is set; 0 if z has no entry.
 */
int osh_scoring_estimator_find_proj_idx(struct osh_material_runtime const *tables,
                                        unsigned int z,
                                        size_t *proj_idx_out) {
    size_t i;

    for (i = 0; i < tables->nprojectiles; ++i) {
        if (tables->projectile_z[i] == z) {
            *proj_idx_out = i;
            return 1;
        }
    }
    return 0;
}

/* Compute the 0-based bin index for a differential axis value.
 * Returns diff_nbins (out-of-range sentinel) when val is outside [lo, hi). */
static size_t diff_axis_bin(struct osh_scoring_page_runtime const *page, double val) {
    double frac; /* normalised position within [lo, hi) */
    size_t bin;  /* candidate bin index */

    if (!(val > page->diff_lo) || !(val < page->diff_hi)) {
        return page->diff_nbins; /* sentinel: out of range */
    }
    if (page->diff_log) {
        frac = log10(val / page->diff_lo) / log10(page->diff_hi / page->diff_lo);
    } else {
        frac = (val - page->diff_lo) / (page->diff_hi - page->diff_lo);
    }
    bin = (size_t) floor(frac * (double) page->diff_nbins);
    return (bin < page->diff_nbins) ? bin : page->diff_nbins - 1u;
}

/* Compute LET [MeV/cm] in an explicit medium at step midpoint.
 * Uses the SP table; returns 0 when tables are unavailable or particle is neutral.
 * ovr_medium selects which row in mat_tables to look up (may differ from st->medium).
 * ovr_rho is the density used to convert SP [MeV*cm2/g] to LET [MeV/cm]. */
static double compute_step_let_medium(struct osh_scoring_runtime const *rt,
                                      struct particle const *part,
                                      struct step const *st,
                                      size_t ovr_medium,
                                      double ovr_rho) {
    double mean_energy;
    double e_per_nuc;
    size_t proj_idx;

    if (part->z == 0 || part->a == 0 || !rt->mat_tables) {
        return 0.0;
    }
    mean_energy = 0.5 * (st->p[3] + st->q[3]);
    e_per_nuc = mean_energy / (double) part->a;
    if (!osh_scoring_estimator_find_proj_idx(rt->mat_tables, (unsigned int) part->z, &proj_idx)) {
        return 0.0;
    }
    return osh_material_runtime_sp_lookup(rt->mat_tables, ovr_medium, proj_idx, e_per_nuc) * ovr_rho;
}

/* Compute LET in the transport medium [MeV/cm] at step midpoint.
 * Uses the SP table at mean step energy; falls back to de/ds when tables are unavailable. */
static double
compute_step_let(struct osh_scoring_runtime const *rt, struct particle const *part, struct step const *st) {
    double mean_energy; /* kinetic energy at step midpoint [MeV] */
    double e_per_nuc;   /* mean_energy / A [MeV/u] */
    size_t proj_idx;    /* projectile row in SP table */

    if (part->z == 0 || part->a == 0 || !(st->rho > 0.0)) {
        return 0.0;
    }
    if (rt->mat_tables && st->medium >= 0) {
        mean_energy = 0.5 * (st->p[3] + st->q[3]);
        e_per_nuc = mean_energy / (double) part->a;
        if (osh_scoring_estimator_find_proj_idx(rt->mat_tables, (unsigned int) part->z, &proj_idx)) {
            return osh_material_runtime_sp_lookup(rt->mat_tables, (size_t) st->medium, proj_idx, e_per_nuc) * st->rho;
        }
    }
    return (st->ds > 0.0) ? st->de / st->ds : 0.0;
}

/* Compute the differential LET/DEDX axis value [MeV/cm] with an optional
 * per-axis Settings override. */
static double compute_step_let_with_override(struct osh_scoring_runtime const *rt,
                                             struct particle const *part,
                                             struct step const *st,
                                             struct osh_scoring_page_override const *ovr) {
    double let;
    double rho;

    if (!ovr || (!ovr->has_medium && !ovr->has_density_g_cm3)) {
        return compute_step_let(rt, part, st);
    }

    if (ovr->has_medium && ovr->medium >= 0) {
        if (ovr->has_density_g_cm3) {
            rho = ovr->density_g_cm3;
        } else if (rt->mat_tables) {
            rho = rt->mat_tables->rho[(size_t) ovr->medium];
        } else {
            rho = st->rho;
        }
        return compute_step_let_medium(rt, part, st, (size_t) ovr->medium, rho);
    }

    rho = ovr->density_g_cm3;
    if (rt->mat_tables && st->medium >= 0) {
        return compute_step_let_medium(rt, part, st, (size_t) st->medium, rho);
    }

    let = compute_step_let(rt, part, st);
    if (!(let > 0.0) || !(st->rho > 0.0)) {
        return 0.0;
    }
    return let * (rho / st->rho);
}

/* Compute (z_eff/beta)^2 at step midpoint.
 * Returns 0 for neutrals or when tables are unavailable. */
static double
compute_step_qeff(struct osh_scoring_runtime const *rt, struct particle const *part, struct step const *st) {
    double mean_energy; /* kinetic energy at step midpoint [MeV] */
    double beta;        /* particle velocity / c */
    size_t proj_idx;    /* projectile row in table (provides rest mass) */

    if (part->z == 0) {
        return 0.0;
    }
    if (!rt->mat_tables) {
        return 0.0;
    }
    if (!osh_scoring_estimator_find_proj_idx(rt->mat_tables, (unsigned int) part->z, &proj_idx)) {
        return 0.0;
    }
    mean_energy = 0.5 * (st->p[3] + st->q[3]);
    beta = osh_scoring_estimator_particle_beta(mean_energy, rt->mat_tables->projectile_mass_mev[proj_idx]);
    if (!(beta > 0.0)) {
        return 0.0;
    }
    return osh_scoring_estimator_particle_qeff((int) part->z, beta);
}

/* Compute the diff axis value for a page given the current step.
 * Returns 0 and sets *ok=0 when the value cannot be determined. */
static double diff_step_val(struct osh_scoring_page_runtime const *page,
                            struct osh_scoring_runtime const *rt,
                            struct particle const *part,
                            struct step const *st,
                            int *ok) {
    double mean_ekin; /* kinetic energy at midpoint [MeV] */

    *ok = 1;
    mean_ekin = 0.5 * (st->p[3] + st->q[3]);
    switch (page->diff_kind) {
    case OSH_SCORING_DIFF_EKIN:
        return mean_ekin;
    case OSH_SCORING_DIFF_ENUC:
        if (part->a <= 0) {
            *ok = 0;
            return 0.0;
        }
        return mean_ekin / (double) part->a;
    case OSH_SCORING_DIFF_EAMU:
        if (part->a <= 0) {
            *ok = 0;
            return 0.0;
        }
        return mean_ekin / (double) part->a; /* same as ENUC for integer A */
    case OSH_SCORING_DIFF_LET: {
        double let = page->has_diff_sset ? compute_step_let_with_override(rt, part, st, &page->diff_sset)
                                         : compute_step_let(rt, part, st);
        if (!(let > 0.0)) {
            *ok = 0;
            return 0.0;
        }
        return let;
    }
    case OSH_SCORING_DIFF_QEFF: {
        double qeff = compute_step_qeff(rt, part, st);
        if (!(qeff > 0.0)) {
            *ok = 0;
            return 0.0;
        }
        return qeff;
    }
    default:
        *ok = 0;
        return 0.0;
    }
}

/* Same as diff_axis_bin but operates on the second differential axis (diff2_*) fields. */
static size_t diff2_axis_bin(struct osh_scoring_page_runtime const *page, double val) {
    double frac; /* normalised position within [diff2_lo, diff2_hi) */
    size_t bin;  /* candidate bin index */

    if (!(val > page->diff2_lo) || !(val < page->diff2_hi)) {
        return page->diff2_nbins; /* out-of-range sentinel */
    }
    if (page->diff2_log) {
        frac = log10(val / page->diff2_lo) / log10(page->diff2_hi / page->diff2_lo);
    } else {
        frac = (val - page->diff2_lo) / (page->diff2_hi - page->diff2_lo);
    }
    bin = (size_t) floor(frac * (double) page->diff2_nbins);
    return (bin < page->diff2_nbins) ? bin : page->diff2_nbins - 1u;
}

/* Same as diff_step_val but dispatches on page->diff2_kind. */
static double diff2_step_val(struct osh_scoring_page_runtime const *page,
                             struct osh_scoring_runtime const *rt,
                             struct particle const *part,
                             struct step const *st,
                             int *ok) {
    double mean_ekin;

    *ok = 1;
    mean_ekin = 0.5 * (st->p[3] + st->q[3]);
    switch (page->diff2_kind) {
    case OSH_SCORING_DIFF_EKIN:
        return mean_ekin;
    case OSH_SCORING_DIFF_ENUC:
    case OSH_SCORING_DIFF_EAMU:
        if (part->a <= 0) {
            *ok = 0;
            return 0.0;
        }
        return mean_ekin / (double) part->a;
    case OSH_SCORING_DIFF_LET: {
        double let = page->has_diff2_sset ? compute_step_let_with_override(rt, part, st, &page->diff2_sset)
                                          : compute_step_let(rt, part, st);
        if (!(let > 0.0)) {
            *ok = 0;
            return 0.0;
        }
        return let;
    }
    case OSH_SCORING_DIFF_QEFF: {
        double qeff = compute_step_qeff(rt, part, st);
        if (!(qeff > 0.0)) {
            *ok = 0;
            return 0.0;
        }
        return qeff;
    }
    default:
        *ok = 0;
        return 0.0;
    }
}

/* Resolve the per-page differential-axis bin offsets for this step.
 *
 * LET and QEFF appear here as axis values (Diff1Type/Diff2Type) for spectra of
 * ENERGY/FLUENCE/DOSE/DOSEGY.  They are not score quantities in this helper.
 * Returns 1 with *db / *db2 set, or 0 to skip the page when the axis value is
 * undefined or out of range. */
int osh_scoring_estimator_resolve_diff_bins(struct osh_scoring_page_runtime const *page,
                                            struct osh_scoring_runtime const *rt,
                                            struct particle const *part,
                                            struct step const *st,
                                            size_t *db,
                                            size_t *db2) {
    int dv_ok;
    double dv;

    *db = 0u;
    *db2 = 0u;
    if (page->diff_nbins > 0u) {
        dv = diff_step_val(page, rt, part, st, &dv_ok);
        if (!dv_ok) {
            return 0;
        }
        *db = diff_axis_bin(page, dv);
        if (*db >= page->diff_nbins) {
            return 0;
        }
    }
    if (page->diff2_nbins > 0u) {
        dv = diff2_step_val(page, rt, part, st, &dv_ok);
        if (!dv_ok) {
            return 0;
        }
        *db2 = diff2_axis_bin(page, dv);
        if (*db2 >= page->diff2_nbins) {
            return 0;
        }
    }
    return 1;
}

/* Gather the transport-medium stopping power once per step for dose-to-medium pages. */
struct osh_scoring_dose_sp_ctx osh_scoring_estimator_dose_sp_gather(struct osh_scoring_runtime const *rt,
                                                                    struct particle const *part,
                                                                    struct step const *st) {
    struct osh_scoring_dose_sp_ctx c;
    struct osh_material_runtime const *mat;
    double mean_energy;

    mat = rt->mat_tables;
    c.have_proj = 0;
    c.proj_idx = 0;
    c.sp_tr = 0.0;
    c.e_per_nuc = 0.0;
    if (mat && part->z > 0 && part->a > 0 && st->medium >= 0) {
        mean_energy = 0.5 * (st->p[3] + st->q[3]);
        c.e_per_nuc = mean_energy / (double) part->a;
        if (osh_scoring_estimator_find_proj_idx(mat, (unsigned int) part->z, &c.proj_idx)) {
            c.have_proj = 1;
            c.sp_tr = osh_material_runtime_sp_lookup(mat, (size_t) st->medium, c.proj_idx, c.e_per_nuc);
        }
    }
    return c;
}

/* Dose-to-medium stopping-power ratio S(ovr)/S(tr) for a page's medium override;
 * 1.0 when there is no override (density-only overrides are Fano-invariant). */
double osh_scoring_estimator_dose_sp_ratio(struct osh_scoring_dose_sp_ctx const *c,
                                           struct osh_scoring_runtime const *rt,
                                           struct osh_scoring_page_runtime const *page) {
    double sp_ovr;

    if (c->have_proj && page->has_sset && page->sset.has_medium && page->sset.medium >= 0) {
        sp_ovr = osh_material_runtime_sp_lookup(rt->mat_tables, (size_t) page->sset.medium, c->proj_idx, c->e_per_nuc);
        if (c->sp_tr > 0.0) {
            return sp_ovr / c->sp_tr;
        }
        return 0.0;
    }
    return 1.0;
}

/* Gather the transport-medium LET once per step for the two-pass LET scorers.
 *
 * The fallback here is de/score_len (the geometric scoring chord), which is the
 * score-path convention.  The differential-axis LET path (compute_step_let) uses
 * de/ds instead; the two are kept separate on purpose and only differ for the
 * no-table case on a bent condensed-history step. */
struct osh_scoring_step_let_ctx osh_scoring_estimator_step_let_gather(struct osh_scoring_runtime const *rt,
                                                                      struct particle const *part,
                                                                      struct step const *st,
                                                                      double score_len) {
    struct osh_scoring_step_let_ctx c;
    struct osh_material_runtime const *mat;
    double mean_energy;

    mat = rt->mat_tables;
    c.have_proj = 0;
    c.proj_idx = 0;
    c.e_per_nuc = 0.0;
    c.sp_transport = 0.0;
    c.let_default = st->de / score_len; /* geometric fallback */
    if (mat && st->medium >= 0) {
        mean_energy = 0.5 * (st->p[3] + st->q[3]);
        c.e_per_nuc = mean_energy / (double) part->a;
        if (osh_scoring_estimator_find_proj_idx(mat, (unsigned int) part->z, &c.proj_idx)) {
            c.have_proj = 1;
            c.sp_transport = osh_material_runtime_sp_lookup(mat, (size_t) st->medium, c.proj_idx, c.e_per_nuc);
            c.let_default = c.sp_transport * st->rho;
        }
    }
    return c;
}

/* LET [MeV/cm] a page scores for this step, honouring its Settings override.
 *
 * A per-page override changes only the LET value in the DLET/TLET numerator, not
 * the dose/track weight in the denominator.  Density-only overrides rescale the
 * transport-medium stopping power (Fano-invariant otherwise). */
double osh_scoring_estimator_step_let_apply(struct osh_scoring_step_let_ctx const *ctx,
                                            struct osh_scoring_runtime const *rt,
                                            struct osh_scoring_page_runtime const *page) {
    struct osh_scoring_page_override const *sset;
    double rho_ovr;

    if (!(ctx->have_proj && page->has_sset)) {
        return ctx->let_default;
    }
    sset = &page->sset;
    if (sset->has_medium && sset->medium >= 0) {
        if (sset->has_density_g_cm3) {
            rho_ovr = sset->density_g_cm3;
        } else {
            rho_ovr = rt->mat_tables->rho[sset->medium];
        }
        return osh_material_runtime_sp_lookup(rt->mat_tables, (size_t) sset->medium, ctx->proj_idx, ctx->e_per_nuc)
               * rho_ovr;
    }
    if (sset->has_density_g_cm3) {
        return ctx->sp_transport * sset->density_g_cm3;
    }
    return ctx->let_default;
}
