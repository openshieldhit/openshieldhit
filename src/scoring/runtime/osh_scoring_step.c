#include "scoring/runtime/osh_scoring_step.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_vect.h"
#include "common/raytrace/osh_raytrace.h"
#include "common/raytrace/osh_raytrace_cyl.h"
#include "material/runtime/osh_material_runtime.h"

static int axis_index(struct osh_scoring_geometry_runtime const *geo, char const *label);
static void step_scoring_segment(struct step const *st, double dir_out[3], double *len_out);
static int find_proj_idx(struct osh_material_runtime const *tables, unsigned int z, size_t *proj_idx_out);
static double
compute_step_let(struct osh_scoring_runtime const *rt, struct particle const *part, struct step const *st);
static enum osh_status mesh_geometry_to_grid(struct osh_scoring_geometry_runtime const *geo,
                                             struct osh_raytrace_grid *grid,
                                             double *voxel_volume_inv_out);
static enum osh_status cyl_geometry_to_grid(struct osh_scoring_geometry_runtime const *geo,
                                            struct osh_raytrace_grid *grid);
static enum osh_status score_group_energy(struct osh_scoring_runtime *rt,
                                          struct osh_scoring_geometry_score_group const *group,
                                          struct osh_voxel_crossing const *crossings,
                                          size_t ncross,
                                          struct particle const *part,
                                          struct step const *st,
                                          double score_len);
static enum osh_status score_group_fluence(struct osh_scoring_runtime *rt,
                                           struct osh_scoring_geometry_score_group const *group,
                                           struct osh_voxel_crossing const *crossings,
                                           size_t ncross,
                                           struct particle const *part,
                                           struct step const *st,
                                           double score_len);
static enum osh_status score_group_dose(struct osh_scoring_runtime *rt,
                                        struct osh_scoring_geometry_score_group const *group,
                                        struct osh_voxel_crossing const *crossings,
                                        size_t ncross,
                                        struct particle const *part,
                                        struct step const *st,
                                        double score_len);
static enum osh_status score_group_dlet(struct osh_scoring_runtime *rt,
                                        struct osh_scoring_geometry_score_group const *group,
                                        struct osh_voxel_crossing const *crossings,
                                        size_t ncross,
                                        struct particle const *part,
                                        struct step const *st,
                                        double score_len);
static enum osh_status score_group_tlet(struct osh_scoring_runtime *rt,
                                        struct osh_scoring_geometry_score_group const *group,
                                        struct osh_voxel_crossing const *crossings,
                                        size_t ncross,
                                        struct particle const *part,
                                        struct step const *st,
                                        double score_len);
static enum osh_status score_group_dqeff(struct osh_scoring_runtime *rt,
                                         struct osh_scoring_geometry_score_group const *group,
                                         struct osh_voxel_crossing const *crossings,
                                         size_t ncross,
                                         struct particle const *part,
                                         struct step const *st,
                                         double score_len);
static enum osh_status score_group_tqeff(struct osh_scoring_runtime *rt,
                                         struct osh_scoring_geometry_score_group const *group,
                                         struct osh_voxel_crossing const *crossings,
                                         size_t ncross,
                                         struct particle const *part,
                                         struct step const *st,
                                         double score_len);

enum osh_status
osh_scoring_score_step(struct osh_scoring_runtime *rt, struct particle const *part, struct step const *st) {
    size_t i;
    size_t j;
    size_t cap;
    size_t ncross;
    double score_dir[3];
    double score_len;
    struct osh_raytrace_grid grid = {0};
    struct osh_voxel_crossing *crossings;
    enum osh_status rc;
    int hit;
    int is_cyl;

    if (!rt || !part || !st) {
        return OSH_EINVAL;
    }
    if (!(st->ds > 0.0)) {
        return OSH_EINVAL;
    }

    step_scoring_segment(st, score_dir, &score_len);
    if (!(score_len > 0.0)) {
        return OSH_EINVAL;
    }

    for (i = 0; i < rt->ngeometries; ++i) {
        struct osh_scoring_geometry_runtime const *geo = &rt->geometries[i];
        double p_local[3];
        double dir_local[3];
        double const *p_trace;
        double const *dir_trace;
        double voxel_volume_inv;
        double const *lut;
        size_t nr_cyl;
        size_t g;

        if (geo->ngroups == 0u) {
            continue;
        }

        /* Universe→local rotation (shared by Mesh and Cyl: same t[16] layout). */
        if (geo->has_rotation) {
            osh_vect_trans_point_affine(st->p, p_local, geo->t);
            osh_vect_trans_vector_affine(score_dir, dir_local, geo->t);
            p_trace = p_local;
            dir_trace = dir_local;
        } else {
            p_trace = st->p;
            dir_trace = score_dir;
        }

        is_cyl = (geo->geo_kind == OSH_SCORING_GEO_CYL);
        if (!is_cyl) {
            rc = mesh_geometry_to_grid(geo, &grid, &voxel_volume_inv);
            if (rc != OSH_OK) {
                return rc;
            }
            cap = grid.n[0] + grid.n[1] + grid.n[2];
            if (cap == 0u) {
                continue;
            }
            if (cap > rt->crossing_cap || !rt->crossing_buf) {
                return OSH_ESTATE;
            }
            crossings = rt->crossing_buf;
            hit = osh_raytrace_traverse(&grid, p_trace, dir_trace, score_len, crossings, &ncross);
            if (!hit || ncross == 0u) {
                continue;
            }
            for (j = 0; j < ncross; ++j) {
                crossings[j].vol_inv = voxel_volume_inv;
            }
        } else {
            rc = cyl_geometry_to_grid(geo, &grid);
            if (rc != OSH_OK) {
                return rc;
            }
            cap = 2u * grid.n[0] + grid.n[2];
            if (cap == 0u) {
                continue;
            }
            if (cap > rt->crossing_cap || !rt->crossing_buf) {
                return OSH_ESTATE;
            }
            crossings = rt->crossing_buf;
            hit = osh_raytrace_cyl_traverse(&grid, p_trace, dir_trace, score_len, crossings, &ncross);
            if (!hit || ncross == 0u) {
                continue;
            }
            lut = geo->cyl_vol_inv;
            nr_cyl = geo->cyl_nr;
            for (j = 0; j < ncross; ++j) {
                crossings[j].vol_inv = lut[crossings[j].idx % nr_cyl];
            }
        }

        for (g = 0; g < geo->ngroups; ++g) {
            switch (geo->groups[g].score_kind) {
            case OSH_SCORING_SCORE_ENERGY:
                rc = score_group_energy(rt, &geo->groups[g], crossings, ncross, part, st, score_len);
                break;
            case OSH_SCORING_SCORE_FLUENCE:
                rc = score_group_fluence(rt, &geo->groups[g], crossings, ncross, part, st, score_len);
                break;
            case OSH_SCORING_SCORE_DOSE:
            case OSH_SCORING_SCORE_DOSEGY:
                rc = score_group_dose(rt, &geo->groups[g], crossings, ncross, part, st, score_len);
                break;
            case OSH_SCORING_SCORE_DLET:
                rc = score_group_dlet(rt, &geo->groups[g], crossings, ncross, part, st, score_len);
                break;
            case OSH_SCORING_SCORE_TLET:
                rc = score_group_tlet(rt, &geo->groups[g], crossings, ncross, part, st, score_len);
                break;
            case OSH_SCORING_SCORE_DQEFF:
                rc = score_group_dqeff(rt, &geo->groups[g], crossings, ncross, part, st, score_len);
                break;
            case OSH_SCORING_SCORE_TQEFF:
                rc = score_group_tqeff(rt, &geo->groups[g], crossings, ncross, part, st, score_len);
                break;
            default:
                rc = OSH_ENOTSUP;
                break;
            }
            if (rc != OSH_OK) {
                return rc;
            }
        }
    }

    return OSH_OK;
}

enum osh_status
osh_scoring_score_point(struct osh_scoring_runtime *rt, struct particle const *part, struct position const *pos) {
    (void) rt;
    (void) part;
    (void) pos;
    return OSH_ENOTSUP;
}

/**
 * @brief Compute the chord direction and length for the scoring raytrace.
 *
 * Uses the geometric chord p→q; falls back to the step velocity and physical
 * track length ds when the chord is degenerate (strongly bent CH step).
 */
static void step_scoring_segment(struct step const *st, double dir_out[3], double *len_out) {
    double dx;
    double dy;
    double dz;
    double chord_len;

    dx = st->q[0] - st->p[0];
    dy = st->q[1] - st->p[1];
    dz = st->q[2] - st->p[2];
    chord_len = sqrt((dx * dx) + (dy * dy) + (dz * dz));

    if (chord_len > 1.0e-12) {
        dir_out[0] = dx / chord_len;
        dir_out[1] = dy / chord_len;
        dir_out[2] = dz / chord_len;
        *len_out = chord_len;
        return;
    }

    /*
     * Degenerate fallback: a strongly bent condensed-history step can in
     * principle return to its starting point, giving a zero chord but a
     * non-zero physical track length.  Preserve legacy behaviour rather than
     * rejecting the step outright.
     */
    dir_out[0] = st->v[0];
    dir_out[1] = st->v[1];
    dir_out[2] = st->v[2];
    *len_out = st->ds;
}

/** @brief Relativistic speed β = v/c from total kinetic energy and rest mass [MeV/c²]. */
static inline double particle_beta(double e_kin_mev, double rest_mass_mev) {
    double gamma_inv = rest_mass_mev / (e_kin_mev + rest_mass_mev);
    return sqrt(1.0 - (gamma_inv * gamma_inv));
}

/**
 * @brief Barkas effective charge: z_eff = z·(1 − exp(−125·β·z^(−2/3))).
 * @pre   z > 0
 */
static inline double particle_zeff(int z, double beta) {
    double cbrt_z = cbrt((double) z);
    return z * (1.0 - exp(-125.0 * beta / (cbrt_z * cbrt_z)));
}

/**
 * @brief Dose/track quality factor (z_eff/β)² [dimensionless].
 * @pre   z > 0 (delegates to particle_zeff)
 */
static inline double particle_qeff(int z, double beta) {
    double zeff = particle_zeff(z, beta);
    return (zeff * zeff) / (beta * beta);
}

/**
 * @brief Find the column index for projectile atomic number z in the SP table.
 * @return 1 if found and proj_idx_out is set; 0 if z has no entry.
 */
static int find_proj_idx(struct osh_material_runtime const *tables, unsigned int z, size_t *proj_idx_out) {
    size_t i;
    for (i = 0; i < tables->nprojectiles; ++i) {
        if (tables->projectile_z[i] == z) {
            *proj_idx_out = i;
            return 1;
        }
    }
    return 0;
}

/** @brief Return the index of the named axis in geo, or -1 if not found. */
static int axis_index(struct osh_scoring_geometry_runtime const *geo, char const *label) {
    size_t i;

    for (i = 0; i < geo->naxes; ++i) {
        if (strcmp(geo->axes[i].label, label) == 0) {
            return (int) i;
        }
    }
    return -1;
}

/**
 * @brief Fill an osh_raytrace_grid from a mesh scoring geometry and compute 1/voxel_volume.
 *
 * Requires exactly three axes labelled "X", "Y", "Z" with positive bin sizes.
 */
static enum osh_status mesh_geometry_to_grid(struct osh_scoring_geometry_runtime const *geo,
                                             struct osh_raytrace_grid *grid,
                                             double *voxel_volume_inv_out) {
    int ix;
    int iy;
    int iz;
    double dx;
    double dy;
    double dz;

    if (!geo || !grid || !voxel_volume_inv_out) {
        return OSH_EINVAL;
    }
    if (geo->geo_kind != OSH_SCORING_GEO_MESH) {
        return OSH_ENOTSUP;
    }
    if (geo->naxes != 3u) {
        return OSH_EINVAL;
    }

    ix = axis_index(geo, "X");
    iy = axis_index(geo, "Y");
    iz = axis_index(geo, "Z");
    if (ix < 0 || iy < 0 || iz < 0) {
        return OSH_EINVAL;
    }
    if (geo->axes[ix].nbins <= 0 || geo->axes[iy].nbins <= 0 || geo->axes[iz].nbins <= 0) {
        return OSH_EINVAL;
    }

    dx = (geo->axes[ix].hi - geo->axes[ix].lo) / (double) geo->axes[ix].nbins;
    dy = (geo->axes[iy].hi - geo->axes[iy].lo) / (double) geo->axes[iy].nbins;
    dz = (geo->axes[iz].hi - geo->axes[iz].lo) / (double) geo->axes[iz].nbins;
    if (!(dx > 0.0) || !(dy > 0.0) || !(dz > 0.0)) {
        return OSH_EINVAL;
    }

    grid->origin[0] = geo->axes[ix].lo;
    grid->origin[1] = geo->axes[iy].lo;
    grid->origin[2] = geo->axes[iz].lo;
    grid->spacing[0] = dx;
    grid->spacing[1] = dy;
    grid->spacing[2] = dz;
    grid->n[0] = (size_t) geo->axes[ix].nbins;
    grid->n[1] = (size_t) geo->axes[iy].nbins;
    grid->n[2] = (size_t) geo->axes[iz].nbins;
    grid->tile_order = OSH_RAYTRACE_GRID_TILE_ORDER_DEFAULT;
    *voxel_volume_inv_out = 1.0 / (dx * dy * dz);
    return OSH_OK;
}

/**
 * @brief Fill an osh_raytrace_grid from a cylindrical (R,Z) scoring geometry.
 *
 * Requires exactly two axes labelled "R" and "Z"; declaration order in
 * detect.dat does not matter.
 * Field convention: origin/spacing/n[0] = r_min/dr/nr;
 * origin/spacing/n[2] = z_min/dz/nz; origin/spacing/n[1] = 0/0/1 (unused).
 */
static enum osh_status cyl_geometry_to_grid(struct osh_scoring_geometry_runtime const *geo,
                                            struct osh_raytrace_grid *grid) {
    int ir;
    int iz;
    double dr;
    double dz;

    if (!geo || !grid) {
        return OSH_EINVAL;
    }
    if (geo->geo_kind != OSH_SCORING_GEO_CYL) {
        return OSH_ENOTSUP;
    }
    if (geo->naxes != 2u) {
        return OSH_EINVAL;
    }

    ir = axis_index(geo, "R");
    iz = axis_index(geo, "Z");
    if (ir < 0 || iz < 0) {
        return OSH_EINVAL;
    }
    if (geo->axes[ir].nbins <= 0 || geo->axes[iz].nbins <= 0) {
        return OSH_EINVAL;
    }

    dr = (geo->axes[ir].hi - geo->axes[ir].lo) / (double) geo->axes[ir].nbins;
    dz = (geo->axes[iz].hi - geo->axes[iz].lo) / (double) geo->axes[iz].nbins;
    if (!(dr > 0.0) || !(dz > 0.0)) {
        return OSH_EINVAL;
    }

    grid->origin[0] = geo->axes[ir].lo;
    grid->origin[1] = 0.0;
    grid->origin[2] = geo->axes[iz].lo;
    grid->spacing[0] = dr;
    grid->spacing[1] = 0.0;
    grid->spacing[2] = dz;
    grid->n[0] = (size_t) geo->axes[ir].nbins;
    grid->n[1] = 1u;
    grid->n[2] = (size_t) geo->axes[iz].nbins;
    grid->tile_order = OSH_RAYTRACE_GRID_TILE_ORDER_DEFAULT;
    return OSH_OK;
}

/* ---- Differential scoring helpers ---------------------------------------- */

/* Compute the 0-based bin index for a differential axis value.
 * Returns diff_nbins (out-of-range sentinel) when val is outside [lo, hi). */
static inline size_t diff_axis_bin(struct osh_scoring_page_runtime const *page, double val) {
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
 * ovr_rho is the density used to convert SP [MeV·cm²/g] → LET [MeV/cm]. */
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
    if (!find_proj_idx(rt->mat_tables, (unsigned int) part->z, &proj_idx)) {
        return 0.0;
    }
    return osh_material_runtime_sp_lookup(rt->mat_tables, ovr_medium, proj_idx, e_per_nuc) * ovr_rho;
}

/* Compute the differential LET/DEDX axis value [MeV/cm] with an optional
 * per-axis Settings override.
 *
 * Semantics:
 *   no override                -> rho_transport * S_transport(E)
 *   density override only      -> rho_override  * S_transport(E)
 *   medium override only       -> rho_medium    * S_medium(E)
 *   medium + density override  -> rho_override  * S_medium(E)
 *
 * When a different medium is requested we need SP tables for that medium. If
 * only the density is overridden and no tables are available, fall back to the
 * transport-medium LET/de-ds and scale it by rho_override / rho_transport. */
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
        rho = ovr->has_density_g_cm3 ? ovr->density_g_cm3
                                     : (rt->mat_tables ? rt->mat_tables->rho[(size_t) ovr->medium] : st->rho);
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

/* Compute LET in the transport medium [MeV/cm] at step midpoint.
 * Uses the SP table at mean step energy; falls back to de/ds when tables are unavailable.
 * No per-page Settings override — callers that need override apply it themselves. */
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
        if (find_proj_idx(rt->mat_tables, (unsigned int) part->z, &proj_idx)) {
            return osh_material_runtime_sp_lookup(rt->mat_tables, (size_t) st->medium, proj_idx, e_per_nuc) * st->rho;
        }
    }
    return (st->ds > 0.0) ? st->de / st->ds : 0.0;
}

/* Compute (z_eff/β)² at step midpoint.
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
    if (!find_proj_idx(rt->mat_tables, (unsigned int) part->z, &proj_idx)) {
        return 0.0;
    }
    mean_energy = 0.5 * (st->p[3] + st->q[3]);
    beta = particle_beta(mean_energy, rt->mat_tables->projectile_mass_mev[proj_idx]);
    if (!(beta > 0.0)) {
        return 0.0;
    }
    return particle_qeff((int) part->z, beta);
}

/* Compute the diff axis value for a page given the current step.
 * Returns 0 and sets *ok=0 when the value cannot be determined
 * (e.g. LET for a neutral particle). */
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

/* Same as diff_axis_bin but operates on the second differential axis (diff2_*) fields.
 * Returns diff2_nbins as an out-of-range sentinel when val is outside [diff2_lo, diff2_hi). */
static inline size_t diff2_axis_bin(struct osh_scoring_page_runtime const *page, double val) {
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

/* Same as diff_step_val but dispatches on page->diff2_kind.
 * Returns 0 and sets *ok=0 when the value cannot be determined
 * (e.g. LET for a neutral particle, or unknown kind). */
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

/**
 * @brief Accumulate energy deposition [MeV] into the ENERGY scorer pages.
 *
 * Distributes st->de proportionally to path length in each crossed voxel.
 */
static enum osh_status score_group_energy(struct osh_scoring_runtime *rt,
                                          struct osh_scoring_geometry_score_group const *group,
                                          struct osh_voxel_crossing const *crossings,
                                          size_t ncross,
                                          struct particle const *part,
                                          struct step const *st,
                                          double score_len) {
    size_t i;
    size_t j;
    size_t db;
    size_t db2;
    int dv_ok;
    double dv;
    double frac;
    struct osh_scoring_page_runtime *page;

    for (i = 0; i < group->npages; ++i) {
        page = &rt->pages[group->first_page + i];
        if (!osh_scoring_page_passes_filters(page, part, st)) {
            continue;
        }
        /* Determine diff bin indices.  For non-differential pages diff_nbins == 0
         * so both db and db2 stay 0, and the index offsets below collapse to 0. */
        db = 0u;
        db2 = 0u;
        if (page->diff_nbins > 0u) {
            /* Skip the page when the axis value cannot be determined (e.g. LET
             * for a neutral particle) or lies outside the configured [lo, hi). */
            dv = diff_step_val(page, rt, part, st, &dv_ok);
            if (!dv_ok) {
                continue;
            }
            db = diff_axis_bin(page, dv);
            if (db >= page->diff_nbins) {
                continue;
            }
        }
        if (page->diff2_nbins > 0u) {
            dv = diff2_step_val(page, rt, part, st, &dv_ok);
            if (!dv_ok) {
                continue;
            }
            db2 = diff2_axis_bin(page, dv);
            if (db2 >= page->diff2_nbins) {
                continue;
            }
        }
        for (j = 0; j < ncross; ++j) {
            if (crossings[j].idx >= page->diff_stride) {
                return OSH_ESTATE;
            }
            /* Flat index: spatial_idx + db * diff_stride + db2 * diff2_stride.
             * When differential axes are inactive the extra terms evaluate to 0. */
            frac = crossings[j].path_len / score_len;
            page->data[crossings[j].idx + (db * page->diff_stride) + (db2 * page->diff2_stride)] += st->de * frac;
        }
    }
    return OSH_OK;
}

/**
 * @brief Accumulate fluence [1/cm²] into the FLUENCE scorer pages.
 *
 * Scores track_length * vol_inv per voxel crossing.  vol_inv is pre-filled into
 * each crossing by the caller (uniform for Mesh; per-R-bin LUT for Cyl).
 */
static enum osh_status score_group_fluence(struct osh_scoring_runtime *rt,
                                           struct osh_scoring_geometry_score_group const *group,
                                           struct osh_voxel_crossing const *crossings,
                                           size_t ncross,
                                           struct particle const *part,
                                           struct step const *st,
                                           double score_len) {
    size_t i;
    size_t j;
    size_t db;  /* diff1 bin index for current page (0 when no diff axis) */
    size_t db2; /* diff2 bin index for current page (0 when no second diff axis) */
    int dv_ok;  /* flag: differential axis value is physically valid for this step */
    double dv;  /* value on the diff axis (energy, LET, or Qeff at step midpoint) */
    struct osh_scoring_page_runtime *page;
    (void) score_len; /* unused for non-differential fluence; kept for diff LET/QEFF */

    for (i = 0; i < group->npages; ++i) {
        page = &rt->pages[group->first_page + i];
        if (!osh_scoring_page_passes_filters(page, part, st)) {
            continue;
        }
        /* Determine diff bin indices.  For non-differential pages diff_nbins == 0
         * so both db and db2 stay 0, and the index offsets below collapse to 0. */
        db = 0u;
        db2 = 0u;
        if (page->diff_nbins > 0u) {
            /* Skip the page when the axis value cannot be determined (e.g. LET
             * for a neutral particle) or lies outside the configured [lo, hi). */
            dv = diff_step_val(page, rt, part, st, &dv_ok);
            if (!dv_ok) {
                continue;
            }
            db = diff_axis_bin(page, dv);
            if (db >= page->diff_nbins) {
                continue;
            }
        }
        if (page->diff2_nbins > 0u) {
            dv = diff2_step_val(page, rt, part, st, &dv_ok);
            if (!dv_ok) {
                continue;
            }
            db2 = diff2_axis_bin(page, dv);
            if (db2 >= page->diff2_nbins) {
                continue;
            }
        }
        for (j = 0; j < ncross; ++j) {
            if (crossings[j].idx >= page->diff_stride) {
                return OSH_ESTATE;
            }
            /* Flat index: spatial_idx + db * diff_stride + db2 * diff2_stride.
             * When differential axes are inactive the extra terms evaluate to 0. */
            page->data[crossings[j].idx + (db * page->diff_stride) + (db2 * page->diff2_stride)] +=
                crossings[j].path_len * crossings[j].vol_inv;
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
static enum osh_status score_group_dose(struct osh_scoring_runtime *rt,
                                        struct osh_scoring_geometry_score_group const *group,
                                        struct osh_voxel_crossing const *crossings,
                                        size_t ncross,
                                        struct particle const *part,
                                        struct step const *st,
                                        double score_len) {
    size_t i;
    size_t j;
    size_t db;  /* diff1 bin index for current page (0 when no diff axis) */
    size_t db2; /* diff2 bin index for current page (0 when no second diff axis) */
    int dv_ok;  /* flag: differential axis value is physically valid for this step */
    double dv;  /* value on the diff axis (energy, LET, or Qeff at step midpoint) */
    double base_scale;
    double dose_scale;
    double mean_energy;
    double e_per_nuc;
    double sp_tr;
    double sp_ovr;
    size_t proj_idx;
    int have_proj;
    struct osh_scoring_page_runtime *page;
    struct osh_scoring_page_override const *sset; /* per-page settings override pointer */
    struct osh_material_runtime const *mat_tables = rt->mat_tables;
    e_per_nuc = 0.0; /* initialized to satisfy MSVC C4701; overwritten when table-based projectile data is available */

    if (!(st->rho > 0.0)) {
        return OSH_OK;
    }
    /* base_scale = de / (rho * score_len)  — vol_inv applied per crossing below.
     * Accumulates in [MeV/g]; postprocess converts to [Gy] (× OSH_MEVG2GY). */
    base_scale = st->de / (score_len * st->rho);

    /* Precompute projectile index and transport SP once per step. */
    have_proj = 0;
    proj_idx = 0;
    sp_tr = 0.0;
    if (mat_tables && part->z > 0 && part->a > 0 && st->medium >= 0) {
        mean_energy = 0.5 * (st->p[3] + st->q[3]);
        e_per_nuc = mean_energy / (double) part->a;
        if (find_proj_idx(mat_tables, (unsigned int) part->z, &proj_idx)) {
            have_proj = 1;
            sp_tr = osh_material_runtime_sp_lookup(mat_tables, (size_t) st->medium, proj_idx, e_per_nuc);
        }
    }

    for (i = 0; i < group->npages; ++i) {
        page = &rt->pages[group->first_page + i];
        if (!osh_scoring_page_passes_filters(page, part, st)) {
            continue;
        }
        dose_scale = base_scale;
        if (have_proj && page->has_sset) {
            sset = &page->sset;
            if (sset->has_medium && sset->medium >= 0) {
                /* Dose-to-medium: multiply by stopping-power ratio S(ovr)/S(tr). */
                sp_ovr = osh_material_runtime_sp_lookup(mat_tables, (size_t) sset->medium, proj_idx, e_per_nuc);
                dose_scale *= (sp_tr > 0.0) ? sp_ovr / sp_tr : 0.0;
            }
            /* Density-only override: Fano theorem — dose is density-independent,
             * so no correction is needed for pure density overrides. */
        }
        /* Determine diff bin indices.  For non-differential pages diff_nbins == 0
         * so both db and db2 stay 0, and the index offsets below collapse to 0. */
        db = 0u;
        db2 = 0u;
        if (page->diff_nbins > 0u) {
            /* Skip the page when the axis value cannot be determined (e.g. LET
             * for a neutral particle) or lies outside the configured [lo, hi). */
            dv = diff_step_val(page, rt, part, st, &dv_ok);
            if (!dv_ok) {
                continue;
            }
            db = diff_axis_bin(page, dv);
            if (db >= page->diff_nbins) {
                continue;
            }
        }
        if (page->diff2_nbins > 0u) {
            dv = diff2_step_val(page, rt, part, st, &dv_ok);
            if (!dv_ok) {
                continue;
            }
            db2 = diff2_axis_bin(page, dv);
            if (db2 >= page->diff2_nbins) {
                continue;
            }
        }
        for (j = 0; j < ncross; ++j) {
            if (crossings[j].idx >= page->diff_stride) {
                return OSH_ESTATE;
            }
            /* Flat index: spatial_idx + db * diff_stride + db2 * diff2_stride.
             * When differential axes are inactive the extra terms evaluate to 0. */
            page->data[crossings[j].idx + (db * page->diff_stride) + (db2 * page->diff2_stride)] +=
                crossings[j].path_len * crossings[j].vol_inv * dose_scale;
        }
    }
    return OSH_OK;
}

/**
 * @brief Accumulate dose-averaged LET [MeV/cm] via a two-pass accumulator.
 *
 * Uses S(medium,E)·ρ from the SP tables when available; falls back to de/score_len.
 * Per-page Settings overrides apply S(ovr,E)·ρ_ovr (medium) or S(tr,E)·ρ_ovr (density-only).
 */
static enum osh_status score_group_dlet(struct osh_scoring_runtime *rt,
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
    double w;
    size_t proj_idx;
    int have_proj;
    struct osh_scoring_page_runtime *page;
    struct osh_scoring_page_override const *sset; /* per-page settings override pointer */
    struct osh_material_runtime const *mat_tables = rt->mat_tables;
    e_per_nuc = 0.0; /* initialized to satisfy MSVC C4701; overwritten when table-based projectile data is available */

    /* LET is only defined for charged particles in a material (not vacuum). */
    if (part->z == 0 || part->a == 0 || !(st->rho > 0.0)) {
        return OSH_OK;
    }

    /* Compute table-based LET in the transport medium: S(medium,E)·ρ [MeV/cm].
     * sp_transport is kept separately so density-only overrides can reuse it.
     * Fall back to geometric de/score_len when tables are unavailable. */
    have_proj = 0;
    proj_idx = 0;
    sp_transport = 0.0;
    let_default = st->de / score_len; /* geometric fallback */
    if (mat_tables && st->medium >= 0) {
        mean_energy = 0.5 * (st->p[3] + st->q[3]);
        e_per_nuc = mean_energy / (double) part->a;
        if (find_proj_idx(mat_tables, (unsigned int) part->z, &proj_idx)) {
            have_proj = 1;
            sp_transport = osh_material_runtime_sp_lookup(mat_tables, (size_t) st->medium, proj_idx, e_per_nuc);
            let_default = sp_transport * st->rho;
        }
    }

    for (i = 0; i < group->npages; ++i) {
        page = &rt->pages[group->first_page + i];
        if (!osh_scoring_page_passes_filters(page, part, st)) {
            continue;
        }
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
        for (j = 0; j < ncross; ++j) {
            if (crossings[j].idx >= page->diff_stride) {
                return OSH_ESTATE;
            }
            /* Two-pass dose-averaged LET:
             *   data  = sum(LET · dose_weight)   [MeV/cm · MeV]
             *   data2 = sum(dose_weight)          [MeV]
             * osh_scoring_postprocess() divides data by data2 to yield LETd. */
            w = st->de * crossings[j].path_len / score_len;
            page->data[crossings[j].idx] += let_step * w;
            page->data2[crossings[j].idx] += w;
        }
    }
    return OSH_OK;
}

/**
 * @brief Accumulate track-averaged LET [MeV/cm] via a two-pass accumulator.
 *
 * Same table lookup as score_group_dlet(); uses track-length ds_vox as the weight
 * rather than dose weight.
 */
static enum osh_status score_group_tlet(struct osh_scoring_runtime *rt,
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
    double ds_vox;
    size_t proj_idx;
    int have_proj;
    struct osh_scoring_page_runtime *page;
    struct osh_scoring_page_override const *sset; /* per-page settings override pointer */
    struct osh_material_runtime const *mat_tables = rt->mat_tables;
    e_per_nuc = 0.0; /* initialized to satisfy MSVC C4701; overwritten when table-based projectile data is available */

    /* LET is only defined for charged particles in a material (not vacuum). */
    if (part->z == 0 || part->a == 0 || !(st->rho > 0.0)) {
        return OSH_OK;
    }

    /* Table-based LET in transport medium; geometric fallback when unavailable. */
    have_proj = 0;
    proj_idx = 0;
    sp_transport = 0.0;
    let_default = st->de / score_len;
    if (mat_tables && st->medium >= 0) {
        mean_energy = 0.5 * (st->p[3] + st->q[3]);
        e_per_nuc = mean_energy / (double) part->a;
        if (find_proj_idx(mat_tables, (unsigned int) part->z, &proj_idx)) {
            have_proj = 1;
            sp_transport = osh_material_runtime_sp_lookup(mat_tables, (size_t) st->medium, proj_idx, e_per_nuc);
            let_default = sp_transport * st->rho;
        }
    }

    for (i = 0; i < group->npages; ++i) {
        page = &rt->pages[group->first_page + i];
        if (!osh_scoring_page_passes_filters(page, part, st)) {
            continue;
        }
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
        for (j = 0; j < ncross; ++j) {
            if (crossings[j].idx >= page->diff_stride) {
                return OSH_ESTATE;
            }
            /* Two-pass track-averaged LET:
             *   data  = sum(LET · ds_vox)   [MeV/cm · cm = MeV]
             *   data2 = sum(ds_vox)         [cm]
             * osh_scoring_postprocess() divides data by data2 to yield LETt. */
            ds_vox = st->ds * crossings[j].path_len / score_len;
            page->data[crossings[j].idx] += let_step * ds_vox;
            page->data2[crossings[j].idx] += ds_vox;
        }
    }
    return OSH_OK;
}

/**
 * @brief Accumulate dose-averaged (z_eff/β)² [dimensionless] via a two-pass accumulator.
 *
 * Requires mat_tables for the rest mass needed to compute β.  Neutrals and
 * particles without SP table entries contribute nothing.
 *
 * @ref Kalholm et al. Medical Physics. 2023 Jan;50(1):651-9.
 *      https://doi.org/10.1002/mp.16029
 *
 */
static enum osh_status score_group_dqeff(struct osh_scoring_runtime *rt,
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
    double w;
    size_t proj_idx;
    struct osh_scoring_page_runtime *page;
    struct osh_material_runtime const *mat_tables = rt->mat_tables;

    if (part->z == 0) {
        return OSH_OK;
    }
    if (!mat_tables) {
        return OSH_OK;
    }
    if (!find_proj_idx(mat_tables, (unsigned int) part->z, &proj_idx)) {
        return OSH_OK;
    }

    mean_energy = 0.5 * (st->p[3] + st->q[3]);
    beta = particle_beta(mean_energy, mat_tables->projectile_mass_mev[proj_idx]);
    if (!(beta > 0.0)) {
        return OSH_OK;
    }
    qeff = particle_qeff((int) part->z, beta);

    for (i = 0; i < group->npages; ++i) {
        page = &rt->pages[group->first_page + i];
        if (!osh_scoring_page_passes_filters(page, part, st)) {
            continue;
        }
        for (j = 0; j < ncross; ++j) {
            if (crossings[j].idx >= page->diff_stride) {
                return OSH_ESTATE;
            }
            /* Two-pass dose-averaged (z_eff/β)²:
             *   data  = sum(qeff · dose_weight)
             *   data2 = sum(dose_weight)           [MeV]
             * osh_scoring_postprocess() divides to yield the dose-averaged value. */
            w = st->de * crossings[j].path_len / score_len;
            page->data[crossings[j].idx] += qeff * w;
            page->data2[crossings[j].idx] += w;
        }
    }
    return OSH_OK;
}

/**
 * @brief Accumulate track-averaged (z_eff/β)² [dimensionless] via a two-pass accumulator.
 *
 * Same as score_group_dqeff() but uses track-length ds_vox as the weight.
 *
 * @ref Kalholm et al. Medical Physics. 2023 Jan;50(1):651-9.
 *      https://doi.org/10.1002/mp.16029
 *
 */
static enum osh_status score_group_tqeff(struct osh_scoring_runtime *rt,
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
    double ds_vox;
    size_t proj_idx;
    struct osh_scoring_page_runtime *page;
    struct osh_material_runtime const *mat_tables = rt->mat_tables;

    if (part->z == 0) {
        return OSH_OK;
    }
    if (!mat_tables) {
        return OSH_OK;
    }
    if (!find_proj_idx(mat_tables, (unsigned int) part->z, &proj_idx)) {
        return OSH_OK;
    }

    mean_energy = 0.5 * (st->p[3] + st->q[3]);
    beta = particle_beta(mean_energy, mat_tables->projectile_mass_mev[proj_idx]);
    if (!(beta > 0.0)) {
        return OSH_OK;
    }
    qeff = particle_qeff((int) part->z, beta);

    for (i = 0; i < group->npages; ++i) {
        page = &rt->pages[group->first_page + i];
        if (!osh_scoring_page_passes_filters(page, part, st)) {
            continue;
        }
        for (j = 0; j < ncross; ++j) {
            if (crossings[j].idx >= page->diff_stride) {
                return OSH_ESTATE;
            }
            /* Two-pass track-averaged (z_eff/β)²:
             *   data  = sum(qeff · ds_vox)
             *   data2 = sum(ds_vox)           [cm]
             * osh_scoring_postprocess() divides to yield the track-averaged value. */
            ds_vox = st->ds * crossings[j].path_len / score_len;
            page->data[crossings[j].idx] += qeff * ds_vox;
            page->data2[crossings[j].idx] += ds_vox;
        }
    }
    return OSH_OK;
}
