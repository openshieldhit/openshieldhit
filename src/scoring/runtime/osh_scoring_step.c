#include "scoring/runtime/osh_scoring_step.h"

#include <math.h>
#include <stdlib.h>

#include "common/osh_vect.h"
#include "common/raytrace/osh_raytrace.h"
#include "common/raytrace/osh_raytrace_cyl.h"
#include "scoring/runtime/osh_scoring_estimator.h"
#include "scoring/runtime/osh_scoring_geometry_runtime_internal.h"

static void step_scoring_segment(struct step const *st, double dir_out[3], double *len_out);

/* Book one group's step deposit via its estimator handler.  An unknown or
 * non-step estimator returns OSH_ENOTSUP, matching the old switch default. */
static enum osh_status dispatch_step(struct osh_scoring_runtime const *rt,
                                     struct osh_scoring_accumulator *acc_set,
                                     struct osh_scoring_geometry_score_group const *group,
                                     struct osh_voxel_crossing const *crossings,
                                     size_t ncross,
                                     struct particle const *part,
                                     struct step const *st,
                                     double score_len) {
    struct osh_scoring_estimator const *est;

    est = osh_scoring_estimator_for(group->score_kind);
    if (!est || !est->score_step) {
        return OSH_ENOTSUP;
    }
    return est->score_step(rt, acc_set, group, crossings, ncross, part, st, score_len);
}

enum osh_status osh_scoring_score_step(struct osh_scoring_runtime const *rt,
                                       struct osh_scoring_accumulator *acc_set,
                                       struct osh_scoring_scratch *scratch,
                                       struct particle const *part,
                                       struct step const *st) {
    size_t i;
    size_t cap;
    size_t ncross;
    double score_dir[3];
    double score_len;
    struct osh_raytrace_grid grid = {0};
    struct osh_voxel_crossing *crossings;
    enum osh_status rc;
    int hit;

    /* ---- Validate caller-owned state ------------------------------------ */

    if (!rt || !part || !st) {
        return OSH_EINVAL;
    }
    if (rt->npages > 0u && !acc_set) {
        return OSH_EINVAL;
    }
    /* Traversal writes into the caller-owned scratch; require it whenever there
     * is geometry to traverse, alongside the other argument checks. */
    if (rt->ngeometries > 0u && !scratch) {
        return OSH_EINVAL;
    }
    if (!(st->ds > 0.0)) {
        return OSH_EINVAL;
    }

    /* ---- Choose the geometric chord used for scoring -------------------- */

    step_scoring_segment(st, score_dir, &score_len);
    if (!(score_len > 0.0)) {
        return OSH_EINVAL;
    }

    /* ---- Visit each scoring geometry ------------------------------------ */

    for (i = 0; i < rt->ngeometries; ++i) {
        struct osh_scoring_geometry_runtime const *geo = &rt->geometries[i];
        double p_local[3];
        double dir_local[3];
        double const *p_trace;
        double const *dir_trace;
        struct osh_voxel_crossing one;
        size_t g;

        if (geo->ngroups == 0u) {
            continue;
        }

        /* ---- Zone geometry: already classified by transport ------------- */

        if (geo->geo_kind == OSH_SCORING_GEO_ZONE) {
            /* Zone scoring does not raytrace the step.  Transport has already
             * assigned st->zone; map that transport zone id to this scorer's
             * dense Zone-bin index and book one whole-step crossing. */
            if (!osh_scoring_geometry_zone_bin_index(geo, st->zone, &one.idx)) {
                continue;
            }
            one.path_len = score_len;
            one.vol_inv = 1.0; /* unused by estimators; volume is applied in postprocess (geo->bin_vol_inv) */
            for (g = 0; g < geo->ngroups; ++g) {
                rc = dispatch_step(rt, acc_set, &geo->groups[g], &one, 1u, part, st, score_len);
                if (rc != OSH_OK) {
                    return rc;
                }
            }
            continue;
        }

        /* ---- Mesh/Cyl geometry: raytrace in local scoring coordinates ---- */

        /* Mesh and Cyl geometries raytrace in the scoring geometry's local
         * coordinate frame.  Zone was handled above because it has no local
         * coordinate lookup. */
        if (geo->has_rotation) {
            osh_vect_trans_point_affine(st->p, p_local, geo->t);
            osh_vect_trans_vector_affine(score_dir, dir_local, geo->t);
            p_trace = p_local;
            dir_trace = dir_local;
        } else {
            p_trace = st->p;
            dir_trace = score_dir;
        }

        /* ---- Cartesian mesh traversal ----------------------------------- */

        if (geo->geo_kind == OSH_SCORING_GEO_MESH) {
            /* Cartesian mesh: build an X/Y/Z raytrace grid.  Per-bin volume is
             * applied later in postprocess (geo->bin_vol_inv), not per crossing. */
            rc = osh_scoring_geometry_mesh_to_grid(geo, &grid);
            if (rc != OSH_OK) {
                return rc;
            }
            cap = grid.n[0] + grid.n[1] + grid.n[2];
            if (cap == 0u) {
                continue;
            }
            if (cap > scratch->crossing_cap || !scratch->crossing_buf) {
                return OSH_ESTATE;
            }
            crossings = scratch->crossing_buf;
            hit = osh_raytrace_traverse(&grid, p_trace, dir_trace, score_len, crossings, &ncross);
            if (!hit || ncross == 0u) {
                continue;
            }
        } else if (geo->geo_kind == OSH_SCORING_GEO_CYL) {
            /* ---- Cylindrical mesh traversal ----------------------------- */

            /* Cylindrical mesh: trace in R/Z.  Per-R-bin volume is applied later
             * in postprocess (geo->bin_vol_inv), not per crossing. */
            rc = osh_scoring_geometry_cyl_to_grid(geo, &grid);
            if (rc != OSH_OK) {
                return rc;
            }
            cap = 2u * grid.n[0] + grid.n[2];
            if (cap == 0u) {
                continue;
            }
            if (cap > scratch->crossing_cap || !scratch->crossing_buf) {
                return OSH_ESTATE;
            }
            crossings = scratch->crossing_buf;
            hit = osh_raytrace_cyl_traverse(&grid, p_trace, dir_trace, score_len, crossings, &ncross);
            if (!hit || ncross == 0u) {
                continue;
            }
        } else {
            return OSH_ENOTSUP;
        }

        /* ---- Dispatch estimator groups over the located crossings -------- */

        for (g = 0; g < geo->ngroups; ++g) {
            rc = dispatch_step(rt, acc_set, &geo->groups[g], crossings, ncross, part, st, score_len);
            if (rc != OSH_OK) {
                return rc;
            }
        }
    }

    return OSH_OK;
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
