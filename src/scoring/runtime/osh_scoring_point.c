#include "scoring/runtime/osh_scoring_point.h"

#include <math.h>

#include "common/osh_vect.h"
#include "common/raytrace/osh_raytrace.h"
#include "scoring/runtime/osh_scoring_step_internal.h"

static int zone_bin_index(struct osh_scoring_geometry_runtime const *geo, int zone, size_t *idx_out);

/*
 * Point scorer: deposit a particle's energy at a single location with no track
 * length (c.f. issue #179).  Reuses the step scorer's per-group ENERGY and DOSE
 * deposit helpers (osh_scoring_step_internal.h) so filters, differential axes
 * and dose-to-medium overrides behave identically to score_step — the only
 * difference is a single located voxel instead of a raytraced set of crossings.
 */
enum osh_status osh_scoring_score_point(struct osh_scoring_runtime const *rt,
                                        struct osh_scoring_accumulator *acc_set,
                                        struct osh_scoring_scratch *scratch,
                                        struct particle const *part,
                                        struct step const *st) {
    size_t i;
    size_t g;
    struct osh_raytrace_grid grid = {0};
    enum osh_status rc;

    if (!rt || !part || !st) {
        return OSH_EINVAL;
    }
    if (rt->npages > 0u && !acc_set) {
        return OSH_EINVAL;
    }
    /* Kept in the signature for parity with score_step and forward-compatibility;
     * require it when there is geometry, as score_step does. */
    if (rt->ngeometries > 0u && !scratch) {
        return OSH_EINVAL;
    }
    /* A point deposit carries only energy (st->de) at st->p; nothing else. */
    if (!(st->de > 0.0)) {
        return OSH_OK;
    }

    for (i = 0; i < rt->ngeometries; ++i) {
        struct osh_scoring_geometry_runtime const *geo = &rt->geometries[i];
        double p_local[3];
        double const *p_at;
        double voxel_volume_inv;
        struct osh_voxel_crossing one;
        size_t idx;
        double r_cyl;
        size_t r_bin;
        size_t z_bin;

        if (geo->ngroups == 0u) {
            continue;
        }

        if (geo->geo_kind == OSH_SCORING_GEO_ZONE) {
            if (!zone_bin_index(geo, st->zone, &idx)) {
                continue;
            }
            voxel_volume_inv = geo->zone_vol_inv ? geo->zone_vol_inv[idx] : 1.0;
        } else {
            /* Universe->local rotation (same t[16] layout as score_step). */
            if (geo->has_rotation) {
                osh_vect_trans_point_affine(st->p, p_local, geo->t);
                p_at = p_local;
            } else {
                p_at = st->p;
            }

            if (geo->geo_kind == OSH_SCORING_GEO_CYL) {
                rc = cyl_geometry_to_grid(geo, &grid);
                if (rc != OSH_OK) {
                    return rc;
                }
                /* Cylinder axis is local Z, R the transverse radius.  Uniform dr/dz
                 * bins; flat index z_bin*nr + r_bin matches score_step's layout, and
                 * 1/V is the precomputed per-R-bin value. */
                r_cyl = sqrt(p_at[0] * p_at[0] + p_at[1] * p_at[1]);
                if (r_cyl < grid.origin[0] || p_at[2] < grid.origin[2]) {
                    continue; /* inside the bore or below the axial stack */
                }
                r_bin = (size_t) ((r_cyl - grid.origin[0]) / grid.spacing[0]);
                z_bin = (size_t) ((p_at[2] - grid.origin[2]) / grid.spacing[2]);
                if (r_bin >= grid.n[0] || z_bin >= grid.n[2]) {
                    continue; /* beyond the R or Z extent */
                }
                idx = z_bin * grid.n[0] + r_bin;
                voxel_volume_inv = (geo->cyl_vol_inv && r_bin < geo->cyl_nr) ? geo->cyl_vol_inv[r_bin] : 0.0;
            } else {
                rc = mesh_geometry_to_grid(geo, &grid, &voxel_volume_inv);
                if (rc != OSH_OK) {
                    return rc;
                }
                if (!osh_raytrace_locate(&grid, p_at, &idx)) {
                    continue; /* point lies outside this geometry's mesh */
                }
            }
        }

        /* One unit-"length" crossing at the located voxel.  With score_len == 1
         * the energy helper deposits st->de and the dose helper deposits
         * st->de / (rho * V) [MeV/g] — identical accumulator maths to score_step,
         * without a raytrace. */
        one.idx = idx;
        one.path_len = 1.0;
        one.vol_inv = voxel_volume_inv;

        for (g = 0; g < geo->ngroups; ++g) {
            switch (geo->groups[g].score_kind) {
            case OSH_SCORING_SCORE_ENERGY:
                rc = score_group_energy(rt, acc_set, &geo->groups[g], &one, 1u, part, st, 1.0);
                break;
            case OSH_SCORING_SCORE_DOSE:
            case OSH_SCORING_SCORE_DOSEGY:
                /* Neutral point deposits (e.g. a de-excitation gamma) release
                 * their energy here but carry it away without depositing dose
                 * locally, so they book energy but never dose. */
                rc = (part->charge == 0) ? OSH_OK
                                         : score_group_dose(rt, acc_set, &geo->groups[g], &one, 1u, part, st, 1.0);
                break;
            default:
                /* Fluence / LET / QEFF need a track length; wired in later via
                 * precomputed LETd/LETt-vs-Ekin tables.  Skip for now. */
                rc = OSH_OK;
                break;
            }
            if (rc != OSH_OK) {
                return rc;
            }
        }
    }
    return OSH_OK;
}

static int zone_bin_index(struct osh_scoring_geometry_runtime const *geo, int zone, size_t *idx_out) {
    size_t iz;

    if (!geo || !idx_out || zone < 0) {
        return 0;
    }
    for (iz = 0u; iz < geo->nzone_indices; ++iz) {
        if (geo->zone_indices[iz] == (size_t) zone) {
            *idx_out = iz;
            return 1;
        }
    }
    return 0;
}
