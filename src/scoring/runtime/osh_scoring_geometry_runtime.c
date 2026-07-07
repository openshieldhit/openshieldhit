#include <string.h>

#include "common/raytrace/osh_raytrace.h"
#include "scoring/runtime/osh_scoring_geometry_runtime_internal.h"

/*
 * Scoring-geometry runtime helpers: turn a compiled scoring geometry into the
 * raytrace grid the hot path traverses, and map transport zone ids to dense
 * Zone-scorer bins.  These are shared by the step and point geometry drivers
 * (osh_scoring_step.c / osh_scoring_point.c); declared in
 * osh_scoring_geometry_runtime_internal.h and kept here — where the header
 * implies — rather than inline in the step driver.
 */

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

int osh_scoring_geometry_zone_bin_index(struct osh_scoring_geometry_runtime const *geo, int zone, size_t *idx_out) {
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

/**
 * @brief Fill an osh_raytrace_grid from a mesh scoring geometry.
 *
 * Requires exactly three axes labelled "X", "Y", "Z" with positive bin sizes.
 * Per-bin volume is not returned here; volume-normalised estimators divide by
 * geo->bin_vol_inv in postprocess.
 */
enum osh_status osh_scoring_geometry_mesh_to_grid(struct osh_scoring_geometry_runtime const *geo,
                                                  struct osh_raytrace_grid *grid) {
    int ix;
    int iy;
    int iz;
    double dx;
    double dy;
    double dz;

    if (!geo || !grid) {
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
enum osh_status osh_scoring_geometry_cyl_to_grid(struct osh_scoring_geometry_runtime const *geo,
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
