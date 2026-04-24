#include <float.h>
#include <math.h>

#include "common/osh_voxel_order.h"
#include "common/raytrace/osh_raytrace.h"

static int grid_is_valid(struct osh_raytrace_grid const *grid) {
    int i;

    if (!grid) {
        return 0;
    }

    for (i = 0; i < 3; ++i) {
        if (grid->n[i] == 0u || grid->spacing[i] <= 0.0) {
            return 0;
        }
    }

    return 1;
}

static size_t grid_flat_index(struct osh_raytrace_grid const *grid, size_t ix, size_t iy, size_t iz) {
    if (grid->tile_order == OSH_VOXEL_ORDER_MORTON8) {
        size_t Tx = (grid->n[0] + 7u) >> 3u;
        size_t Ty = (grid->n[1] + 7u) >> 3u;
        return osh_voxel_tile_idx(ix, iy, iz, Tx, Ty);
    }

    return ix + grid->n[0] * (iy + grid->n[1] * iz);
}

static size_t grid_entry_voxel_axis(struct osh_raytrace_grid const *grid, int axis, double coord, double dir) {
    double rel = (coord - grid->origin[axis]) / grid->spacing[axis];
    double lower;
    double upper;
    double tol;
    size_t boundary_idx;
    size_t idx;

    if (rel <= 0.0) {
        return 0u;
    }
    if (rel >= (double) grid->n[axis]) {
        return grid->n[axis] - 1u;
    }

    lower = floor(rel);
    upper = lower + 1.0;
    tol = 64.0 * DBL_EPSILON * (1.0 + fabs(rel));

    if (fabs(rel - lower) <= tol || fabs(upper - rel) <= tol) {
        boundary_idx = (fabs(rel - lower) <= tol) ? (size_t) lower : (size_t) upper;
        if (dir < 0.0 && boundary_idx > 0u) {
            return boundary_idx - 1u;
        }
        if (boundary_idx >= grid->n[axis]) {
            return grid->n[axis] - 1u;
        }
        return boundary_idx;
    }

    idx = (size_t) lower;
    if (idx >= grid->n[axis]) {
        return grid->n[axis] - 1u;
    }
    return idx;
}

int osh_raytrace_locate(struct osh_raytrace_grid const *grid, double const p[3], size_t *idx_out) {
    size_t vox[3];
    double rel;
    int i;

    if (!grid_is_valid(grid) || !p || !idx_out) {
        return 0;
    }

    for (i = 0; i < 3; ++i) {
        rel = (p[i] - grid->origin[i]) / grid->spacing[i];
        if (rel < 0.0 || rel >= (double) grid->n[i]) {
            return 0;
        }

        vox[i] = (size_t) rel;
        if (vox[i] >= grid->n[i]) {
            return 0;
        }
    }

    *idx_out = grid_flat_index(grid, vox[0], vox[1], vox[2]);
    return 1;
}

int osh_raytrace_first_crossing(struct osh_raytrace_grid const *grid,
                                double const p[3],
                                double const v[3],
                                double ds,
                                struct osh_voxel_crossing *out) {
    double t_entry;
    double t_exit;
    double t_lo;
    double t_hi;
    double grid_hi;
    double entry_p;
    size_t vox[3];
    double t_next;
    double t_face;
    int i;

    if (!grid_is_valid(grid) || !p || !v || !out || ds <= 0.0) {
        return 0;
    }

    t_entry = 0.0;
    t_exit = ds;

    for (i = 0; i < 3; ++i) {
        grid_hi = grid->origin[i] + (double) grid->n[i] * grid->spacing[i];

        if (v[i] > 0.0) {
            t_lo = (grid->origin[i] - p[i]) / v[i];
            t_hi = (grid_hi - p[i]) / v[i];
        } else if (v[i] < 0.0) {
            t_lo = (grid_hi - p[i]) / v[i];
            t_hi = (grid->origin[i] - p[i]) / v[i];
        } else {
            if (p[i] < grid->origin[i] || p[i] >= grid_hi) {
                return 0;
            }
            continue;
        }

        if (t_lo > t_entry) {
            t_entry = t_lo;
        }
        if (t_hi < t_exit) {
            t_exit = t_hi;
        }

        if (t_entry >= t_exit) {
            return 0;
        }
    }

    for (i = 0; i < 3; ++i) {
        entry_p = p[i] + t_entry * v[i];
        vox[i] = grid_entry_voxel_axis(grid, i, entry_p, v[i]);
    }

    t_next = t_exit;
    for (i = 0; i < 3; ++i) {
        if (v[i] > 0.0) {
            t_face = (grid->origin[i] + ((double) vox[i] + 1.0) * grid->spacing[i] - p[i]) / v[i];
        } else if (v[i] < 0.0) {
            t_face = (grid->origin[i] + (double) vox[i] * grid->spacing[i] - p[i]) / v[i];
        } else {
            continue;
        }

        if (t_face < t_next) {
            t_next = t_face;
        }
    }

    if (t_next <= t_entry) {
        return 0;
    }

    out->idx = grid_flat_index(grid, vox[0], vox[1], vox[2]);
    out->path_len = t_next - t_entry;
    return 1;
}
