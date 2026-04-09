/*
 * Voxel traversal using the parametric algorithm of Siddon (1985), extended
 * to 3-D by Jacobs et al. and reformulated here in DDA-walk form.
 *
 * Reference: R.L. Siddon, "Fast calculation of the exact radiological path
 * length for a three-dimensional CT array", Med. Phys. 12(2), 252-255, 1985.
 *
 * Difference from the DDA reference (osh_raytrace_simple_msh.c)
 * ----------------------------------------------------------
 * The DDA assumes the ray start lies inside the grid.  Siddon's contribution
 * is the parametric bounding-box clip that correctly locates the true entry
 * and exit points of the ray into the grid before the walk begins.  Once the
 * entry point is known, the per-axis initialisation and the voxel walk are
 * identical to the DDA.
 *
 * The clip works by treating the grid as the intersection of three axis-
 * aligned slabs.  For each slab the parametric entry (t_lo) and exit (t_hi)
 * are computed.  The ray enters the grid at t_entry = max(t_lo_i) and exits
 * at t_exit = min(t_hi_i).  If t_entry >= t_exit the ray misses entirely.
 *
 * Numerical note: tmax[] is computed from entry_p (the position at t_entry)
 * rather than from the original p.  This keeps the subtraction small and
 * avoids cancellation errors when t_entry is large.
 */

#include <float.h>
#include <math.h>

#include "common/raytrace/osh_raytrace.h"

static int same_crossing(double a, double b) {
    double tol = 64.0 * DBL_EPSILON * (1.0 + fmax(fabs(a), fabs(b)));
    return fabs(a - b) <= tol;
}

int osh_raytrace_traverse(struct osh_raytrace_grid const *grid,
                          double const p[3],
                          double const v[3],
                          double ds,
                          struct osh_voxel_crossing *crossings,
                          size_t *n_out) {
    double t_entry;   /* parametric entry into the grid bounding box */
    double t_exit;    /* parametric exit from the grid bounding box */
    double t_lo;      /* near slab crossing for the current axis */
    double t_hi;      /* far slab crossing for the current axis */
    double grid_hi;   /* upper bound of the grid along the current axis [cm] */
    double entry_p;   /* ray position at t_entry along the current axis [cm] */
    int vox[3];       /* starting voxel index along each axis */
    int step[3];      /* +1 or -1 per axis; 0 when v[i]==0 */
    double tmax[3];   /* parametric t to next boundary crossing per axis */
    double tdelta[3]; /* parametric t between consecutive crossings per axis */
    double t;
    double t_next;
    int axis;
    size_t n;
    int i;

    /* --- Parametric bounding-box clip ---------------------------------- *
     *                                                                     *
     * For each axis, compute the parametric interval [t_lo, t_hi] within  *
     * which the ray is inside the grid slab.  The overall entry and exit   *
     * are the intersection of all three slab intervals with [0, ds].       */
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
            /* Ray is parallel to this axis: check if p lies within the slab. */
            if (p[i] < grid->origin[i] || p[i] >= grid_hi) {
                *n_out = 0;
                return 0;
            }
            continue;
        }

        if (t_lo > t_entry)
            t_entry = t_lo;
        if (t_hi < t_exit)
            t_exit = t_hi;

        if (t_entry >= t_exit) {
            *n_out = 0;
            return 0;
        }
    }

    /* --- Starting voxel at the entry point ----------------------------- */
    for (i = 0; i < 3; ++i) {
        entry_p = p[i] + t_entry * v[i];
        vox[i] = (int) ((entry_p - grid->origin[i]) / grid->spacing[i]);
        /* Clamp: floating-point imprecision may place entry_p just outside. */
        if (vox[i] < 0)
            vox[i] = 0;
        if (vox[i] >= (int) grid->n[i])
            vox[i] = (int) grid->n[i] - 1;
    }

    /* --- Per-axis DDA parameters, initialised from the entry point ----- *
     *                                                                     *
     * tmax[i] is computed relative to entry_p so that the subtraction    *
     * (next_plane - entry_p) is small and numerically stable even when    *
     * t_entry is large.  The result is then offset by t_entry to keep     *
     * all t values in the same frame as the original ray parameter.       */
    for (i = 0; i < 3; ++i) {
        entry_p = p[i] + t_entry * v[i];

        if (v[i] > 0.0) {
            step[i] = 1;
            tdelta[i] = grid->spacing[i] / v[i];
            tmax[i] = t_entry + (grid->origin[i] + (vox[i] + 1) * grid->spacing[i] - entry_p) / v[i];
        } else if (v[i] < 0.0) {
            step[i] = -1;
            tdelta[i] = -grid->spacing[i] / v[i];
            tmax[i] = t_entry + (grid->origin[i] + vox[i] * grid->spacing[i] - entry_p) / v[i];
        } else {
            step[i] = 0;
            tdelta[i] = HUGE_VAL;
            tmax[i] = HUGE_VAL;
        }
    }

    t = t_entry;
    n = 0;

    for (;;) {
        if (tmax[0] <= tmax[1] && tmax[0] <= tmax[2])
            axis = 0;
        else if (tmax[1] <= tmax[2])
            axis = 1;
        else
            axis = 2;

        /* Use t_exit (not ds) as the clip: the bounding-box clip already
         * accounts for ds, and t_exit may be tighter if the ray leaves the
         * grid before the full step is exhausted. */
        t_next = (tmax[axis] < t_exit) ? tmax[axis] : t_exit;

        crossings[n].idx = (size_t) vox[0] + grid->n[0] * ((size_t) vox[1] + grid->n[1] * (size_t) vox[2]);
        crossings[n].path_len = t_next - t;
        ++n;

        t = t_next;
        if (t >= t_exit)
            break;

        /* Advance every axis that reaches the same boundary plane. */
        for (axis = 0; axis < 3; ++axis) {
            if (!same_crossing(tmax[axis], t_next))
                continue;

            vox[axis] += step[axis];
            tmax[axis] += tdelta[axis];
            if (vox[axis] < 0 || vox[axis] >= (int) grid->n[axis]) {
                *n_out = n;
                return (n > 0) ? 1 : 0;
            }
        }
    }

    *n_out = n;
    return (n > 0) ? 1 : 0;
}
