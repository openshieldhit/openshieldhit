/*
 * Voxel traversal using the alpha-parametric formulation of Jacobs (1998).
 *
 * Reference: F. Jacobs, E. Sundermann, B. De Sutter, M. Christiaens,
 * I. Lemahieu, "A fast algorithm to trace rays in the presence of very
 * large numbers of objects", IEEE Trans. Nucl. Sci. 45(3):907-912, 1998.
 *
 * Difference from osh_raytrace_simple_msh.c and osh_raytrace_siddon_msh.c
 * ---------------------------------------------------------------
 * Rather than using the absolute parametric distance t in [0, ds], this
 * implementation follows Jacobs' original dimensionless alpha in [0, 1]:
 *
 *   pos(alpha) = p + alpha * d,   d[i] = v[i] * ds.
 *
 * The bounding-box clip (identical to Siddon's) produces alpha_entry and
 * alpha_exit.  Per-axis alpha increments are au[i] = spacing[i] / |d[i]|
 * and ai[i] holds the alpha of the next boundary crossing along each axis.
 *
 * Output path_len for each voxel is (delta_alpha) * ds, converting the
 * dimensionless interval back to a physical length.
 *
 * Starting voxel
 * --------------
 * An EPS/|d[i]| offset is added to alpha_entry before computing the
 * starting voxel index.  This places the probe point strictly inside the
 * entry voxel when alpha_entry lands exactly on a plane crossing.
 */

#include <float.h>
#include <math.h>

#include "common/osh_voxel_order.h"
#include "common/raytrace/osh_raytrace.h"

#define JACOBS_EPS 1.0e-9 /* probe offset to step inside the entry voxel */

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
    double alpha_entry = 0.0; /* alpha at entry into the grid bounding box */
    double alpha_exit = 1.0;  /* alpha at exit from the grid bounding box  */
    double d[3];              /* full displacement d[i] = v[i] * ds           (Jacobs: DXS, DYS, DZS) */
    double au[3];             /* alpha increment per voxel crossing on axis i  (Jacobs eq. 28: AX, AY, AZ) */
    int iu[3];                /* step direction: +1 or -1 per axis             (Jacobs eq. 29: XINC, YINC, ZINC) */
    double ai[3];             /* alpha to the next boundary crossing per axis  (Jacobs: AX2, AY2, AZ2) */
    int ii[3];                /* current voxel index along each axis           (Jacobs: KX, KY, KZ) */
    double ac;                /* alpha at the start of the current segment     (Jacobs: AC) */
    double a_next;            /* alpha at the end of the current segment */
    int axis;
    size_t n;
    int i;

    /* --- Displacement vector and bounding-box clip (Jacobs eq. 3-8) --- *
     *                                                                    *
     * For each axis compute the alpha interval [a_lo, a_hi] within which *
     * the ray is inside the grid slab.  The overall entry and exit are   *
     * the intersection of all slab intervals with [0, 1].                */
    for (i = 0; i < 3; ++i) {
        double grid_lo = grid->origin[i];
        double grid_hi = grid->origin[i] + (double) grid->n[i] * grid->spacing[i];
        double a_lo, a_hi;

        d[i] = v[i] * ds;

        if (d[i] > 0.0) {
            a_lo = (grid_lo - p[i]) / d[i]; /* alpha at near slab */
            a_hi = (grid_hi - p[i]) / d[i]; /* alpha at far slab  */
        } else if (d[i] < 0.0) {
            a_lo = (grid_hi - p[i]) / d[i];
            a_hi = (grid_lo - p[i]) / d[i];
        } else {
            /* Ray parallel to this axis: check that p lies inside the slab. */
            if (p[i] < grid_lo || p[i] >= grid_hi) {
                *n_out = 0;
                return 0;
            }
            continue; /* no contribution to alpha_entry / alpha_exit */
        }

        if (a_lo > alpha_entry)
            alpha_entry = a_lo;
        if (a_hi < alpha_exit)
            alpha_exit = a_hi;

        if (alpha_entry >= alpha_exit) {
            *n_out = 0;
            return 0;
        }
    }

    /* --- Starting voxel at alpha_entry (Jacobs eq. 19) ---------------- *
     *                                                                    *
     * Add JACOBS_EPS / |d[i]| to alpha_entry so the probe point lies    *
     * strictly inside the entry voxel and not on its boundary.           */
    for (i = 0; i < 3; ++i) {
        double alpha_eps = alpha_entry + (d[i] != 0.0 ? JACOBS_EPS / fabs(d[i]) : 0.0);
        double entry_x = p[i] + alpha_eps * d[i];
        int idx = (int) ((entry_x - grid->origin[i]) / grid->spacing[i]);
        if (idx < 0)
            idx = 0;
        if (idx >= (int) grid->n[i])
            idx = (int) grid->n[i] - 1;
        ii[i] = idx;
    }

    /* --- Per-axis DDA parameters (Jacobs eq. 28-29) ------------------- */
    for (i = 0; i < 3; ++i) {
        if (d[i] > 0.0) {
            iu[i] = 1;
            au[i] = grid->spacing[i] / d[i];
            ai[i] = (grid->origin[i] + (ii[i] + 1) * grid->spacing[i] - p[i]) / d[i];
        } else if (d[i] < 0.0) {
            iu[i] = -1;
            au[i] = -grid->spacing[i] / d[i];
            ai[i] = (grid->origin[i] + ii[i] * grid->spacing[i] - p[i]) / d[i];
        } else {
            iu[i] = 0;
            au[i] = HUGE_VAL;
            ai[i] = HUGE_VAL;
        }
    }

    /* Pre-compute tile counts for Morton layout (unused for row-major).
     * Declared here (after statements) to keep them const and let the compiler
     * see them as invariant across the traversal loop — exception to the
     * project's C89-style top-of-function declaration convention. */
    size_t const Tx = (grid->n[0] + 7u) >> 3u;
    size_t const Ty = (grid->n[1] + 7u) >> 3u;

    ac = alpha_entry;
    n = 0;

    /* --- Walk through the grid (Jacobs eq. 30-34) --------------------- */
    for (;;) {
        /* Choose the axis whose next boundary crossing comes first. */
        if (ai[0] <= ai[1] && ai[0] <= ai[2])
            axis = 0;
        else if (ai[1] <= ai[2])
            axis = 1;
        else
            axis = 2;

        /* Clip to alpha_exit. */
        a_next = (ai[axis] < alpha_exit) ? ai[axis] : alpha_exit;

        crossings[n].idx = (grid->tile_order == OSH_VOXEL_ORDER_MORTON8)
                               ? osh_voxel_tile_idx((size_t) ii[0], (size_t) ii[1], (size_t) ii[2], Tx, Ty)
                               : (size_t) ii[0] + grid->n[0] * ((size_t) ii[1] + grid->n[1] * (size_t) ii[2]);
        crossings[n].path_len = (a_next - ac) * ds; /* convert alpha fraction to length */
        ++n;

        ac = a_next;
        if (ac >= alpha_exit)
            break;

        /* Advance every axis tied at the same crossing alpha. */
        for (axis = 0; axis < 3; ++axis) {
            if (!same_crossing(ai[axis], a_next))
                continue;

            ii[axis] += iu[axis];
            if (ii[axis] < 0 || ii[axis] >= (int) grid->n[axis]) {
                *n_out = n;
                return (n > 0) ? 1 : 0;
            }

            ai[axis] += au[axis]; /* Jacobs eq. 34 */
        }
    }

    *n_out = n;
    return (n > 0) ? 1 : 0;
}
