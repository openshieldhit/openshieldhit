/*
 * Voxel traversal using the 3-D digital differential analyser (DDA).
 *
 * This is the pedagogical reference implementation, following the algorithm
 * described by Amanatides and Woo in "A Fast Voxel Traversal Algorithm for
 * Ray Tracing" (Eurographics, 1987).
 *
 * Design goals
 * ------------
 * - Exact: every plane crossing is found analytically; no step-size error.
 * - O(N): one iteration per voxel crossed, where N = n[0]+n[1]+n[2].
 * - Readable: control flow is linear, no precomputed auxiliary arrays.
 *
 * Limitation: the ray start point p must lie inside the grid.  Rays that
 * begin outside are not clipped to the grid entry; the function returns 0
 * crossings in that case.  The Siddon implementation (osh_raytrace_siddon_msh.c)
 * handles outside starts correctly by computing parametric entry/exit bounds
 * before the walk.
 *
 * Algorithm outline
 * -----------------
 * For each axis i, precompute:
 *   step[i]   - direction of travel (+1 or -1); 0 if v[i]==0 (no crossings).
 *   tmax[i]   - parametric t at which the ray first crosses a plane on axis i.
 *   tdelta[i] - parametric distance between consecutive plane crossings.
 *
 * At each iteration, choose the axis with the smallest tmax (nearest crossing),
 * record the current voxel with the path length from the previous crossing,
 * then advance one voxel along that axis.
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
    int vox[3];       /* current voxel index along each axis */
    int step[3];      /* +1 or -1 per axis; 0 when v[i]==0 */
    double tmax[3];   /* parametric t to next boundary crossing per axis */
    double tdelta[3]; /* parametric t between consecutive crossings per axis */
    double t;         /* parametric position along the ray */
    double t_next;    /* parametric position of the next recorded crossing */
    int axis;         /* axis chosen for the current step */
    size_t n;         /* number of crossings written so far */
    int i;

    /* Determine the starting voxel.  Return immediately if p is outside. */
    for (i = 0; i < 3; ++i) {
        vox[i] = (int) ((p[i] - grid->origin[i]) / grid->spacing[i]);
        if (vox[i] < 0 || vox[i] >= (int) grid->n[i]) {
            *n_out = 0;
            return 0;
        }
    }

    /*
     * For each axis, set up the DDA parameters.
     *
     * When v[i] > 0 the next boundary in the +x direction is at the high face
     * of the current voxel: origin[i] + (vox[i]+1)*spacing[i].
     * When v[i] < 0 the next boundary is at the low face: origin[i]+vox[i]*spacing[i].
     * When v[i] == 0 the ray is parallel to that axis; tmax and tdelta are
     * set to HUGE_VAL so this axis is never chosen as the minimum.
     */
    for (i = 0; i < 3; ++i) {
        if (v[i] > 0.0) {
            step[i] = 1;
            tdelta[i] = grid->spacing[i] / v[i];
            tmax[i] = (grid->origin[i] + (vox[i] + 1) * grid->spacing[i] - p[i]) / v[i];
        } else if (v[i] < 0.0) {
            step[i] = -1;
            tdelta[i] = -grid->spacing[i] / v[i];
            tmax[i] = (grid->origin[i] + vox[i] * grid->spacing[i] - p[i]) / v[i];
        } else {
            step[i] = 0;
            tdelta[i] = HUGE_VAL;
            tmax[i] = HUGE_VAL;
        }
    }

    t = 0.0;
    n = 0;

    for (;;) {
        /* Choose the axis whose next boundary crossing comes first. */
        if (tmax[0] <= tmax[1] && tmax[0] <= tmax[2])
            axis = 0;
        else if (tmax[1] <= tmax[2])
            axis = 1;
        else
            axis = 2;

        /* Clip the crossing to the end of the requested step length. */
        t_next = (tmax[axis] < ds) ? tmax[axis] : ds;

        crossings[n].idx = (size_t) vox[0] + grid->n[0] * ((size_t) vox[1] + grid->n[1] * (size_t) vox[2]);
        crossings[n].path_len = t_next - t;
        ++n;

        t = t_next;
        if (t >= ds)
            break;

        /* Advance every axis that hits the same plane crossing. This avoids
         * zero-length fake voxels when the ray crosses an edge or corner. */
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

        /* Stop if the ray exits the grid on any advanced axis. */
        if (vox[0] < 0 || vox[0] >= (int) grid->n[0] || vox[1] < 0 || vox[1] >= (int) grid->n[1] || vox[2] < 0
            || vox[2] >= (int) grid->n[2])
            break;
    }

    *n_out = n;
    return (n > 0) ? 1 : 0;
}
