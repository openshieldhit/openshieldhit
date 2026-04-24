#include "gemca/runtime/osh_gemca_runtime_voxel.h"

#include <assert.h>
#include <stddef.h>

#include "common/osh_ray.h"
#include "common/raytrace/osh_raytrace.h"
#include "gemca/osh_gemca2.h"
#include "gemca/runtime/osh_gemca_runtime.h"

/*
 * Safe upper bound on crossings for clinical CT grids (≤ 512×512×512).
 * A ray can cross at most nx+ny+nz voxels: 512+512+512 = 1536 < 2048.
 * Increase if ever using larger grids.
 */
#define _CROSSINGS_CAP 2048

static int clamp_hu(int16_t hu) {
    if (hu < -1000) {
        return -1000;
    }
    if (hu > 1600) {
        return 1600;
    }
    return (int) hu;
}

double dist_voxel_body_rt(struct osh_gemca_runtime const *rt,
                          int body_idx,
                          struct ray const *r,
                          struct osh_step_segment *step_segments,
                          size_t step_segments_cap,
                          size_t *n_out,
                          int *bin_out) {
    struct gemca_rt_body const *body;
    struct osh_voxel_crossing crossings[_CROSSINGS_CAP];
    struct ray r_local;
    size_t crossings_needed;
    size_t n_crossings;
    int bin0;
    int hu_clamped;

    if (n_out) {
        *n_out = 0;
    }
    if (bin_out) {
        *bin_out = -1;
    }

    assert(rt->hu_bin_lut);
    if (!rt->hu_bin_lut) {
        return OSH_GEMCA_INFINITY;
    }

    body = &rt->bodies[body_idx];

    /* osh_raytrace_traverse() requires capacity >= n[0] + n[1] + n[2]. */
    crossings_needed = body->ct_grid.n[0];
    if (crossings_needed > _CROSSINGS_CAP) {
        return OSH_GEMCA_INFINITY;
    }
    if (body->ct_grid.n[1] > _CROSSINGS_CAP - crossings_needed) {
        return OSH_GEMCA_INFINITY;
    }
    crossings_needed += body->ct_grid.n[1];
    if (body->ct_grid.n[2] > _CROSSINGS_CAP - crossings_needed) {
        return OSH_GEMCA_INFINITY;
    }

    /* Transform ray into body-local frame where ct_grid is defined. */
    osh_ray_transform(r, &r_local, body->t);

    if (!osh_raytrace_traverse(&body->ct_grid, r_local.p, r_local.cp, 1e300, crossings, &n_crossings)
        || n_crossings == 0) {
        return OSH_GEMCA_INFINITY;
    }

    if (!body->hu) {
        return OSH_GEMCA_INFINITY;
    }

    /* Current M5 policy: one voxel is the current medium.  Stop at the first
     * voxel boundary even if the next voxel belongs to the same material bin. */
    hu_clamped = clamp_hu(body->hu[crossings[0].idx]);
    bin0 = (int) rt->hu_bin_lut[hu_clamped + 1000];

    if (step_segments && step_segments_cap > 0u) {
        step_segments[0].ds = crossings[0].path_len;
        step_segments[0].rho = 0.0;
    }

    if (n_out) {
        *n_out = (!step_segments || step_segments_cap > 0u) ? 1u : 0u;
    }
    if (bin_out) {
        *bin_out = bin0;
    }
    return crossings[0].path_len;
}
