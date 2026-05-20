#include "gemca/runtime/osh_gemca_runtime_voxel.h"

#include <assert.h>
#include <stddef.h>

#include "common/osh_ray.h"
#include "common/raytrace/osh_raytrace.h"
#include "gemca/osh_gemca2.h"
#include "gemca/runtime/osh_gemca_runtime.h"

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
    struct osh_voxel_crossing crossing;
    struct ray r_local;
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
    if (!body->hu) {
        return OSH_GEMCA_INFINITY;
    }

    /* Transform ray into body-local frame where ct_grid is defined. */
    osh_ray_transform(r, &r_local, body->t);

    if (!osh_raytrace_first_crossing(&body->ct_grid, r_local.p, r_local.cp, 1e300, &crossing)) {
        return OSH_GEMCA_INFINITY;
    }

    /* Current M5 policy: one voxel is the current medium.  Stop at the first
     * voxel boundary even if the next voxel belongs to the same material bin. */
    hu_clamped = clamp_hu(body->hu[crossing.idx]);

    bin0 = (int) rt->hu_bin_lut[hu_clamped + 1000];

    if (step_segments && step_segments_cap > 0u) {
        step_segments[0].ds = crossing.path_len;
    }

    if (n_out) {
        *n_out = (!step_segments || step_segments_cap > 0u) ? 1u : 0u;
    }
    if (bin_out) {
        *bin_out = bin0;
    }
    return crossing.path_len;
}
