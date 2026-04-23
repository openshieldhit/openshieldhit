#include "gemca/runtime/osh_gemca_runtime_voxel.h"

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
    if (hu < -1000)
        return -1000;
    if (hu > 1600)
        return 1600;
    return (int) hu;
}

double dist_voxel_body_rt(struct osh_gemca_runtime const *rt,
                          int body_idx,
                          struct ray const *r,
                          struct gemca_rt_voxel_segment *segs,
                          size_t segs_cap,
                          size_t *n_out,
                          int *bin_out) {
    struct gemca_rt_body const *body;
    struct osh_voxel_crossing crossings[_CROSSINGS_CAP];
    struct ray r_local;
    size_t n_crossings, i, n_segs;
    double total_ds;
    int bin0, bini, hu_clamped;

    if (n_out)
        *n_out = 0;
    if (bin_out)
        *bin_out = -1;

    if (!rt->hu_bin_lut || !rt->hu_rho_lut) {
        return OSH_GEMCA_INFINITY;
    }

    body = &rt->bodies[body_idx];

    /* Transform ray into body-local frame where ct_grid is defined. */
    osh_ray_transform(r, &r_local, body->t);

    if (!osh_raytrace_traverse(&body->ct_grid, r_local.p, r_local.cp, 1e300, crossings, &n_crossings)
        || n_crossings == 0) {
        return OSH_GEMCA_INFINITY;
    }

    /* Determine the starting bin from the first voxel. */
    hu_clamped = clamp_hu(body->hu[crossings[0].idx]);
    bin0 = (int) rt->hu_bin_lut[hu_clamped + 1000];

    /* Walk crossings: stop at bin change or segs_cap. */
    total_ds = 0.0;
    n_segs = 0;
    for (i = 0; i < n_crossings; i++) {
        hu_clamped = clamp_hu(body->hu[crossings[i].idx]);
        bini = (int) rt->hu_bin_lut[hu_clamped + 1000];

        if (bini != bin0)
            break; /* bin change → new zone */

        if (segs && n_segs < segs_cap) {
            segs[n_segs].ds = crossings[i].path_len;
            segs[n_segs].rho = (double) rt->hu_rho_lut[hu_clamped + 1000];
        }

        total_ds += crossings[i].path_len;
        n_segs++;

        if (segs && n_segs >= segs_cap)
            break; /* buffer full → transport re-enters */
    }

    if (n_out)
        *n_out = n_segs;
    if (bin_out)
        *bin_out = bin0;
    return total_ds;
}
