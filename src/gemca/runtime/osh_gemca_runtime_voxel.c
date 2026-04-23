#include "gemca/runtime/osh_gemca_runtime_voxel.h"

#include <stddef.h>

#include "common/osh_ray.h"
#include "gemca/runtime/osh_gemca_runtime.h"

/*
 * TODO M4: replace this stub with Jacobs voxel traversal.
 *
 *  1. Transform ray to body-local frame (reuse transform_to_local_rt pattern).
 *  2. Call osh_raytrace_traverse() on body->ct_grid with a local crossings
 *     buffer (stack-allocated, size OSH_GEMCA_VOXEL_SEGS_STACK is sufficient).
 *  3. If n_crossings == 0 (ray misses grid): return OSH_GEMCA_INFINITY.
 *  4. Record the starting bin:
 *       int bin0 = osh_gemca_voxel_hu2idx(hu[crossings[0].idx]);
 *  5. Walk crossings up to min(n_crossings, segs_cap):
 *       double total_ds = 0.0;
 *       for each crossing i:
 *           int bin_i = osh_gemca_voxel_hu2idx(hu[crossings[i].idx]);
 *           if (bin_i != bin0) break;          // bin change — stop, new zone
 *           double rhocorr_i = osh_gemca_voxel_hu2rho(hu[crossings[i].idx])
 *                              / rho_bin_nominal; // per-voxel, never averaged
 *           if (segs && i < segs_cap):
 *               segs[i].ds      = crossings[i].path_len;
 *               segs[i].rhocorr = rhocorr_i;    // exact for this voxel
 *           total_ds += crossings[i].path_len;
 *           // segs_cap reached → stop; transport re-enters on next step
 *  6. Write n_out, bin_out if non-NULL.
 *  7. Return total_ds.
 *
 * Scoring callers walk segs[0..n_out-1] and may sub-step within each segment
 * against their own scoring grid.  rhocorr is guaranteed constant over ds for
 * each segment, so dose weighting (rhocorr_i × sub_ds_j) is always exact.
 */
double dist_voxel_body_rt(struct osh_gemca_runtime const *rt,
                          int body_idx,
                          struct ray const *r,
                          struct gemca_rt_voxel_segment *segs,
                          size_t segs_cap,
                          size_t *n_out,
                          int *bin_out) {
    (void) rt;
    (void) body_idx;
    (void) r;
    (void) segs;
    (void) segs_cap;
    if (n_out)
        *n_out = 0;
    if (bin_out)
        *bin_out = -1;
    return OSH_GEMCA_INFINITY;
}
