#ifndef OSH_GEMCA_RUNTIME_VOXEL_H
#define OSH_GEMCA_RUNTIME_VOXEL_H

#include <stddef.h>

#include "common/osh_ray.h"
#include "gemca/runtime/osh_gemca_runtime.h"

/**
 * One voxel crossing within a transport step.
 *
 * ds  — path length through this voxel [cm].
 * rho — actual voxel density [g/cm³] from the HU calibration LUT.
 *       Constant over ds.  Areal density contribution: rho × ds [g/cm²].
 *       Never averaged across voxels.
 */
struct gemca_rt_voxel_segment {
    double ds;
    double rho;
};

/**
 * Recommended stack buffer size for transport callers.
 *
 * 256 × 16 bytes = 4 KB.  Covers a worst-case diagonal through a clinical
 * 512×512 CT (~210 voxels at 45°) in a single step.  If the actual traversal
 * exceeds this, dist_voxel_body_rt() stops at the buffer limit — the transport
 * re-enters naturally on the next step, just like any other step-length limit.
 */
#define OSH_GEMCA_VOXEL_SEGS_STACK 256

/**
 * Fill segs[] with per-voxel crossings for a voxel body and return the total
 * step distance.
 *
 * Traverses the voxel grid with Jacobs and stops at whichever comes first:
 *   - a material-bin change (new zone starts)
 *   - grid exit
 *   - segs_cap reached (acts as a step-length limiter; transport re-enters)
 *
 * Precondition: the compiled runtime must provide both `hu_bin_lut` and
 * `hu_rho_lut` for voxel transport. A debug build asserts this invariant.
 * Release builds return `OSH_GEMCA_INFINITY` defensively if the LUTs are
 * missing.
 *
 * segs and segs_cap may be NULL/0 for callers that only need the distance
 * (e.g. the CSG distance evaluator).  n_out and bin_out may also be NULL.
 *
 * Returns OSH_GEMCA_INFINITY if the ray misses the grid entirely.
 */
double dist_voxel_body_rt(struct osh_gemca_runtime const *rt,
                          int body_idx,
                          struct ray const *r,
                          struct gemca_rt_voxel_segment *segs, /* caller-owned; may be NULL */
                          size_t segs_cap,                     /* 0 if segs is NULL */
                          size_t *n_out, /* voxel crossings traversed; may exceed segs_cap; may be NULL */
                          int *bin_out); /* material bin for all segments; may be NULL */

#endif /* OSH_GEMCA_RUNTIME_VOXEL_H */
