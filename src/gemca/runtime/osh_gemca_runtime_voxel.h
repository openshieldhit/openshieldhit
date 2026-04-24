#ifndef OSH_GEMCA_RUNTIME_VOXEL_H
#define OSH_GEMCA_RUNTIME_VOXEL_H

#include <stddef.h>

#include "common/osh_ray.h"
#include "common/osh_step_segment.h"
#include "gemca/runtime/osh_gemca_runtime.h"

/**
 * Fill step_segments[] with per-voxel crossings for a voxel body and return
 * the total step distance.
 *
 * Traverses the voxel grid with Jacobs (Siddon-like incremental DDA) and stops
 * at whichever comes first:
 *   - a material-bin change (the next voxel belongs to a different zone)
 *   - grid exit
 *   - step_segments_cap reached (acts as a step-length limiter; transport
 *     re-enters the same zone on the next call, just like any other limit)
 *
 * Each filled entry has:
 *   ds  — path length through that voxel [cm]
 *   rho — 0.0 (density is resolved later by the material runtime via
 *          osh_material_runtime_get_rho(), which uses its own HU→density LUT)
 *
 * Precondition: the compiled runtime must provide `hu_bin_lut` for voxel
 * transport.  A debug build asserts this invariant.  Release builds return
 * `OSH_GEMCA_INFINITY` defensively if the LUT is missing.
 *
 * step_segments and step_segments_cap may be NULL/0 for callers that only need
 * the distance (e.g. the CSG distance evaluator).  n_out and bin_out may also
 * be NULL.
 *
 * Returns OSH_GEMCA_INFINITY if the ray misses the grid entirely.
 */
double dist_voxel_body_rt(struct osh_gemca_runtime const *rt,
                          int body_idx,
                          struct ray const *r,
                          struct osh_step_segment *step_segments, /* caller-owned; may be NULL */
                          size_t step_segments_cap,               /* 0 if step_segments is NULL */
                          size_t *n_out, /* segments filled (≤ step_segments_cap); may be NULL */
                          int *bin_out); /* material bin for all segments; may be NULL */

#endif /* OSH_GEMCA_RUNTIME_VOXEL_H */
