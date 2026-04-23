#ifndef OSH_STEP_SEGMENT_H
#define OSH_STEP_SEGMENT_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * A step segment is one constant-density piece of a transport step.
 *
 * A transport step advances a particle from its current position to a new
 * position without crossing a material-zone boundary.  Inside an analytic
 * (non-voxelized) zone the step has exactly one segment — the zone is uniform.
 * Inside a CT voxel zone the step may contain many segments: each segment
 * corresponds to one voxel, all belonging to the same material bin, but
 * potentially with different densities.
 *
 * The step function (osh_transport_ion_step) receives the segment list
 * produced by the geometry layer before the physics is evaluated.  It uses:
 *   ds  — to accumulate the geometric step length
 *   rho — to accumulate the areal density Σ(rho_i × ds_i) for energy loss
 *           and straggling, and to weight per-segment dose deposition in scoring
 *
 * The segment list is built in the wavefront outer loop:
 *   - analytic zone: one synthetic segment {ds = zone_boundary_ds, rho = zone_rho}
 *   - CT voxel zone: one segment per voxel, filled by dist_voxel_body_rt()
 */
struct osh_step_segment {
    double ds;  /* path length through this constant-density piece [cm] */
    double rho; /* actual material density for this piece [g/cm³] */
};

/*
 * Stack-buffer capacity for step-segment arrays.
 *
 * 256 segments × 16 bytes = 4 KB per slot — fits comfortably on the stack
 * inside the per-slot transport loop.  Covers a worst-case diagonal through a
 * clinical 512×512 CT (~210 voxels at 45°) before a bin change.  When a step
 * hits this limit, dist_voxel_body_rt() stops at the capacity boundary and
 * the transport loop re-enters the same zone on the next step, naturally
 * limiting the step length.
 */
#define OSH_STEP_SEGMENTS_MAX 256

#ifdef __cplusplus
}
#endif

#endif /* OSH_STEP_SEGMENT_H */
