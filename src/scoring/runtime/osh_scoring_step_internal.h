#ifndef OSH_SCORING_STEP_INTERNAL_H
#define OSH_SCORING_STEP_INTERNAL_H

/*
 * Module-internal scoring declarations.
 *
 * NOT part of the public scoring API; do not include outside src/scoring.
 * Geometry helpers are shared by the step and point geometry drivers.
 * Estimator handlers live behind osh_scoring_estimator_internal.h.
 */

#include "common/osh_step.h"
#include "common/raytrace/osh_raytrace.h"
#include "openshieldhit/status.h"
#include "particle/osh_particle.h"
#include "scoring/runtime/osh_scoring_accumulator.h"
#include "scoring/runtime/osh_scoring_runtime.h"

/** Build a raytrace grid + 1/voxel-volume from a Mesh scoring geometry. */
enum osh_status mesh_geometry_to_grid(struct osh_scoring_geometry_runtime const *geo,
                                      struct osh_raytrace_grid *grid,
                                      double *voxel_volume_inv_out);

/**
 * Map a transport step's zone id to the dense Zone-scorer bin index.  Linear scan
 * over the (small) selected-zone list; returns 1 and sets @p *idx_out on a hit, 0
 * when @p zone is negative or is not one of this geometry's selected zones.
 */
int zone_bin_index(struct osh_scoring_geometry_runtime const *geo, int zone, size_t *idx_out);

/**
 * Build a raytrace grid from a Cyl (R,Z) scoring geometry.  Field convention:
 * origin/spacing/n[0] = r_min/dr/nr; origin/spacing/n[2] = z_min/dz/nz; index [1]
 * is unused (n[1] = 1).  Per-voxel 1/V is geo->cyl_vol_inv[r_bin] with
 * r_bin = flat_idx % nr and flat_idx = z_bin * nr + r_bin.
 */
enum osh_status cyl_geometry_to_grid(struct osh_scoring_geometry_runtime const *geo, struct osh_raytrace_grid *grid);

#endif /* OSH_SCORING_STEP_INTERNAL_H */
