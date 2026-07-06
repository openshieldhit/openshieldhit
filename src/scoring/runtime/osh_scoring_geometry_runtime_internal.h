#ifndef OSH_SCORING_GEOMETRY_RUNTIME_INTERNAL_H
#define OSH_SCORING_GEOMETRY_RUNTIME_INTERNAL_H

/*
 * Module-internal scoring geometry helpers.
 *
 * NOT part of the public scoring API; do not include outside src/scoring.
 * These helpers are shared by the step and point geometry drivers.
 */

#include "common/osh_step.h"
#include "common/raytrace/osh_raytrace.h"
#include "openshieldhit/status.h"
#include "particle/osh_particle.h"
#include "scoring/runtime/osh_scoring_accumulator.h"
#include "scoring/runtime/osh_scoring_runtime.h"

/** Build a raytrace grid from a Mesh scoring geometry.  Per-bin volume is applied
 * separately in postprocess (geo->bin_vol_inv), not returned here. */
enum osh_status mesh_geometry_to_grid(struct osh_scoring_geometry_runtime const *geo, struct osh_raytrace_grid *grid);

/**
 * Map a transport step's zone id to the dense Zone-scorer bin index.  Linear scan
 * over the (small) selected-zone list; returns 1 and sets @p *idx_out on a hit, 0
 * when @p zone is negative or is not one of this geometry's selected zones.
 */
int zone_bin_index(struct osh_scoring_geometry_runtime const *geo, int zone, size_t *idx_out);

/**
 * Build a raytrace grid from a Cyl (R,Z) scoring geometry.  Field convention:
 * origin/spacing/n[0] = r_min/dr/nr; origin/spacing/n[2] = z_min/dz/nz; index [1]
 * is unused (n[1] = 1).  Flat bin index is z_bin * nr + r_bin; per-bin 1/V is
 * geo->bin_vol_inv (built at compile), applied in postprocess, not here.
 */
enum osh_status cyl_geometry_to_grid(struct osh_scoring_geometry_runtime const *geo, struct osh_raytrace_grid *grid);

#endif /* OSH_SCORING_GEOMETRY_RUNTIME_INTERNAL_H */
