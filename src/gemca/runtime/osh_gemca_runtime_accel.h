#ifndef OSH_GEMCA_RUNTIME_ACCEL_H
#define OSH_GEMCA_RUNTIME_ACCEL_H

#include <stddef.h>
#include <stdint.h>

#include "gemca/runtime/osh_gemca_runtime.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @defgroup osh_gemca_accel Geometry acceleration structure (prototype)
 *
 * @brief Uniform-grid spatial index over the compiled GEMCA runtime.
 *
 * @details
 * Built once by osh_gemca_accel_build() after osh_gemca_compile() has
 * produced the flat runtime arrays.  Three pieces of derived data:
 *
 *   1. Per-body axis-aligned bounding boxes (AABBs) in the universe frame,
 *      derived from each body's surface set and coordinate transform.
 *      Conservatively inflated; bodies whose extent cannot be bounded
 *      (e.g. infinite planes) are flagged unbounded and never culled.
 *
 *   2. A uniform grid over the union of all bounded zone AABBs.  Each cell
 *      stores a candidate list of zones whose AABB overlaps the cell.  For
 *      every (cell, zone) pair the zone's RPN program is *specialized* by
 *      constant folding: PUSH instructions for bodies whose AABB does not
 *      touch the cell are folded to constant "outside", and the boolean
 *      operators are simplified accordingly.  A zone whose program folds to
 *      constant false is dropped from the cell.  Zones with unbounded AABBs
 *      are kept in a global "always" list evaluated with their full program.
 *
 *   3. Per-zone flat lists of unique leaf bodies, used by the boundary
 *      distance query: the distance to the next zone boundary equals the
 *      minimum over leaf-body boundary distances (all CSG operators combine
 *      distances with minpos), so the RPN walk can be replaced by a flat
 *      minimum with ray/AABB slab culling.
 *
 * The accelerator is an optional, semantics-preserving cache: every query
 * answered through it returns the same result as the linear-scan fallback.
 * Building it can fail (allocation, pathological geometry); callers must
 * treat a NULL accel pointer as "use the fallback path".
 *
 * @{
 */

/* A zone whose distance evaluation touches fewer leaf bodies than this keeps
 * the plain RPN walk: grid bookkeeping costs more than it saves there. */
#define OSH_GEMCA_ACCEL_DIST_MIN_LEAVES 12

/* Geometries smaller than these thresholds skip the accelerator entirely:
 * a linear scan over a handful of small zones beats any spatial index. */
#define OSH_GEMCA_ACCEL_MIN_ZONES 16
#define OSH_GEMCA_ACCEL_MIN_ZONE_INSNS 25

/** Axis-aligned bounding box; unbounded extents use +/-HUGE_VAL. */
struct gemca_rt_aabb {
    double lo[3];
    double hi[3];
};

/** One candidate zone in a grid cell, with its cell-specialized program. */
struct gemca_accel_cand {
    uint32_t zone_idx;   /**< Dense zone index (ascending within a cell). */
    uint32_t insn_begin; /**< Offset of the specialized program in insns[]. */
    uint32_t ninsns;     /**< Length of the specialized program. */
};

struct osh_gemca_accel {
    /* Grid geometry */
    double lo[3];       /**< Domain minimum corner. */
    double hi[3];       /**< Domain maximum corner. */
    double inv_cell[3]; /**< 1 / cell size per axis. */
    int n[3];           /**< Cell counts per axis. */
    size_t ncells;      /**< n[0]*n[1]*n[2]. */

    /* Per-cell candidate zones (CSR layout, cell-major) */
    uint32_t *cand_begin;           /**< ncells+1 offsets into cands[]. */
    struct gemca_accel_cand *cands; /**< Candidate pool. */
    struct gemca_rt_insn *insns;    /**< Specialized instruction pool. */

    /* Zones whose AABB is unbounded: always evaluated with the full program */
    uint32_t *always; /**< Ascending zone indices. */
    size_t nalways;

    /* Per-body universe-frame AABBs (conservatively inflated) */
    struct gemca_rt_aabb *body_aabb; /**< nbodies entries. */
    uint8_t *body_bounded;           /**< nbodies entries; 1 = finite on all axes. */
    uint8_t *zone_bounded;           /**< nzones entries; 1 = zone AABB finite (zone is in the grid). */

    /* Per-zone unique leaf bodies for the flat distance evaluation.  Within a
     * zone's slice, unbounded bodies come first (they can never be culled and
     * the grid walk evaluates them once, upfront). */
    uint32_t *zone_leaf_begin;     /**< nzones+1 offsets into zone_leaves[]. */
    uint32_t *zone_leaves;         /**< Leaf body index pool. */
    uint32_t *zone_leaf_unbounded; /**< nzones entries: count of leading unbounded leaves. */
};

/**
 * @brief Build the acceleration structure for a compiled runtime.
 *
 * @details
 * On success rt->accel points to a fully built accelerator.  On any failure
 * (allocation, no bounded zones, size caps exceeded) rt->accel is left NULL
 * and OSH_OK is still returned: the accelerator is an optional optimization
 * and the runtime falls back to linear scans.  Set the environment variable
 * OSH_GEMCA_ACCEL=0 to skip building it (for A/B benchmarking).
 *
 * @param[in,out] rt  Compiled runtime (surfaces/bodies/zones populated).
 *
 * @returns OSH_OK (failures degrade to the fallback path), OSH_EINVAL if rt
 *          is NULL.
 */
enum osh_status osh_gemca_accel_build(struct osh_gemca_runtime *rt);

/**
 * @brief Release an acceleration structure.
 *
 * @param[in] accel  Accelerator to free; may be NULL.
 */
void osh_gemca_accel_free(struct osh_gemca_accel *accel);

/** @} */ /* end defgroup osh_gemca_accel */

#ifdef __cplusplus
}
#endif

#endif /* OSH_GEMCA_RUNTIME_ACCEL_H */
