#ifndef OSH_GEMCA_RUNTIME_H
#define OSH_GEMCA_RUNTIME_H

#include <stddef.h>

#include "openshieldhit/geometry.h"
#include "openshieldhit/status.h"

/* Forward declaration: full definition is in common/osh_ray.h, included by
 * consumers that need struct ray fields.  The header only uses struct ray as a
 * pointer parameter, so a forward declaration suffices here. */
struct ray;
#include "gemca/osh_gemca2.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @defgroup osh_gemca_runtime Geometry Runtime (GEMCA runtime layer)
 *
 * @brief Compiled, cache-friendly representation of GEMCA geometry for use in
 *        the hot transport kernel.
 *
 * @details
 * The cold @ref osh_gemca_prepared (produced by osh_gemca_load()) stores geometry
 * as a set of pointer-linked structs: zones reference cgnode ASTs, cgnodes
 * reference bodies, bodies reference dynamically-allocated surface arrays.
 * This layout is convenient for parsing but causes pointer-chasing cache misses
 * in the transport inner loop, which calls get_zone and get_distance once per
 * step per particle.
 *
 * The runtime layer compiles the cold workspace once into three flat contiguous
 * arrays (surfaces, bodies, zones) and replaces the recursive AST traversal
 * with a flat RPN (reverse Polish notation) instruction evaluator.  All
 * transport-relevant geometry data is reachable without heap indirection.
 *
 * Zone evaluation uses a small integer stack whose depth is bounded by the
 * depth of the CSG tree (not by the number of bodies).  Each zone's instruction
 * array is prepended with one or more @ref GEMCA_RT_GUARD_BODY instructions
 * auto-derived at setup time; a guard rejects the ray for that zone in O(1)
 * surface evaluations before the full CSG expression is checked.
 *
 * Voxel bodies are represented as bounding RPP surfaces for zone-membership
 * purposes (@ref GEMCA_RT_PUSH_VOXEL_BODY).  A zone that contains a voxel body
 * records the body index in @ref gemca_rt_zone.voxel_body_idx; the distance
 * query dispatches to the voxel traversal helper and returns the current voxel
 * exit distance.
 *
 * @{
 */

/* ---- RPN opcodes --------------------------------------------------------- */

/**
 * @brief Opcodes for the flat RPN instruction array of a compiled zone.
 *
 * @details
 * A zone's CSG expression tree is compiled into a flat post-order sequence of
 * instructions at setup time.  Evaluation uses a small value stack:
 *   - For zone membership: the stack holds int (0/1) values.
 *   - For boundary distance: the stack holds (double dist, int is_inside) pairs.
 * Both evaluations consume the same instruction array; only the action per
 * opcode differs.
 *
 * Instruction semantics:
 *
 *   GUARD_BODY(idx)     — Before touching the value stack, check whether the
 *                         ray is inside body[idx].  If not, return 0 (outside
 *                         zone) immediately.  Multiple GUARD instructions may
 *                         appear at the front of the array.  In the distance
 *                         evaluator, GUARD instructions are no-ops (zone
 *                         membership is already confirmed before distance is
 *                         queried).
 *
 *   PUSH_BODY(idx)      — Evaluate body[idx] membership (and distance for the
 *                         distance evaluator); push the result.
 *
 *   PUSH_VOXEL_BODY(idx) — Identical to PUSH_BODY for membership.  For
 *                         distance queries, dispatches to the voxel traversal
 *                         helper and returns the current voxel-exit distance.
 *
 *   UNION               — Pop a, b (a=left, b=right); push (a || b).
 *   INTERSECT           — Pop a, b; push (a && b).
 *   DIFF                — Pop a, b; push (a && !b)  (left minus right).
 */
enum gemca_rt_op {
    GEMCA_RT_GUARD_BODY,      /**< Fast-reject guard: early exit if outside body[operand]. */
    GEMCA_RT_PUSH_BODY,       /**< Push body[operand] membership/distance result. */
    GEMCA_RT_PUSH_VOXEL_BODY, /**< Push voxel body[operand]; distance queries use voxel traversal. */
    GEMCA_RT_UNION,           /**< Boolean OR of top two stack entries. */
    GEMCA_RT_INTERSECT,       /**< Boolean AND of top two stack entries. */
    GEMCA_RT_DIFF             /**< Boolean (left && !right) of top two stack entries. */
};

/* ---- Compile-time limits ------------------------------------------------- */

/**
 * @brief Fixed parameter array size for flat surfaces.
 *
 * @details
 * The widest surface type (arbitrary plane) requires four parameters
 * [A, B, C, D] for Ax + By + Cz + D = 0.  All other surface types use at most
 * three.  Unused slots are zero-padded at setup time.  Fixed-size storage
 * eliminates the `double *p` indirection present in the cold @ref surface
 * struct and keeps surface data contiguous in the flat array.
 */
#define OSH_GEMCA_RT_SURF_NPAR 4

/**
 * @brief Maximum depth of the RPN evaluation stack.
 *
 * @details
 * Stack depth is bounded by the depth of the CSG tree, not by the number of
 * bodies.  A perfectly balanced binary tree of 256 leaf bodies requires at most
 * ceil(log2(256)) + 1 = 9 stack slots.  32 is a generous upper bound for any
 * geometry expected in practice.
 */
#define OSH_GEMCA_RT_MAX_STACK 32

/* ---- Flat runtime types -------------------------------------------------- */

/**
 * @brief A surface in the flat runtime representation.
 *
 * @details
 * Replaces the dynamically allocated cold @ref surface struct (with its
 * `double *p` indirection) with a fixed-size layout.  The `type` field uses
 * the same OSH_GEMCA_SURF_* identifiers as the cold representation; parameter
 * semantics are identical.  Unused parameter slots are zero-padded.
 */
struct gemca_rt_surface {
    double p[OSH_GEMCA_RT_SURF_NPAR]; /**< Surface parameters, zero-padded for unused slots. */
    int type;                         /**< Surface type: OSH_GEMCA_SURF_* */
};

/**
 * @brief A body in the flat runtime representation.
 *
 * @details
 * Replaces the cold @ref body struct and its `struct surface **surfs` double
 * indirection with a direct offset into the runtime's flat surface array.
 * Parse-time fields (name, lineno, a[], tokens) are not carried over.
 *
 * Surfaces belonging to this body occupy `osh_gemca_runtime.surfaces[surf_begin ..
 * surf_begin + nsurfs)` as a contiguous slice.
 *
 * For voxel bodies, ct_grid stores the local grid descriptor used by raytrace
 * traversal and hu is a borrowed pointer to the immutable HU volume (owned by
 * app/cold CT storage, not by the runtime).
 */
struct gemca_rt_body {
    double t[16];                     /**< 4x4 row-major transformation matrix (universe → body-local). */
    struct osh_raytrace_grid ct_grid; /**< Local voxel grid descriptor for VOX bodies. */
    int16_t const *hu;                /**< Borrowed HU pointer for VOX bodies; may be NULL for non-VOX bodies. */
    size_t surf_begin;                /**< Start index into osh_gemca_runtime.surfaces[]. */
    int nsurfs;                       /**< Number of surfaces in surfaces[surf_begin..surf_begin+nsurfs). */
    int type;                         /**< Body type: OSH_GEMCA_BODY_* (used to detect voxel bodies). */
    char coord;                       /**< Coordinate system: OSH_COORD_* value. */
};

/**
 * @brief A single RPN instruction in a compiled zone expression.
 *
 * @details
 * For GUARD_BODY, PUSH_BODY, and PUSH_VOXEL_BODY: `operand` is the 0-based
 * index into osh_gemca_runtime.bodies[].  For operator instructions (UNION,
 * INTERSECT, DIFF): `operand` is unused and is set to -1 at setup time.
 */
struct gemca_rt_insn {
    int op;      /**< Opcode: one of enum gemca_rt_op. */
    int operand; /**< Body index for push/guard opcodes; -1 for operators. */
};

/**
 * @brief A compiled zone in the flat runtime representation.
 *
 * @details
 * Holds the flat RPN instruction array compiled from the zone's cgnode AST
 * during osh_gemca_compile().  The instruction array is owned by this
 * struct and freed by osh_gemca_runtime_free().
 *
 * GUARD_BODY instructions (if any) appear at the front of insns[].  They are
 * followed by the main RPN CSG expression.  The evaluator processes the array
 * left-to-right; a failed GUARD causes an immediate early return of 0 (outside
 * zone) without touching the value stack.
 *
 * If `voxel_body_idx >= 0`, the zone contains a voxel body.  Zone membership
 * treats the voxel's surfaces as a regular body (RPP half-planes).  Distance
 * queries for this zone return the current voxel-exit distance.
 *
 * @note GPU migration: `insns` is a host heap pointer and cannot be followed
 * on a GPU device.  Before writing a GPU kernel, add a flat `insns_flat[]`
 * array and a per-zone `insn_begin[]` offset array to `osh_gemca_runtime` (see
 * runtime/README.md for the full plan).  The CPU path is unaffected.
 */
struct gemca_rt_zone {
    struct gemca_rt_insn *insns; /**< Flat RPN instruction array; owned. */
    size_t material_idx;         /**< Dense material index (resolved from cold zone). */
    int ninsns;                  /**< Number of entries in insns[]. */
    int voxel_body_idx;          /**< Index of voxel body in bodies[], or -1 if none. */
};

enum osh_gemca_zone_batch_dispatch {
    OSH_GEMCA_ZONE_BATCH_DISPATCH_SCALAR = 0,
    OSH_GEMCA_ZONE_BATCH_DISPATCH_SCALAR_NOCPU = 1,
    OSH_GEMCA_ZONE_BATCH_DISPATCH_AVX2 = 2
};

/**
 * @brief Compiled runtime representation of the full GEMCA geometry.
 *
 * @details
 * Built once from a cold @ref osh_gemca_prepared by osh_gemca_compile() and
 * passed to the transport kernel for the duration of a run.  The cold workspace
 * is not owned and must outlive this struct.
 *
 * Three contiguous flat arrays replace the cold workspace's pointer-linked
 * layout:
 *
 *   surfaces[] — all surfaces for all bodies concatenated.  Body i owns
 *                surfaces[bodies[i].surf_begin .. +bodies[i].nsurfs).
 *
 *   bodies[]   — one entry per body: transform matrix, coordinate system tag,
 *                and offset into surfaces[].
 *
 *   zones[]    — one entry per zone: compiled RPN instruction array, resolved
 *                material index, and optional voxel body index.
 *
 * Ownership: surfaces[], bodies[], and zones[] (including each zone's insns[])
 * are heap-allocated by setup and freed by osh_gemca_runtime_free().
 *
 * @note GPU migration: surfaces[], bodies[], and zones[] are plain contiguous
 * buffers suitable for device-memory copies as-is.  The only host-specific
 * field is zones[j].insns (a heap pointer).  See runtime/README.md for the
 * planned `insns_flat` / `insn_begin` addition that closes this gap.
 */
struct osh_gemca_runtime {
    struct osh_gemca_prepared const *workspace; /**< Cold storage reference — not owned. */
    struct gemca_rt_surface *surfaces;          /**< Flat surface array (owned). */
    struct gemca_rt_body *bodies;               /**< Flat body array (owned). */
    struct gemca_rt_zone *zones;                /**< Flat zone array; each zone owns its insns[]. */
    uint8_t *hu_bin_lut;     /**< HU→bin index, 2601 entries indexed by hu+1000; NULL for non-VOX runs; owned. */
    size_t nsurfaces;        /**< Number of entries in surfaces[]. */
    size_t nbodies;          /**< Number of entries in bodies[]. */
    size_t nzones;           /**< Number of entries in zones[]. */
    int zone_batch_dispatch; /**< enum osh_gemca_zone_batch_dispatch */
};

/* ---- Lifecycle ------------------------------------------------------------ */

/**
 * @brief Compile a cold gemca workspace into a flat runtime representation.
 *
 * @details
 * Walks the cold @ref osh_gemca_prepared and produces three flat arrays:
 *
 *   1. surfaces[] — all body surfaces copied into a contiguous block with
 *      fixed-size parameter arrays (no `double *p` indirection).
 *
 *   2. bodies[] — per-body transform matrix, coordinate system, surface
 *      slice offset, and body type.  Parse-time fields (names, tokens) are
 *      discarded.
 *
 *   3. zones[] — per-zone RPN instruction array compiled from the cgnode AST
 *      via post-order traversal, with GUARD_BODY instructions auto-prepended
 *      where a suitable guard body can be identified.  The resolved material
 *      index is copied from the cold zone.
 *
 * Zone material indices must be fully resolved before this function is called.
 * Passing unresolved indices produces a runtime whose material lookups will be
 * incorrect.
 *
 * @param[in]  wg             Fully loaded and material-resolved cold workspace.
 * @param[in]  hu_table_type  OSH_HU_TABLE_* selector used to build the optional
 *                            voxel HU→material-bin LUT. Pass
 *                            OSH_HU_TABLE_NONE for non-CT geometry.
 * @param[out] rt             Receives the compiled runtime.  The caller owns
 *                            the inner allocations; call
 *                            osh_gemca_runtime_free() to release them. The
 *                            struct may be stack-allocated; zero-initialise it
 *                            before calling this function.
 *
 * @returns OSH_OK on success, OSH_EINVAL if wg or rt is NULL, OSH_ENOMEM on
 *          allocation failure.
 */
enum osh_status osh_gemca_compile(struct osh_gemca_prepared const *wg,
                                  int hu_table_type,
                                  struct osh_diag_sink const *diag,
                                  struct osh_gemca_runtime *rt);
char const *osh_gemca_runtime_zone_batch_dispatch_name(struct osh_gemca_runtime const *rt);

/**
 * @brief Release all allocations owned by a gemca runtime.
 *
 * @details
 * Frees surfaces[], bodies[], and zones[] including per-zone insns[] buffers.
 * The cold workspace pointer is not freed.  Safe to call on a zero-initialised
 * struct.
 *
 * @param[in] rt  Runtime to release; may be NULL.
 */
void osh_gemca_runtime_free(struct osh_gemca_runtime *rt);

/* ---- Hot query API ------------------------------------------------------- */

/**
 * @brief Return the dense zone index for the zone containing a ray.
 *
 * @details
 * Iterates over zones in order.  For each zone the RPN instruction array is
 * evaluated left-to-right.  GUARD_BODY instructions at the front of the array
 * provide fast rejection: if the ray is outside the guard body, the zone is
 * skipped without evaluating the full CSG expression.
 *
 * When a zone's RPN expression evaluates to 1 (ray inside), its index is
 * returned immediately.
 *
 * @param[in] rt  Compiled gemca runtime.
 * @param[in] r   Ray whose position is tested (must be in OSH_COORD_UNIVERSE).
 *
 * @returns 0-based zone index, or OSH_GEMCA_ZONE_INDEX_INVALID if the ray is
 *          outside all defined zones.
 */
size_t osh_gemca_runtime_get_zone(struct osh_gemca_runtime const *rt, struct ray const *r);

/**
 * @brief Return the distance a ray travels inside zone `zone_idx` before leaving it.
 *
 * @details
 * Evaluates the RPN distance expression at the current ray position to obtain
 * (d, is_inside).  While is_inside is true, advances the ray by d and
 * accumulates the total path length.  This mirrors the step-loop in the cold
 * osh_gemca_get_distance() implementation.
 *
 * For zones with a voxel body (voxel_body_idx >= 0), this returns the distance
 * to the current voxel exit.  The current M5 transport policy treats each voxel
 * as the active medium for one step and re-queries after crossing a voxel
 * boundary.
 *
 * @param[in] rt        Compiled gemca runtime.
 * @param[in] zone_idx  Zone index returned by osh_gemca_runtime_get_zone().
 * @param[in] r         Ray (position and direction; direction will be normalised
 *                      on a local copy — the caller's ray is not modified).
 *
 * @returns Total path length inside the zone, in the same units as the geometry.
 */
double osh_gemca_runtime_get_distance(struct osh_gemca_runtime const *rt, size_t zone_idx, struct ray const *r);

/* ---- Batch query API ----------------------------------------------------- */

/**
 * @brief Surface-membership query for a batch of @p n rays already expressed in
 *        a body's local coordinate system.
 *
 * @details
 * Evaluates the same inside/on-boundary predicate as the scalar runtime
 * surface check used by zone membership.  Inputs are structure-of-arrays so
 * future SIMD implementations can operate directly on contiguous coordinates.
 *
 * Boundary semantics are identical to the scalar path:
 *   - negative implicit-function value => inside
 *   - positive implicit-function value => outside
 *   - on-surface cases are resolved by the ray direction
 *
 * This is a low-level primitive helper intended for batched body and zone
 * evaluators.  The input coordinates must already be in the surface's local
 * body coordinate system; no body transform is applied here.
 *
 * @param[in]  sf         Flat surface to test.
 * @param[in]  x          Ray x-positions, length @p n.
 * @param[in]  y          Ray y-positions, length @p n.
 * @param[in]  z          Ray z-positions, length @p n.
 * @param[in]  ux         Ray x-direction components, length @p n.
 * @param[in]  uy         Ray y-direction components, length @p n.
 * @param[in]  uz         Ray z-direction components, length @p n.
 * @param[in]  n          Number of rays.
 * @param[out] inside_out Per-ray 0/1 inside results, length @p n.
 */
void osh_gemca_runtime_check_surface_batch(struct gemca_rt_surface const *sf,
                                           double const *x,
                                           double const *y,
                                           double const *z,
                                           double const *ux,
                                           double const *uy,
                                           double const *uz,
                                           size_t n,
                                           int *inside_out);

/**
 * @brief Body-membership query for a batch of @p n rays in OSH_COORD_UNIVERSE.
 *
 * @details
 * Applies the body's coordinate transform to the input SoA rays, then checks
 * all surfaces in the body's contiguous surface slice.  The result for each
 * ray is the logical AND of all surface-inside predicates for that body.
 *
 * This is the batched counterpart of the scalar runtime body-membership check
 * used by zone evaluation.  It is intended as a reusable primitive for future
 * batched zone evaluators.
 *
 * @param[in]  rt         Compiled gemca runtime.
 * @param[in]  body_idx   Index into rt->bodies[].
 * @param[in]  x          Ray x-positions in OSH_COORD_UNIVERSE, length @p n.
 * @param[in]  y          Ray y-positions in OSH_COORD_UNIVERSE, length @p n.
 * @param[in]  z          Ray z-positions in OSH_COORD_UNIVERSE, length @p n.
 * @param[in]  ux         Ray x-direction components, length @p n.
 * @param[in]  uy         Ray y-direction components, length @p n.
 * @param[in]  uz         Ray z-direction components, length @p n.
 * @param[in]  n          Number of rays.
 * @param[out] inside_out Per-ray 0/1 inside results, length @p n.
 */
void osh_gemca_runtime_check_body_batch(struct osh_gemca_runtime const *rt,
                                        size_t body_idx,
                                        double const *x,
                                        double const *y,
                                        double const *z,
                                        double const *ux,
                                        double const *uy,
                                        double const *uz,
                                        size_t n,
                                        int *inside_out);

/**
 * @brief Zone lookup for a batch of @p n particles.
 *
 * @details
 * Equivalent to calling osh_gemca_runtime_get_zone() once per particle, but
 * expressed over SoA position/direction arrays so the compiler and CPU
 * prefetcher can better pipeline the zone iterations.
 *
 * Each output element is written unconditionally; particles outside all
 * defined zones receive @ref OSH_GEMCA_ZONE_INDEX_INVALID.
 *
 * @param[in]  rt       Compiled gemca runtime.
 * @param[in]  x        Particle x-positions, length @p n.
 * @param[in]  y        Particle y-positions, length @p n.
 * @param[in]  z        Particle z-positions, length @p n.
 * @param[in]  ux       Particle x-directions, length @p n.
 * @param[in]  uy       Particle y-directions, length @p n.
 * @param[in]  uz       Particle z-directions, length @p n.
 * @param[in]  n        Number of particles.
 * @param[out] zone_out Receives the zone index for each particle, length @p n.
 */
void osh_gemca_runtime_get_zone_batch(struct osh_gemca_runtime const *rt,
                                      double const *x,
                                      double const *y,
                                      double const *z,
                                      double const *ux,
                                      double const *uy,
                                      double const *uz,
                                      size_t n,
                                      size_t *zone_out);

/**
 * @brief Boundary-distance query for a batch of @p n particles.
 *
 * @details
 * Equivalent to calling osh_gemca_runtime_get_distance() once per particle.
 * Particles whose zone_ref.zone_idx is @ref OSH_ZONE_INDEX_INVALID receive
 * a distance of 0.0 and are skipped without touching the geometry.
 *
 * @param[in]  rt        Compiled gemca runtime.
 * @param[in]  x         Particle x-positions, length @p n.
 * @param[in]  y         Particle y-positions, length @p n.
 * @param[in]  z         Particle z-positions, length @p n.
 * @param[in]  ux        Particle x-directions, length @p n.
 * @param[in]  uy        Particle y-directions, length @p n.
 * @param[in]  uz        Particle z-directions, length @p n.
 * @param[in]  zone_refs Zone references from osh_gemca_runtime_get_zone_ref_batch(), length @p n.
 * @param[in]  n         Number of particles.
 * @param[out] dist_out  Receives the boundary distance for each particle, length @p n.
 */
void osh_gemca_runtime_get_distance_batch(struct osh_gemca_runtime const *rt,
                                          double const *x,
                                          double const *y,
                                          double const *z,
                                          double const *ux,
                                          double const *uy,
                                          double const *uz,
                                          struct osh_zone_ref const *zone_refs,
                                          size_t n,
                                          double *dist_out);

/**
 * @brief Zone-reference query for a batch of @p n particles.
 *
 * @details
 * Combines zone lookup with current-voxel HU and material-index resolution into a
 * single pass.  For analytic zones, @p has_hu is 0 and @p material_idx is the
 * zone's ASSIGNMAT index.  For voxel zones, @p has_hu is 1, @p hu holds the
 * clamped HU value of the current voxel, and @p material_idx is the HU-bin
 * index from the HU→bin LUT.
 *
 * Particles outside all zones receive @ref OSH_ZONE_INDEX_INVALID in
 * zone_ref_out[i].zone_idx.
 *
 * @param[in]  rt           Compiled gemca runtime.
 * @param[in]  x            Particle x-positions, length @p n.
 * @param[in]  y            Particle y-positions, length @p n.
 * @param[in]  z            Particle z-positions, length @p n.
 * @param[in]  ux           Particle x-directions, length @p n.
 * @param[in]  uy           Particle y-directions, length @p n.
 * @param[in]  uz           Particle z-directions, length @p n.
 * @param[in]  n            Number of particles.
 * @param[out] zone_ref_out Receives the zone reference for each particle,
 *                          length @p n.
 */
void osh_gemca_runtime_get_zone_ref_batch(struct osh_gemca_runtime const *rt,
                                          double const *x,
                                          double const *y,
                                          double const *z,
                                          double const *ux,
                                          double const *uy,
                                          double const *uz,
                                          size_t n,
                                          struct osh_zone_ref *zone_ref_out);

/** @} */ /* end defgroup osh_gemca_runtime */

#ifdef __cplusplus
}
#endif

#endif /* OSH_GEMCA_RUNTIME_H */
