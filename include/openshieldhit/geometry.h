#ifndef OPENSHIELDHIT_GEOMETRY_H
#define OPENSHIELDHIT_GEOMETRY_H

#include <stddef.h>
#include <stdint.h>

#include "openshieldhit/diag.h"
#include "openshieldhit/geometry_defs.h"
#include "openshieldhit/status.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Sentinel zone index returned when a ray is outside all zones.
 */
#define OSH_ZONE_INDEX_INVALID ((size_t) -1)

/**
 * @brief Geometry location of a particle: zone identity, material assignment,
 *        and optional raw HU value for voxel zones.
 *
 * @details
 * Filled by @ref osh_gemca_runtime_get_zone_ref_batch for each particle.
 * For analytic zones @p has_hu is 0 and @p material_idx is the zone's
 * ASSIGNMAT index.  For voxel zones @p has_hu is 1, @p hu holds the raw
 * (clamped to [-1000, 1600]) Hounsfield value of the current voxel, and
 * @p material_idx is resolved via the HU→bin LUT.
 *
 * @p zone_idx equals @ref OSH_ZONE_INDEX_INVALID when the particle is
 * outside all defined zones.
 */
struct osh_zone_ref {
    size_t zone_idx;     /**< Zone index; OSH_ZONE_INDEX_INVALID if outside geometry. */
    size_t material_idx; /**< Resolved material index: ASSIGNMAT for analytic, HU-bin for voxel. */
    int16_t hu;          /**< Clamped HU value [−1000, 1600]; meaningful only when has_hu != 0. */
    char has_hu;         /**< Non-zero for voxel zones; zero for analytic zones. */
};

/**
 * @file geometry.h
 * @brief Public cold-data geometry model for `osh_core`.
 *
 * @details
 * This header defines the geometry input structs that library users and
 * frontend adapters fill before calling @ref osh_geometry_workspace_prepare().
 *
 * Design rules (mirror those of beam.h):
 * - Types here are the public cold model: user-visible names, body types,
 *   raw argument arrays, and raw zone boolean expressions.
 * - Derived state (surface lists, transformation matrices, CSG AST nodes,
 *   compiled RPN bytecode) lives in core-private prepared/runtime structs and
 *   must not be written by callers.
 * - The `prepared` pointer in @ref osh_geometry_workspace is reserved for
 *   core-internal use; callers must treat it as opaque.
 * - Parse-time scratch (line numbers, token arrays, file names) is not part
 *   of the public model.
 * - Voxel bodies (@ref OSH_GEOMETRY_BODY_VOX) are recognised by the type
 *   code but their full support is not yet exposed through this API.
 */

/** @cond PRIVATE */
struct osh_gemca_prepared; /* defined in src/gemca/osh_gemca2.h */

/** @endcond */

/**
 * @brief One cold body primitive description.
 *
 * @details
 * Carries only the user-facing description needed to define the body.
 * Derived state (surfaces, transformation matrix, bounding box) is computed
 * by @ref osh_geometry_workspace_prepare() and stored privately.
 *
 * The @p a array and @p name string are owned by this struct and freed by
 * @ref osh_geometry_workspace_free().
 *
 * For DCM-backed voxel bodies, @p hu points to a heap allocation of HU values
 * [int16] in the compiled voxel layout. For row-major this is `nx*ny*nz`;
 * for tiled layouts such as Morton-8 it may include zero-padded boundary
 * slots, so @p n_hu is the storage length rather than the physical voxel
 * count. The workspace owns and frees this memory. For non-voxel bodies,
 * @p hu is NULL and @p n_hu is 0.
 */
struct osh_geometry_body {
    char *name;  /**< User-given body name (null-terminated, unique within the workspace). */
    double *a;   /**< Raw argument list as given in the geometry input; length is @p na. */
    int16_t *hu; /**< Optional HU voxel array owned by the workspace (DCM/VOX bodies). */
    size_t n_hu; /**< Number of entries in @p hu (0 when hu is NULL). */
    int na;      /**< Number of entries in @p a[]. */
    int type;    /**< Body type: one of the @ref OSH_GEOMETRY_BODY_* codes. */
    char coord;  /**< Coordinate system of body parameters: one of the @ref OSH_GEOMETRY_COORD_* codes. */
};

/**
 * @brief One cold zone description.
 *
 * @details
 * A zone is a named region of space defined by a boolean combination of
 * bodies.  The boolean expression is stored as the raw input string so that
 * callers (parsers, JSON adapters, …) do not need to know the internal
 * tokenisation format.  @ref osh_geometry_workspace_prepare() parses and
 * compiles @p expr into the internal CSG representation.
 *
 * All string fields are owned by this struct and freed by
 * @ref osh_geometry_workspace_free().
 */
struct osh_geometry_zone {
    char *name;          /**< User-given zone name (null-terminated). */
    char *material_name; /**< Name of the material assigned to this zone (null-terminated). */
    char *expr;          /**< Raw boolean body expression (e.g. @c "+BODY1 -BODY2 | +BODY3"). */
};

/**
 * @brief Public cold geometry workspace.
 *
 * @details
 * The input object populated by frontend adapters or application parsers
 * before calling @ref osh_geometry_workspace_prepare().
 *
 * The @p bodies and @p zones arrays are flat value arrays (not pointer
 * arrays) owned by this struct.  Resizing them requires reallocating via
 * the internal helper @c osh_geometry_workspace_add_body() /
 * @c osh_geometry_workspace_add_zone() (not yet exposed).
 *
 * The @p prepared pointer is reserved for core-internal prepared state and
 * must be treated as opaque by callers.  It is NULL until
 * @ref osh_geometry_workspace_prepare() succeeds.
 */
struct osh_geometry_workspace {
    struct osh_geometry_body *bodies;    /**< Flat array of body descriptions; length is @p nbodies. */
    struct osh_geometry_zone *zones;     /**< Flat array of zone descriptions; length is @p nzones. */
    size_t nbodies;                      /**< Number of entries in @p bodies[]. */
    size_t nzones;                       /**< Number of entries in @p zones[]. */
    struct osh_gemca_prepared *prepared; /**< Internal prepared state; owned by core. */
};

/**
 * @brief Allocate a geometry workspace with defaults.
 *
 * @param[out] ws_out  Receives a pointer to the newly allocated workspace.
 *                     The workspace is zero-initialised; @p prepared is NULL.
 *
 * @returns OSH_OK on success, or an OSH_E* code on failure.
 */
enum osh_status osh_geometry_workspace_create(struct osh_geometry_workspace **ws_out);

/**
 * @brief Build internal prepared state from a cold geometry workspace.
 *
 * @details
 * Validates body types and argument counts, computes transformation
 * matrices and surface lists, tokenises and compiles zone boolean
 * expressions into the internal CSG representation, and resolves
 * geometry-material linkage.
 *
 * Must be called before passing the workspace to the geometry runtime.
 * Calling it a second time on an already-prepared workspace is an error
 * (returns an OSH_E* code without modifying state).
 *
 * @param[in,out] ws    Workspace to prepare.
 * @param[in]     diag  Borrowed diagnostics sink for preparation messages, or NULL.
 *
 * @returns OSH_OK on success, or an OSH_E* code on failure.
 */
enum osh_status osh_geometry_workspace_prepare(struct osh_geometry_workspace *ws, struct osh_diag_sink const *diag);

/**
 * @brief Free a geometry workspace and all geometry it owns.
 *
 * @details
 * Frees all body arguments and names, all zone strings, internal prepared
 * state, and the workspace struct itself.  Passing NULL is a no-op.
 *
 * @param[in] ws  Workspace to free.
 */
enum osh_status osh_geometry_workspace_free(struct osh_geometry_workspace *ws);

#ifdef __cplusplus
}
#endif

#endif /* OPENSHIELDHIT_GEOMETRY_H */
