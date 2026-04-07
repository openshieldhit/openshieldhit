#ifndef _OSH_GEMCA2
#define _OSH_GEMCA2

#include <float.h>
#include <stdio.h>

#include "common/osh_rc.h"
#include "transport/osh_transport.h" /* needed for struct ray */

/** @defgroup gemca Geometry Engine (GEMCA)
 *
 * @brief Constructive solid geometry (CSG) engine for zone and body lookup.
 *
 * @details Geometry is described as a set of bodies (primitive solids) and zones (CSG
 *          combinations of bodies). Bodies are combined into zones via boolean operators
 *          (union `|`, intersection `+`, difference `-`) stored as an abstract syntax
 *          tree (AST) of cgnodes. User-facing references to bodies, zones, and materials
 *          are always strings; dense 0-based array indices are internal and appear only
 *          in debug output or resolver code.
 *
 * @{
 */

/** @name Internal node-type tags (do not use outside gemca/) */
/** @{ */
#define _OSH_GEMCA_CGNODE_BODY      0 /**< cgnode is a leaf holding a body pointer */
#define _OSH_GEMCA_CGNODE_COMPOSITE 1 /**< cgnode is an interior node with left/right children and an operator */
/** @} */

/** @name Numeric constants */
/** @{ */
#ifdef INFINITY
#define OSH_GEMCA_INFINITY INFINITY /**< Positive infinity, used as a sentinel "no intersection" distance */
#else
#define OSH_GEMCA_INFINITY 1e300
#endif

#define OSH_GEMCA_SMALL    1e-12 /**< Epsilon for floating-point comparisons */
#define OSH_GEMCA_STEPLIM  1e-8  /**< Minimum step size to avoid getting stuck on a surface due to numerical precision */

#define OSH_GEMCA_ZONE_INDEX_INVALID ((size_t) -1) /**< Sentinel returned when no zone contains the ray */
/** @} */

/**
 * @brief Top-level workspace holding all parsed geometry.
 *
 * @details Allocated and populated by osh_gemca_workspace_init() and osh_gemca_load().
 *          Released by osh_gemca_workspace_free(). Zones and bodies are owned by this
 *          struct and must not be freed independently.
 */
struct gemca_workspace {
    struct body **bodies; /**< Array of pointers to all body primitives */
    struct zone **zones;  /**< Array of pointers to all zones */
    size_t nbodies;       /**< Number of entries in bodies[] */
    size_t nzones;        /**< Number of entries in zones[] */
    char *filename;       /**< Path to the geo.dat file that was loaded */
};

/**
 * @brief A body primitive: a solid defined by one or more surfaces.
 *
 * @details A body is the leaf element of the CSG tree. It is described by a set of
 *          half-space surfaces; a point is inside the body if it satisfies all surface
 *          conditions simultaneously. The transformation matrix `t` maps from
 *          OSH_COORD_UNIVERSE to the body's local coordinate system.
 */
struct body {
    double t[16];           /**< 4x4 row-major transformation matrix (universe → body-local coords) */
    struct surface **surfs; /**< Array of pointers to the surfaces that bound this body */
    size_t lineno;          /**< Line number in geo.dat where this body was defined */
    char *name;             /**< User-given name (string key, never converted to an integer ID) */
    char *filename_vox;     /**< Path to voxel file, only set for voxel body types */
    double *a;              /**< Raw argument list as given in geo.dat */
    int na;                 /**< Number of entries in a[] */
    int type;               /**< Body type identifier (OSH_GEMCA_BODY_* defines) */
    int nsurfs;             /**< Number of entries in surfs[] */
    char coord;             /**< Coordinate system of body parameters (OSH_COORD_* value) */
};

/**
 * @brief A node in the CSG abstract syntax tree.
 *
 * @details Leaf nodes (type == _OSH_GEMCA_CGNODE_BODY) hold a pointer to a body.
 *          Interior nodes (type == _OSH_GEMCA_CGNODE_COMPOSITE) combine two child nodes
 *          with a boolean operator stored in `op`. The tree is evaluated recursively to
 *          determine zone membership and distance to the nearest boundary.
 */
struct cgnode {
    double bb_max[3];     /**< Bounding box maximum corner (TODO: not yet used) */
    double bb_min[3];     /**< Bounding box minimum corner (TODO: not yet used) */
    struct cgnode *left;  /**< Left child; non-NULL only for composite nodes */
    struct cgnode *right; /**< Right child; non-NULL only for composite nodes */
    struct body *body;    /**< Body pointer; non-NULL only for leaf nodes */
    int type;             /**< Node type: _OSH_GEMCA_CGNODE_BODY or _OSH_GEMCA_CGNODE_COMPOSITE */
    char op;              /**< CSG operator for composite nodes: `|` union, `+` intersection, `-` difference */
};

/**
 * @brief A zone: a named CSG region with an associated material.
 *
 * @details A zone is defined by a boolean expression over bodies, stored as a cgnode
 *          AST rooted at `node`. User-facing names for zones and materials are strings;
 *          `material_idx` is resolved after both geo.dat and mat.dat have been parsed.
 *
 * @note Per-zone transport cutoff overrides (tcut, ncut) are planned but not yet
 *       implemented. When added they will take precedence over the global defaults in
 *       beam_workspace at simulation time.
 */
struct zone {
    struct cgnode node;   /**< Root of the CSG AST describing this zone's geometry */
    size_t material_idx;  /**< Dense material array index; (size_t)-1 until post-parse resolution */
    size_t lineno;        /**< Line number in geo.dat where this zone was first defined */
    size_t ntokens;       /**< Number of tokens in the raw zone expression */
    char **tokens;        /**< Tokenised zone expression (reversed, for stack-based parsing) */
    char *name;           /**< User-given zone name (string key) */
    char *material_name;  /**< User-facing material name assigned to this zone */
};

/**
 * @brief A half-space surface primitive.
 *
 * @details Surfaces are the leaf elements of a body. A body is defined as the
 *          intersection of all its half-spaces. The meaning of p[] depends on `type`.
 */
struct surface {
    double *p;    /**< Surface parameters; interpretation depends on type */
    double _dist; /**< Cached distance to surface along current ray; negative if no crossing */
    int np;       /**< Number of entries in p[] */
    int type;     /**< Surface type identifier (OSH_GEMCA_SURF_* defines) */
};

/**
 * @brief Allocate and zero-initialise a gemca workspace.
 *
 * @param[out] wg  Receives a pointer to the newly allocated workspace.
 *
 * @returns OSH_OK on success, OSH_ENOMEM if allocation fails.
 */
enum osh_status osh_gemca_workspace_init(struct gemca_workspace **wg);

/**
 * @brief Free a gemca workspace and all geometry it owns.
 *
 * @param[in] wg  Workspace to free. Passing NULL is a no-op.
 *
 * @returns OSH_OK.
 */
enum osh_status osh_gemca_workspace_free(struct gemca_workspace *wg);

/**
 * @brief Parse a geo.dat file and populate a gemca workspace.
 *
 * @param[in]  filename  Path to the geometry file.
 * @param[out] g         Workspace to populate (must have been initialised with osh_gemca_workspace_init()).
 *
 * @returns OSH_OK on success, OSH_EPARSE on parse error, OSH_ENOMEM on allocation failure.
 */
enum osh_status osh_gemca_load(char const *filename, struct gemca_workspace *g);

/**
 * @brief Return the dense zone index for the zone containing a ray.
 *
 * @param[in] g  Gemca workspace.
 * @param[in] r  Ray whose position is tested.
 *
 * @returns 0-based index into g->zones[], or OSH_GEMCA_ZONE_INDEX_INVALID if the
 *          ray is outside all defined zones. Callers must check for
 *          OSH_GEMCA_ZONE_INDEX_INVALID before indexing g->zones[].
 */
size_t osh_gemca_get_zone_index(struct gemca_workspace *g, struct ray *r);

/**
 * @brief Return the distance a ray travels inside zone `z` before leaving it.
 *
 * @param[in] z  Zone the ray is currently inside.
 * @param[in] r  Ray (position and direction).
 *
 * @returns Total path length inside the zone, in the same units as the geometry.
 */
double osh_gemca_get_distance(struct zone *z, struct ray const *r);

/** @name Debug printing */
/** @{ */
void osh_gemca_print_gemca(struct gemca_workspace const *g); /**< Print full workspace summary to debug log */
void osh_gemca_print_body(struct body const *b);             /**< Print body parameters to debug log */
void osh_gemca_print_zone(struct zone const *z);             /**< Print zone parameters to debug log */
void osh_gemca_print_surface(struct surface const *s);       /**< Print surface parameters to debug log */
void osh_gemca_print_cgnodes(struct cgnode const *self);     /**< Recursively print CSG tree to debug log */
/** @} */

/** @} */ /* end defgroup gemca */

#endif /* gemca2.h */
