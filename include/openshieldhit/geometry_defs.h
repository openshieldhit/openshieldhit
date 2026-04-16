#ifndef OPENSHIELDHIT_GEOMETRY_DEFS_H
#define OPENSHIELDHIT_GEOMETRY_DEFS_H

/**
 * @file geometry_defs.h
 * @brief Public constants for the cold geometry model.
 *
 * @details
 * Body-type codes and coordinate-system codes used by
 * @ref osh_geometry_body.  The numeric values mirror the internal
 * OSH_GEMCA_BODY_* and OSH_COORD_* constants so that the internal
 * GEMCA layer can consume them without conversion during the migration.
 *
 * When the geometry layer is fully migrated these will be the
 * authoritative definitions; the internal aliases will be removed.
 */

/* ---- Body type codes (osh_geometry_body.type) ---------------------------- */

#define OSH_GEOMETRY_BODY_NONE 0  /**< Not set / unknown body type. */
#define OSH_GEOMETRY_BODY_SPH  1  /**< Sphere. */
#define OSH_GEOMETRY_BODY_WED  2  /**< Wedge. */
#define OSH_GEOMETRY_BODY_ARB  3  /**< Arbitrary convex polyhedron. */
#define OSH_GEOMETRY_BODY_BOX  4  /**< Axis-aligned box. */
#define OSH_GEOMETRY_BODY_VOX  5  /**< Voxel body (reserved — not yet usable via public API). */
#define OSH_GEOMETRY_BODY_RPP  6  /**< Rectangular parallelepiped. */
#define OSH_GEOMETRY_BODY_RCC  7  /**< Right circular cylinder. */
#define OSH_GEOMETRY_BODY_REC  8  /**< Right elliptical cylinder. */
#define OSH_GEOMETRY_BODY_TRC  9  /**< Truncated right cone. */
#define OSH_GEOMETRY_BODY_ELL  10 /**< Ellipsoid. */
#define OSH_GEOMETRY_BODY_YZP  11 /**< Infinite plane perpendicular to X (YZ plane). */
#define OSH_GEOMETRY_BODY_XZP  12 /**< Infinite plane perpendicular to Y (XZ plane). */
#define OSH_GEOMETRY_BODY_XYP  13 /**< Infinite plane perpendicular to Z (XY plane). */
#define OSH_GEOMETRY_BODY_PLA  14 /**< Arbitrary oriented infinite plane. */
#define OSH_GEOMETRY_BODY_ROT  15 /**< Transform helper: rotation of a referenced body. */
#define OSH_GEOMETRY_BODY_CPY  16 /**< Transform helper: copy/reference of an existing body. */
#define OSH_GEOMETRY_BODY_MOV  17 /**< Transform helper: translation of a referenced body. */

/* ---- Coordinate-system codes (osh_geometry_body.coord) ------------------- */

/**
 * @defgroup osh_geometry_coord Coordinate system codes
 * @{
 *
 * These values identify the coordinate system in which a body's raw
 * argument array (@ref osh_geometry_body.a) is expressed.
 */

#define OSH_GEOMETRY_COORD_UNKNOWN  0 /**< Not set or not applicable. */
#define OSH_GEOMETRY_COORD_UNIVERSE 1 /**< Simulation universe coordinates (default). */
#define OSH_GEOMETRY_COORD_PZALIGN  2 /**< Particle Z-aligned system (used by some straggling models). */
#define OSH_GEOMETRY_COORD_VOXELCT  3 /**< CT coordinate system (voxel bodies). */
#define OSH_GEOMETRY_COORD_BZALIGN  4 /**< Body-corner-at-origin, aligned along Z. */
#define OSH_GEOMETRY_COORD_BCALIGN  5 /**< Body-centre-at-origin. */

/** @} */

#endif /* OPENSHIELDHIT_GEOMETRY_DEFS_H */
