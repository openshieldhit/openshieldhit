#ifndef OSH_RAY_H
#define OSH_RAY_H

#include "common/osh_coord.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Ray types and helpers.
 *
 * A ray is the minimal description of a particle in flight: a point in
 * space, a unit direction vector, and a coordinate-system tag.  Energy is
 * carried in p[3] when the ray represents a particle state (as opposed to a
 * pure geometric ray where p[3] is unused).
 *
 * Three representations are provided for different contexts:
 *
 *   struct ray    — legacy form (direction field named cp).  Present only
 *                   for GEMCA compatibility.  Scheduled for removal once the
 *                   GEMCA layer is refactored; use struct ray_v instead.
 *
 *   struct ray_v  — canonical form (direction field named v).  Used
 *                   throughout the beam, transport, and scoring layers.
 *
 *   struct ray_c  — direction expressed as spherical-coordinate cosines
 *                   (cos θ, sin φ, cos φ).  Used internally by the beam
 *                   model where the spherical representation is natural.
 *
 * These structs live here rather than in osh_coord.h because they are
 * particle-transport primitives, not pure coordinate-system machinery.
 * Compare with struct step (osh_step.h), which extends struct ray_v with
 * physics quantities (ds, de, rho) and history context (gen, prim_idx, wt).
 */

/* ---- Struct definitions -------------------------------------------------- */

/**
 * @brief Legacy ray: position + direction via cosine-vector field cp.
 *
 * @details
 * Retained for compatibility with the GEMCA geometry layer, which stores
 * direction as cp (a unit vector, not direction cosines).  New code should
 * use struct ray_v.  Will be removed once GEMCA is refactored.
 */
struct ray {
    double p[3];   /* position [cm] */
    double cp[3];  /* unit direction vector (field name is historical) */
    int system;    /* coordinate system (OSH_COORD_*) */
};

/**
 * @brief Canonical ray: position + energy + unit direction vector v.
 *
 * @details
 * The primary particle-state representation used by transport, beam, and
 * scoring.  p[3] carries total kinetic energy [MeV] when the ray represents
 * a particle in flight; it is unused for geometric rays.
 *
 * The v field is always a unit vector; callers are responsible for
 * normalising before storing.
 */
struct ray_v {
    double p[4]; /* x,y,z [cm] and total kinetic energy [MeV] in p[3] */
    double v[3]; /* unit direction vector */
    unsigned char system; /* coordinate system (OSH_COORD_*) */
};

/**
 * @brief Ray with direction expressed as spherical-coordinate cosines.
 *
 * @details
 * Direction is stored as (cos θ, sin φ, cos φ) following the ISO 80000-2
 * spherical convention used by the beam model.  Use osh_coord_c2v() and
 * osh_coord_v2c() (osh_coord.h) to convert between this form and a unit
 * vector.
 */
struct ray_c {
    double p[4]; /* x,y,z [cm] and total kinetic energy [MeV] in p[3] */
    double c[3]; /* direction cosines: (cos θ, sin φ, cos φ) */
    unsigned char system; /* coordinate system (OSH_COORD_*) */
};

/* ---- Functions ----------------------------------------------------------- */

/**
 * @brief Advance a legacy ray along its direction by distance d.
 *
 * @param[in,out] r  Ray to advance; position updated in place.
 * @param[in]     d  Distance [cm]; may be negative.
 */
void osh_ray_move(struct ray *r, double d);

/**
 * @brief Print a legacy ray to stdout for debugging.
 *
 * @param[in] r  Ray to print.
 */
void osh_ray_print(struct ray const *r);

/**
 * @brief Print a ray_c (spherical-cosine form) to stdout for debugging.
 *
 * @param[in] r  Ray to print (passed by value — struct is small).
 */
void osh_ray_c_print(struct ray_c r);

/**
 * @brief Initialise a ray_c to travel along +Z in PZALIGN coordinates.
 *
 * @details
 * Sets position to the origin and direction cosines to (1, 0, 1), which
 * corresponds to θ = 0 (travel along +Z) with φ = 0.  The coordinate system
 * is set to OSH_COORD_PZALIGN.
 *
 * @param[out] r  Ray to initialise.
 */
void osh_ray_c_clear(struct ray_c *r);

/**
 * @brief Apply a 4×4 affine transform to a ray_v in place.
 *
 * @details
 * Transforms both the position and direction of @p r using the row-major
 * matrix @p t, writing the result to @p rt.  The translation column
 * (t[3], t[7], t[11]) is subtracted from the position — this is the
 * GEMCA/SHIELD-HIT sign convention; callers must store the translation
 * negated if the standard convention is required.
 *
 * @param[in]  r   Input ray in PZALIGN (or any source system).
 * @param[out] rt  Output ray; may not alias r.
 * @param[in]  t   Row-major 4×4 affine matrix.
 *
 * @returns 1 always.
 */
int osh_ray_v_transform(struct ray_v const *r, struct ray_v *rt, double const t[16]);

/**
 * @brief Apply a 4×4 affine transform to a legacy struct ray.
 *
 * @details
 * Identical semantics to osh_ray_v_transform() but for the legacy struct ray
 * (cp direction field).  Temporary — will be removed when struct ray is
 * retired in favour of struct ray_v.
 *
 * @param[in]  r   Input legacy ray.
 * @param[out] rt  Output legacy ray; may not alias r.
 * @param[in]  t   Row-major 4×4 affine matrix.
 *
 * @returns 1 always.
 */
int osh_ray_transform(struct ray const *r, struct ray *rt, double const t[16]);

#ifdef __cplusplus
}
#endif

#endif /* OSH_RAY_H */
