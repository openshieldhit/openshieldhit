#include "gemca/osh_gemca2_dist.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "common/osh_coord.h"
#include "common/osh_exit.h"
#include "common/osh_logger.h"
#include "common/osh_vect.h"
#include "gemca/osh_gemca2.h"
#include "gemca/osh_gemca2_calc_surface.h"
#include "gemca/osh_gemca2_defines.h"

#define MIN(a, b) (((a) < (b)) ? (a) : (b))

static inline double _dist_zone(struct cgnode *self, struct ray const *r);
static inline double _dist_body(struct body const *b, struct ray const *r);
static inline double _dist_surface(struct surface const *sf, struct ray const *r);

static inline int _transform_to_local(struct body const *b, struct ray const *r, struct ray *tr);
static inline int _ray_advance(double d, struct ray const *r, struct ray *rr); // TODO: move to osh_transport.h
static inline double _inside_body(struct body const *b, struct ray const *r);

static inline double _dist_plane_xyz(int axis, struct surface const *sf, struct ray const *r);
static inline double _dist_plane(struct surface const *sf, struct ray const *r);
static inline double _dist_sphere(double r2, struct ray const *r);
static inline double _dist_cyl(double r2, struct ray const *r);
static inline double _dist_elipcyl(double ra2, double rb2, struct ray const *r);
static inline double _dist_cone(double ra2, double rb2, struct ray const *r);
static inline double _dist_ellipsoid(double ra2, double rb2, double rc2, struct ray const *r);

static inline double _quadratic_solver(double a, double b, double c);
static inline double _minpos(double a, double b);

/**
 * @brief For a given ray and given zone, get smallest postive distance to zone surface along ray.
 *
 * @details The ray is forward propagated until it leaves the zone given by `z`.
 *          The total distance travelled inside the zone is returned.
 *
 * @param[in] z - current zone the ray is in
 * @param[in] r - a ray
 *
 * @returns distance to next zone boundary
 *
 * @author Niels Bassler
 */
double osh_gemca_get_distance(struct zone *z, struct ray const *r) {

    double d, total_distance;
    struct ray rr;

    // printf("osh_gemca_get_distance(): calculating distance to zone boundary for zone '%s'\n", z->name);

    total_distance = 0.0;
    rr = *r; /* make a copy of the ray */

    osh_vect_norm(rr.cp); /* normalize the direction vector */

    /* per definition, we are always inside this zone *z we begin with */
    z->node._is_inside = 1; // is this really needed?

    while (1) {
        d = _dist_zone(&z->node, &rr); /* find shortest distance to closest body */
        if (!z->node._is_inside) {     /* keep advancing until we left the zone */
            break;
        }
        // printf("  Currently inside zone '%s', advancing %.9e to next boundary\n", z->name, d);
        if (d < 0.0) {
            osh_error(EX_SOFTWARE, "osh_gemca_get_distance(): negative distance to zone boundary");
        }
        if (d < OSH_GEMCA_STEPLIM) {
            // TODO: check if this may introduce scoring artefacts when a step is half in two zones
            d = OSH_GEMCA_STEPLIM; /* avoid getting stuck on surface due to numerical precision */
        }
        total_distance += d;
        /* advance ray once */
        _ray_advance(d, &rr, &rr);
    }
    return total_distance;
}

static inline double _dist_zone(struct cgnode *self, struct ray const *r) {
    double d;
    double d1;
    double d2;
    struct ray tr; /* transformed ray to coordinate system of the body */

    if (self->type == _OSH_GEMCA_CGNODE_BODY) {  /* we are in the leaf node */
        _transform_to_local(self->body, r, &tr); /* Transform ray to body's coordinate system */

        // osh_gemca_print_body(self->body);
        // osh_transport_print_ray(r);
        // osh_transport_print_ray(&tr);

        d = _dist_body(self->body,
                       &tr); /* Calculate distance to body's surface. It may, or may not be at a zone boudnary. */
        self->_is_inside = _inside_body(self->body, &tr); /* Check if ray is inside the body */
        // printf("  is inside body '%s' = %d\n", self->body->name, self->_is_inside);
        // printf("  Shortest distance to body '%s' is %f\n", self->body->name, d);
        return d;
    } else {
        d1 = _dist_zone(self->left, r);  // Distance to left child
        d2 = _dist_zone(self->right, r); // Distance to right child

        /* Apply the appropriate operation based on self->op */
        /* This is following table 3 in Scott D. Roth's algorithm from "Ray Casting for Modeling Solids"
           (Computer Graphics, Vol. 18, No. 3, July 1982) */
        switch (self->op) {
        case '|': // Union
            self->_is_inside = self->left->_is_inside || self->right->_is_inside;
            break;
        case '+': // Intersection
            self->_is_inside = self->left->_is_inside && self->right->_is_inside;
            break;
        case '-': // Difference
            self->_is_inside = self->left->_is_inside && !self->right->_is_inside;
            break;
        default:
            osh_error(EX_SOFTWARE, "_dist_zone(): unknown operator");
            break;
        }
        /* return smallest possible distance */
        // printf("Result of operation %c: _is_inside = %d\n", self->op, self->_is_inside);

        return _minpos(d1, d2);
    }
    return OSH_GEMCA_INFINITY; /* Indicate no intersection found or operation not supported */
}

/* check if ray is inside a body */
static inline double _inside_body(struct body const *b, struct ray const *r) {
    int i;
    struct surface *sf;
    int inside;

    // printf("_inside_body '%s'?\n", b->name);
    // osh_transport_print_ray(r);

    for (i = 0; i < b->nsurfs; i++) {

        sf = b->surfs[i]; /* get surface number i*/
        // printf("  Checking surface %d:   type: %d\n", i, sf->type);
        inside = osh_gemca2_check_surface_side(sf, r); /* Check if ray is inside the surface */
        // printf("  .... RESULT Inside surface %d: %d\n", i, inside);
        if (!inside) {
            return 0;
        }
    }
    return 1;
}

/* calculate closest distance to this body, OSH_GEMCA_INFINITY if no intersection */
static inline double _dist_body(struct body const *b, struct ray const *r) {

    struct surface *sf;
    int i;
    double d, _d;
    struct ray rr;

    d = OSH_GEMCA_INFINITY;
    _d = OSH_GEMCA_INFINITY;

    // printf("******* _dist_body() %s\n", b->name);
    // osh_gemca_print_body(b);

    for (i = 0; i < b->nsurfs; i++) {

        sf = b->surfs[i]; /* get surface number i*/

        d = _dist_surface(sf, r); /* calculate distance to surface along ray (postive direction only) */
        // printf("  Distance to surface %i is %f\n", i, _d);

        if (d > 0.0 && d < _d) { /* check if positive and closer than current minimal distance */
            _d = d;

            // DEBUG: show the next intersection point
            _ray_advance(_d, r, &rr); /* rr will have next intersection point */
            // printf("  Intersection at %f %f %f in body %s\n", rr.p[0], rr.p[1], rr.p[2], b->name);

        } /* end if distance is positve */
    } /* end loop over each surface */
    return _d;
}

/**
 * @brief Calculate the distance to a surface along the ray.
 *
 * @param[in] sf - surface to be checked
 * @param[in] r - ray to be checked. Ray must be in transformed body->system coords, and normalized.
 *
 * @return distance to surface along the ray
 */
static inline double _dist_surface(struct surface const *sf, struct ray const *r) {
    double d = OSH_GEMCA_INFINITY;

    switch (sf->type) {
    case OSH_GEMCA_SURF_SPHERE:
        d = _dist_sphere(sf->p[0], r);
        break;
    case OSH_GEMCA_SURF_ELLIPSOID:
        d = _dist_ellipsoid(sf->p[0], sf->p[1], sf->p[2], r);
        break;
    case OSH_GEMCA_SURF_CYLZ:
        d = _dist_cyl(sf->p[0], r);
        break;
    case OSH_GEMCA_SURF_ELLZ:
        d = _dist_elipcyl(sf->p[0], sf->p[1], r);
        break;
    case OSH_GEMCA_SURF_CONE:
        d = _dist_cone(sf->p[0], sf->p[1], r);
        break;
    case OSH_GEMCA_SURF_PLANEX:
        d = _dist_plane_xyz(0, sf, r);
        break;
    case OSH_GEMCA_SURF_PLANEY:
        d = _dist_plane_xyz(1, sf, r);
        break;
    case OSH_GEMCA_SURF_PLANEZ:
        d = _dist_plane_xyz(2, sf, r);
        // printf(" Surface parameters (transformed): %f %f\n", sf->p[0], sf->p[1]);
        // printf(" Ray coordinates (transformed): %f %f %f - %f %f %f\n", r->p[0], r->p[1], r->p[2], r->cp[0],
        // r->cp[1],
        //        r->cp[2]);
        // printf(" Distance to plane: %f\n", d);
        break;
    case OSH_GEMCA_SURF_PLANE: // TODO: refactor to OSH_GEMCA_SURF_PLANE
        _dist_plane(sf, r);
        break;
    default:
        break;
    }

    return d;
}

/**
 * @brief Transform ray according to surface type and its coordinates.
 *
 * @details
 *
 * @param[in] b - body parameters incl. its transformation matrix
 * @param[in] r - input ray in OSH_COORD_UNIVERSE
 * @param[out] tr - transformed output ray in system given by b->coord
 *
 * @returns 1
 *
 * @author Niels Bassler
 */
static inline int _transform_to_local(struct body const *b, struct ray const *r, struct ray *tr) {

    int i;
    int j;

    // printf("_transform_to_local() %s   coord: %i\n", b->name, b->coord);
    //  For now, just copy all elements of the ray. Later this can be optimized.
    for (i = 0; i < 3; i++) {
        tr->p[i] = r->p[i];
        tr->cp[i] = r->cp[i];
    }
    tr->system = b->coord;
    // printf("_transform_to_local() ray after transform:\n");
    // osh_transport_print_ray(tr);

    // then overwrite the values which may change:
    switch (b->coord) {
    case OSH_COORD_UNIVERSE:
        break;

    case OSH_COORD_BCALIGN:
        /* simple translation */
        for (i = 0; i < 3; i++) {
            j = i * 4;
            tr->p[i] = r->p[i] + b->t[j + 3]; // notice, that in osh_coord.h see comment
            tr->cp[i] = r->cp[i];
            // printf(" tr->p[%i] = %f + %f = %f\n", i, r->p[i],  b->t[j+3], r->p[i] + b->t[j+3]);
        }
        break;

    case OSH_COORD_BZALIGN:
        /* simple translation and rotation, so we have to use osh_coord_trans_ray */
        osh_coord_trans_ray_r(r, tr, b->t);
        break;

    default:
        osh_error(EX_SOFTWARE, "_transform_to_local() unsupported coordinate system :%i", b->coord);
        break;
    }
    return 1;
}

/**
 * @brief Calculates a new ray at a distance `d` along its path.
 *
 * @param[in] d - The distance along the ray to the intersection point.
 * @param[in] r - The original ray.
 * @param[in,out] rr - The new ray at the intersection point, with the same direction as `r`. Must have been allocated
 * before calling this function.
 *
 * @return 0 if successful, -1 if an error occurred.
 *
 * @author Niels Bassler
 *
 */
static inline int _ray_advance(double d, struct ray const *r, struct ray *rr) {
    int i;

    if (d < 0) {
        return -1; /* Invalid distance */
    }

    for (i = 0; i < 3; i++) {
        rr->p[i] = r->p[i] + r->cp[i] * d; /* Calculate new position */
        rr->cp[i] = r->cp[i];              /* Copy the direction */
    }
    rr->system = r->system; /* Preserve the coordinate system */

    return 0;
}

/**
 * @brief For a given ray and a given plane, calculate the signed distance to the plane.
 *
 * @param[in] axis - axis of the plane (0=x, 1=y, 2=z)
 * @param[in] sf - surface parameters
 * @param[in] r - ray which may or may not intersect. ray->cp must be normalized and in body coordinate system.
 *
 * @returns signed distance to intersection, 0.0 if they are parallel in plane, and OSH_GEMCA_INFITIE if they are
 * parallel but not in plane.
 *
 * @author Niels Bassler
 */
static inline double _dist_plane_xyz(int axis, struct surface const *sf, struct ray const *r) {
    double d, _d;

    // Choose the appropriate ray and plane components
    double rp = r->p[axis];
    double rcp = r->cp[axis];

    _d = (sf->p[0] * rcp);
    if (fabs(_d) < OSH_GEMCA_SMALL) { /* ray is parallel to plane */
        _d = sf->p[0] * rp + sf->p[1];
        if (fabs(_d) < OSH_GEMCA_SMALL) { /* check if ray is on the plane */
            return 0.0;
        }
        return OSH_GEMCA_INFINITY;
    }
    d = -(sf->p[0] * rp + sf->p[1]) / _d;

    return d;
}

/**
 * @brief For a given ray and a given plane, calculate the signed distance to the plane along the ray.
 *
 * @details d = - (p dot n + D) / (l dot n) where p is the ray origin, n is the plane normal, D is the plane offset, and
 * l is the ray direction.
 * @param[in] sf - surface holding the plane parameters Ax + By + Cz + D = 0
 * @param[in] r - ray which may or may not intersect. ray->cp must be normalized. Must be in body coordinate system.
 *
 * @returns signed distance to intersection, 0.0 if they are parallel and on the plane,
 *          OSH_GEMCA_INFINITY if they are parallel but not on the plane.
 *
 * @author Niels Bassler
 */
static inline double _dist_plane(struct surface const *sf, struct ray const *r) {
    double dot_ln, dot_pn, d;
    double *n = &(sf->p[0]); /* Plane normal vector (A, B, C) */

    /* Dot product of ray direction with plane normal */
    dot_ln = osh_vect_dot(r->cp, n);
    if (fabs(dot_ln) < OSH_GEMCA_SMALL) {          /* Check if ray is parallel to the plane */
        dot_pn = osh_vect_dot(r->p, n) + sf->p[3]; /* Dot product of ray origin with plane normal */
        if (fabs(dot_pn) < OSH_GEMCA_SMALL) {      /* Ray lies on the plane */
            return 0.0;
        }
        return OSH_GEMCA_INFINITY; /* Ray is parallel but not on the plane */
    }

    /* Compute distance to the plane along the ray*/
    dot_pn = osh_vect_dot(r->p, n) + sf->p[3];
    d = -dot_pn / dot_ln;
    return d;
}

/**
 * @brief For a given ray and a sphere centered at (0,0,0), calculate the closest signed distance.
 *
 * @details The case of a single intersection (tangential) is treated as no intersection.
 *          The ray must be normalized. Returns OSH_GEMCA_INFINITY if there is no valid intersection.
 *
 * @param[in] r2 - squared radius of the sphere
 * @param[in] r - ray which may or may not intersect. ray->cp must be normalized.
 *
 * @returns signed distance to intersection, or OSH_GEMCA_INFINITY if no valid intersection.
 *
 * @author Niels Bassler
 */
static inline double _dist_sphere(double r2, struct ray const *r) {
    double b, c;

    // Compute quadratic coefficients
    b = 2.0 * osh_vect_dot(r->cp, r->p); // 2(l · o)
    c = osh_vect_len2(r->p) - r2;        // ||o||^2 - r^2

    // Use the quadratic solver to find the smallest positive root
    return _quadratic_solver(1.0, b, c);
}

/**
 * @brief For a given ray and a given plane, calculate the closest signed distance to a cylinder
 *
 * @details Cylinder is positioned along Z-axis.
 *
 * @param[in] r2 - radius^2 of cylinder
 * @param[in] r - ray which may or may not intersect. ray->cp must be normalized.
 *
 * @see https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
 *
 * @returns signed distance to intersection, 0.0 if no or only touching the suface (1 intersection).
 *
 * @author Niels Bassler
 */
static inline double _dist_cyl(double r2, struct ray const *r) {

    double a, b, c;

    a = r->cp[0] * r->cp[0] + r->cp[1] * r->cp[1];
    b = 2.0 * (r->cp[0] * r->p[0] + r->cp[1] * r->p[1]);
    c = r->p[0] * r->p[0] + r->p[1] * r->p[1] - r2;

    return _quadratic_solver(a, b, c);
}

/**
 * @brief For a given ray and a given plane, calculate the closest signed distance to an elliptical cylinder
 *
 * @details
 *
 * @param[in] ra2 - radius1^2 of cylinder
 * @param[in] rb2 - radius2^2 of cylinder
 * @param[in] r - ray which may or may not intersect. ray->cp must be normalized.
 *
 * @returns signed distance to intersection, 0.0 if no or only touching the suface (1 intersection).
 *
 * @author Niels Bassler
 */
static inline double _dist_elipcyl(double ra2, double rb2, struct ray const *r) {

    double a, b, c;

    // TODO vectorize me
    a = (r->cp[0] * r->cp[0]) / ra2 + (r->cp[1] * r->cp[1]) / rb2; /* a in gemca */
    b = (r->cp[0] * r->p[0]) / ra2 + (r->cp[1] * r->p[1]) / rb2;   /* b in gemca */
    c = (r->p[0] * r->p[0]) / ra2 + (r->p[1] * r->p[1]) / rb2 - 1; /* c in gemca */

    return _quadratic_solver(a, b, c);
}

/**
 * @brief For a given ray and a cone, calculate the closest positve distance to an elliptical cylinder
 *
 * @details
 *
 * @param[in] ra2 - radius1^2 of cone
 * @param[in] rb2 - radius2^2 of cone
 * @param[in] r - ray which may or may not intersect. ray->cp must be normalized.
 *
 * @returns signed distance to intersection, 0.0 if no or only touching the suface (1 intersection).
 *
 * @author Niels Bassler
 */
static inline double _dist_cone(double ra2, double rb2, struct ray const *r) {

    double a, b, c;
    double t;

    t = (r->p[2] - ra2) / rb2;

    // TODO vectorize me
    a = (r->cp[0] * r->cp[0]) + (r->cp[1] * r->cp[1]) + (r->cp[2] * r->cp[2]) / rb2; /* a in gemca */
    b = (r->cp[0] * r->p[0]) + (r->cp[1] * r->p[1]) - t * r->cp[2];                  /* b in gemca */
    c = (r->p[0] * r->p[0]) + (r->p[1] * r->p[1]) - t * t * rb2;                     /* c in gemca */

    return _quadratic_solver(a, b, c);
}

/**
 * @brief For a given ray and a given plane, calculate the closest signed distance to a cone
 *
 * @details
 *
 * @param[in] radius - radius of cylinder
 * @param[in] r - ray which may or may not intersect. ray->cp must be normalized.
 *
 * @author Niels Bassler
 */
static inline double _dist_ellipsoid(double ra2, double rb2, double rc2, struct ray const *r) {

    double a, b, c;

    // TODO vectorize me
    a = (r->cp[0] * r->cp[0]) / ra2 + (r->cp[1] * r->cp[1]) / rb2 + (r->cp[2] * r->cp[2]) / rc2; /* a in gemca */
    b = (r->cp[0] * r->p[0]) / ra2 + (r->cp[1] * r->p[1]) / rb2 + (r->cp[2] * r->p[2]) / rc2;    /* b in gemca */
    c = (r->p[0] * r->p[0]) / ra2 + (r->p[1] * r->p[1]) / rb2 + (r->p[2] * r->p[2]) / rc2 - 1;   /* c in gemca */

    return _quadratic_solver(a, b, c);
}

/**
 * @brief return smallest positive solution to a*x^2 + b*x + c = 0
 *
 * @param[in] a,b,c -
 *
 * @return the smallest positive number of two, or 0 if all 0 or negative..
 *
 * @author Niels Bassler
 */
static inline double _quadratic_solver(double a, double b, double c) {

    double d;
    double _d;
    double r1, r2;

    if (fabs(a) < OSH_GEMCA_SMALL) {
        if (fabs(b) > OSH_GEMCA_SMALL) {
            r1 = -c / b;
            if (r1 > 0.0)
                return r1;
            else
                return OSH_GEMCA_INFINITY;
        }
        return OSH_GEMCA_INFINITY;
    }

    d = b * b - 4.0 * a * c; /* discriminant */

    if (d < 0.0) {
        return OSH_GEMCA_INFINITY; /* no real roots */
    }

    _d = sqrt(d);
    r1 = (-b + _d) / (2.0 * a);
    r2 = (-b - _d) / (2.0 * a);
    return _minpos(r1, r2);
}

/**
 * @brief return the smallest positive number of two, or 0 if all are negative.
 *
 * @param[in] a,b - two doubles to be checked
 *
 * @return the smallest positive number of two, or 0 if all 0 or negative..
 *
 * @author Niels Bassler
 */
static inline double _minpos(double a, double b) {

    if (a > 0.0) {
        if (b > 0.0) {
            return MIN(a, b);
        } else {
            return a;
        }
    } else { /* a is 0 or not positve */
        if (b > 0.0) {
            return b;
        }
    }
    return 0.0;
}
