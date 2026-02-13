
#include "gemca/osh_gemca2_calc_surface.h"

#include <stdio.h>
#include <stdlib.h>

#include "common/osh_coord.h"
#include "common/osh_exit.h"
#include "common/osh_logger.h"
#include "common/osh_vect.h"
#include "gemca/osh_gemca2.h"
#include "gemca/osh_gemca2_defines.h"
#include "gemca/osh_gemca2_dist.h"

static int _inside_sphere(struct surface const *sf, struct ray const *r);
static int _inside_ellipsoid(struct surface const *sf, struct ray const *r);
static int _inside_cylz(struct surface const *sf, struct ray const *r);
static int _inside_ellz(struct surface const *sf, struct ray const *r);
static int _inside_cone(struct surface const *sf, struct ray const *r);
static int _inside_plane(struct surface const *sf, struct ray const *r);
static int _inside_plane_xyz(int axis, struct surface const *sf, struct ray const *r);

/**
 * @brief Allocate memory for n surfaces in the given body b, initialized to zero.
 *
 * @details
 *
 * @param[in,out] b - a body were the surfaces will be memory allocated
 * @param[in] n - the number of sufraces which will be allocated
 *
 * @returns 1
 *
 * @author Niels Bassler
 */
int osh_gemca2_add_surfaces(struct body *b, int n) {
    int i;

    b->nsurfs = n;
    b->surfs = calloc(n, sizeof(struct surface *));

    for (i = 0; i < n; i++) {
        b->surfs[i] = calloc(1, sizeof(struct surface));
    }

    return 1;
}

/**
 * @brief Set the number of required paramters in s for a given surface type. Allocate memory for parameters.
 *
 * @param[in,out] s - a given surface
 * @param[in] type - a surface type of OSH_GEMCA_SURF_* as given in osh_gemca2_defines.h
 *
 * @returns 1
 *
 * @author Niels Bassler
 */
int osh_gemca2_add_surf_pars(struct surface *s, int type) {

    switch (type) {
    case OSH_GEMCA_SURF_PLANEX:
        s->np = 2; /* [A,D]   Ax + D = 0 */
        break;
    case OSH_GEMCA_SURF_PLANEY:
        s->np = 2; /* [B,D]   By + D = 0 */
        break;
    case OSH_GEMCA_SURF_PLANEZ:
        s->np = 2; /* [C,D]   Cz + D = 0 */
        break;
    case OSH_GEMCA_SURF_PLANE:
        s->np = 4; /* [A,B,C,D]   Ax + By + Cz + D = 0 */
        break;
    case OSH_GEMCA_SURF_SPHERE:
        s->np = 1; /* [R^2]   x^2 + y^2 + z^2 - R^2 = 0  (center point is always at 0,0,0) */
        break;
    case OSH_GEMCA_SURF_ELLIPSOID:
        s->np = 3; /* [A,B,C]^2   x^2/A^2 + y^2/B^2 + z^2/C^2 - 1 = 0  (center point is always at 0,0,0) */
        break;
    case OSH_GEMCA_SURF_CYLZ: /* cylinder, infinite along z */
        s->np = 1;            /* [R^2]   x^2 + y^2 * R^2 = 0  (center point is always at 0,0) in X,Y only */
        break;
    case OSH_GEMCA_SURF_ELLZ: /* elliptical cylinder, infinite along z */
        s->np = 2;            /* [A,B]^2   x^2/A^2 + y^2/B^2 - 1 = 0 (center point is always at 0,0) in X,Y only */
        break;
    case OSH_GEMCA_SURF_CONE: /* cone , infinite along z */
        s->np = 2;            /* [A,B]^2   x^2 + y^2 - (z-A)^2/B^2 = 0 (always along z) */
        break;
    default:
        osh_error(EX_CONFIG, "_add_surf_pars: unknown surface type: %i", type);
        break;
    }
    s->p = calloc(s->np, sizeof(double));
    s->type = type;
    return 1;
}

/**
 * @brief Checks if we are on a "negative" side (= "inside") of a surface.
 *
 * @details All surfaces are setup with a positive and negative side.
 *          The postive side will always face out of a body, so if the signed distance is positve
 *          it means the ray is inside the body.
 *
 * @param[in] sf - a surface to be checked
 * @param[out] r - a transformed ray which will be checked (must be in body->system coords).
 *
 * @returns 1 if inside the surface, 0 if ouside the surface
 *
 * @author Niels Bassler
 */
int osh_gemca2_check_surface_side(struct surface const *sf, struct ray const *r) {

    switch (sf->type) {
    case OSH_GEMCA_SURF_SPHERE:
        return _inside_sphere(sf, r);
        break;

    case OSH_GEMCA_SURF_ELLIPSOID:
        return _inside_ellipsoid(sf, r);
        break;

    case OSH_GEMCA_SURF_CYLZ:
        return _inside_cylz(sf, r);
        break;

    case OSH_GEMCA_SURF_ELLZ:
        return _inside_ellz(sf, r);
        break;

    case OSH_GEMCA_SURF_CONE:
        return _inside_cone(sf, r);
        break;

    case OSH_GEMCA_SURF_PLANEX:
        return _inside_plane_xyz(0, sf, r);
        break;

    case OSH_GEMCA_SURF_PLANEY:
        return _inside_plane_xyz(1, sf, r);
        break;

    case OSH_GEMCA_SURF_PLANEZ:
        return _inside_plane_xyz(2, sf, r);
        break;

    case OSH_GEMCA_SURF_PLANE:
        return _inside_plane(sf, r);
        break;

    case OSH_GEMCA_SURF_NONE:
        osh_error(EX_SOFTWARE, "_check_surface_side(): NONE surface type %i\n", sf->type);
        break;

    default:
        osh_error(EX_SOFTWARE, "_check_surface_side(): unknown surface type %i\n", sf->type);
        break;
    }
    return 1; /* is inside */
}

static int _inside_sphere(struct surface const *sf, struct ray const *r) {
    double d;
    /* outside if x^2 + y^2 + z^2 > r^2 */

    d = osh_vect_len2(r->p) - sf->p[0];
    if (d > OSH_GEMCA_SMALL) {
        return 0; /* point is outside of sphere */
    } else if (d < -OSH_GEMCA_SMALL) {
        return 1; /* point is inside of sphere */
    } else {
        /* handle the case if we are on the surface */
        if (osh_vect_dot(r->p, r->cp) < 0.0) {
            return 1; /* is traveling into the sphere */
        } else {
            return 0; /* is traveling out of the sphere or is tangent*/
        }
    }
}

static int _inside_ellipsoid(struct surface const *sf, struct ray const *r) {
    double d;
    int i;

    d = 0.0;
    /* Calculate "distance" value (sum of squared distances divided by semi-axis lengths) */
    for (i = 0; i < 3; i++) {
        d += (r->p[i] * r->p[i]) / sf->p[i];
    }

    if (d > 1.0 + OSH_GEMCA_SMALL) {
        /* Point is outside the ellipsoid */
        return 0;
    } else if (d < 1.0 - OSH_GEMCA_SMALL) {
        /* Point is inside the ellipsoid */
        return 1;
    } else {
        /* Handle edge case: Point is on the ellipsoid surface */
        d = 0.0;
        for (i = 0; i < 3; i++) {
            d += (r->p[i] / sf->p[i]) * r->cp[i];
        }

        if (d < 0.0) {
            // Ray is traveling into the ellipsoid
            return 1;
        } else {
            // Ray is traveling out of the ellipsoid or tangent
            return 0;
        }
    }
}

static int _inside_cylz(struct surface const *sf, struct ray const *r) {
    double d;
    double dot_product;

    /* Calculate distance-like value for the cylinder in the XY-plane */
    d = r->p[0] * r->p[0] + r->p[1] * r->p[1] - sf->p[0];

    if (d > OSH_GEMCA_SMALL) {
        /* Point is outside the cylinder */
        return 0;
    } else if (d < -OSH_GEMCA_SMALL) {
        /* Point is inside the cylinder */
        return 1;
    } else {
        /* Handle edge case: Point is on the surface of the cylinder */
        dot_product = r->p[0] * r->cp[0] + r->p[1] * r->cp[1];

        if (dot_product < 0.0) {
            /* Ray is traveling into the cylinder */
            return 1;
        } else {
            /* Ray is traveling out of the cylinder or tangent */
            return 0;
        }
    }
}

static int _inside_ellz(struct surface const *sf, struct ray const *r) {
    double d;

    /* Calculate the "distance-like" value for the elliptical cylinder */
    d = (r->p[0] * r->p[0] / sf->p[0]) + (r->p[1] * r->p[1] / sf->p[1]) - 1.0;

    if (d > OSH_GEMCA_SMALL) {
        /* Point is outside the elliptical cylinder */
        return 0;
    } else if (d < -OSH_GEMCA_SMALL) {
        /* Point is inside the elliptical cylinder */
        return 1;
    } else {
        /* Handle edge case: Point is on the surface of the elliptical cylinder */
        d = (r->p[0] / sf->p[0]) * r->cp[0] + (r->p[1] / sf->p[1]) * r->cp[1];

        if (d < 0.0) {
            /* Ray is traveling into the elliptical cylinder */
            return 1;
        } else {
            /* Ray is traveling out of the elliptical cylinder or tangent */
            return 0;
        }
    }
}

static int _inside_cone(struct surface const *sf, struct ray const *r) {
    double d;

    /* Distance-like value for the conical surface */
    d = (r->p[0] * r->p[0]) + (r->p[1] * r->p[1]); /* x^2 + y^2 */
    d -= sf->p[1] * (r->p[2] * r->p[2]);           /* - slope * z^2 */

    if (d > OSH_GEMCA_SMALL) {
        /* Point is outside the cone's lateral surface */
        return 0;
    } else if (d < -OSH_GEMCA_SMALL) {
        /* Point is inside the cone's lateral surface */
        return 1;
    } else {
        /* Handle edge case: Point is on the lateral surface of the cone */
        d = (r->p[0] * r->cp[0]) + (r->p[1] * r->cp[1]) - (sf->p[1] * r->p[2] * r->cp[2]);

        if (d < 0.0) {
            /* Ray is traveling into the cone */
            return 1;
        } else {
            /* Ray is traveling out of the cone or tangent */
            return 0;
        }
    }
}

static int _inside_plane(struct surface const *sf, struct ray const *r) {
    double d;

    d = osh_vect_dot(sf->p, r->p) + sf->p[3];

    if (d > OSH_GEMCA_SMALL) {
        return 0;
    } else if (d < -OSH_GEMCA_SMALL) {
        return 1;
    } else {
        /* we need to decide by direction */
        if (osh_vect_dot(sf->p, r->cp) > 0.0) {
            return 0; /* is traveling out of the surface, so it is outside */
        } else {
            return 1; /* is traveling into the surface, or exactly parallel on the surface, so it is inside */
        }
    }
}

static int _inside_plane_xyz(int axis, struct surface const *sf, struct ray const *r) {
    double d;

    d = sf->p[0] * r->p[axis] + sf->p[1];

    if (d > OSH_GEMCA_SMALL) {
        return 0;
    } else if (d < -OSH_GEMCA_SMALL) {
        return 1;
    } else {
        /* we need to decide by direction */
        if (sf->p[0] * r->cp[axis] > 0.0) {
            return 0; /* is traveling out of the surface, so it is outside */
        } else {
            return 1; /* is traveling into the surface, or exactly parallel on the surface, so it is inside */
        }
    }
}