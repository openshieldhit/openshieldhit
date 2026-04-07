#include "gemca/osh_gemca2_calc_zone.h"

#include <stdlib.h>

#include "common/osh_coord.h"
#include "common/osh_logger.h"
#include "common/osh_vect.h"
#include "gemca/osh_gemca2.h"
#include "gemca/osh_gemca2_calc_surface.h"
#include "gemca/osh_gemca2_defines.h"
#include "transport/osh_transport.h"

static inline int _in_zone(struct zone const *z, struct ray const *r);
static inline int _in_node(struct cgnode const *self, struct ray const *r);
static inline int _in_body(struct body const *b, struct ray const *r);
static inline enum osh_status _transform_to_local(struct body const *b, struct ray const *r, struct ray *tr);

/*
   TODO: Recursive evaluation of the AST can become computationally expensive, especially for complex geometries.
   Consider optimizations such as bounding volume hierarchies (BVHs) or spatial partitioning to quickly exclude
   large portions of geometry from detailed evaluation.
 */

/**
 * @brief For a given ray, check what zone we are in.
 *
 * @param[in] g - a gemca object
 * @param[in] r - a ray
 *
 * @returns Dense internal zone index, or OSH_GEMCA_ZONE_INDEX_INVALID if no zone contains the ray.
 *
 * @author Niels Bassler
 */
size_t osh_gemca_get_zone_index(struct gemca_workspace *g, struct ray *r) {

    size_t i;

    if (!g || !r) {
        return OSH_GEMCA_ZONE_INDEX_INVALID;
    }

    for (i = 0; i < g->nzones; i++) {
        if (_in_zone(g->zones[i], r)) {
            return i;
        }
    }
    return OSH_GEMCA_ZONE_INDEX_INVALID;
}

/**
 * @brief Check if ray is in this zone.
 *
 * @details This function is recursive and will traverse the AST of the zone.
 *
 * @param[in] z - a zone
 * @param[out] r - a ray
 *
 * @returns 1 if ray is in zone, 0 if not.
 *
 * @author Niels Bassler
 */
static inline int _in_zone(struct zone const *z, struct ray const *r) {

    return _in_node(&z->node, r);
}

/**
 * @brief Check if ray is in this node.
 *
 * @param[in] self
 * @param[in] r
 * @return int
 */
static inline int _in_node(struct cgnode const *self, struct ray const *r) {

    int op;

    if (self->type == _OSH_GEMCA_CGNODE_BODY) {
        return _in_body(self->body, r);
    }

    op = (unsigned char) self->op;

    switch (op) { /* TODO: use defines instead of checking on char */

    case '+': /* intersection: both must be inside */
        return _in_node(self->left, r) && _in_node(self->right, r);

    case '-': /* difference: inside left, outside right */
        return _in_node(self->left, r) && !_in_node(self->right, r);

    case '|': /* union: either side suffices */
        return _in_node(self->left, r) || _in_node(self->right, r);

    default:
        osh_error("_in_node(): unknown operator");
        return 0;
    }
}

/**
 * @brief Check if ray is in this body.
 *
 * @details Leaf nodes of the AST are bodies. This function checks if a ray is inside a body.
 *
 * @param[in] - body
 * @param[in] - ray
 *
 * @returns 1 if ray is inside body, 0 if not.
 *
 * @author Niels Bassler
 */
inline int _in_body(struct body const *b, struct ray const *r) {

    int i;
    struct ray tr; /* ray in body-coordinate system */

    if (_transform_to_local(b, r, &tr) != OSH_OK) {
        return 0;
    }

    for (i = 0; i < b->nsurfs; i++) {
        /* see if we are on the good or bad side of the surface */
        if (!(osh_gemca2_check_surface_side(b->surfs[i], &tr))) {
            return 0;
        }
    }
    return 1;
}

/**
 * @brief Transform ray according to surface type and its coordinates.
 *
 * @details This function is used to transform a ray from OSH_COORD_UNIVERSE to the local coordinate system of a body.
 *
 * @param[in] b - body parameters incl. its transformation matrix
 * @param[in] r - input ray in OSH_COORD_UNIVERSE
 * @param[out] tr - transformed output ray in system given by b->coord
 *
 * @returns OSH_OK on success, OSH_ENOTSUP if the coordinate system is not supported.
 *
 * @author Niels Bassler
 */
static inline enum osh_status _transform_to_local(struct body const *b, struct ray const *r, struct ray *tr) {

    int i;
    int j;

    /* TODO: For now, just copy all elements of the ray. Later this can be optimized. */
    for (i = 0; i < 3; i++) {
        tr->p[i] = r->p[i];
        tr->cp[i] = r->cp[i];
    }
    tr->system = (unsigned char) b->coord;

    /* then overwrite the values which may change: */
    switch (b->coord) {
    case OSH_COORD_UNIVERSE:
        break;

    case OSH_COORD_BCALIGN:
        /* simple translation */
        for (i = 0; i < 3; i++) {
            j = i * 4;
            tr->p[i] = r->p[i] + b->t[j + 3]; /* notice, that in osh_coord.h see comment */
            tr->cp[i] = r->cp[i];
        }
        break;

    case OSH_COORD_BZALIGN:
        /* simple translation and rotation, so we have to use osh_coord_trans_ray */
        osh_coord_trans_ray_r(r, tr, b->t);
        break;

    default:
        osh_error("_transform_to_local() unsupported coordinate system :%i", b->coord);
        return OSH_ENOTSUP;
    }
    return OSH_OK;
}
