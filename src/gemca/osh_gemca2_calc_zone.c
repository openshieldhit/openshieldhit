#include "gemca/osh_gemca2_calc_zone.h"

#include <stdio.h>
#include <stdlib.h> // exit()

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
static inline int _transform_to_local(struct body const *b, struct ray const *r, struct ray *tr);

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
 * @returns The zone number we are in.
 *
 * @author Niels Bassler
 */
size_t osh_gemca_get_zone(struct gemca_workspace *g, struct ray *r) {

    size_t i;

    for (i = 0; i < g->nzones; i++) {
        // printf("\n --- _get_zone(), test zone %li '%s' ----------------- \n", g->zones[i]->id, g->zones[i]->name);
        if (_in_zone(g->zones[i], r))
            return g->zones[i]->id;
    }
    return 0; // TODO, -1 for invalid
}

/**
 * @brief For a given ray, check what zone we are in.
 *
 * @param[in] g - a gemca object
 * @param[in] r - a ray
 *
 * @returns The zone number we are in.
 *
 * @author Niels Bassler
 */
size_t osh_gemca_get_zone_index(struct gemca_workspace *g, struct ray *r) {

    size_t i;

    for (i = 0; i < g->nzones; i++) {
        // printf("\n --- _get_zone(), test zone %li '%s' ----------------- \n", g->zones[i]->id, g->zones[i]->name);
        if (_in_zone(g->zones[i], r))
            return i;
    }
    return 0; // TODO, -1 for invalid
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

    int a;
    int b;
    int op;

    if (self->type == _OSH_GEMCA_CGNODE_BODY) {
        return _in_body(self->body, r);
    } else {
        a = _in_node(self->left, r);
        b = _in_node(self->right, r);
        op = self->op;

        switch (op) { // TODO: use defines instead of checking on char

        case '+':
            if (a && b)
                return 1;
            break;

        case '-':
            if (a && !b)
                return 1;
            break;

        case '|':
            if (a || b)
                return 1;
            break;

        default:
            osh_error(EX_SOFTWARE, "_in_node(): unknown operator");
            break;
        }
    }
    return 0;
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

    _transform_to_local(b, r, &tr);

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
 * @returns 1
 *
 * @author Niels Bassler
 */
static inline int _transform_to_local(struct body const *b, struct ray const *r, struct ray *tr) {

    int i;
    int j;

    // TODO: For now, just copy all elements of the ray. Later this can be optimized.
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
