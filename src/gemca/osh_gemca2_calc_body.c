#include "gemca/osh_gemca2_calc_body.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "common/osh_const.h"
#include "common/osh_coord.h"
#include "common/osh_logger.h"
#include "common/osh_vect.h"
#include "gemca/osh_gemca2.h"
#include "gemca/osh_gemca2_calc_surface.h"
#include "gemca/osh_gemca2_defines.h"

static int setup_body(struct body *b);

static int _setup_sph(struct body *b);
static int _setup_wed(struct body *b);
static int _setup_arb(struct body *b);

static int _setup_box(struct body *b);
static int _setup_vox(struct body *b);
static int _setup_rpp(struct body *b);

static int _setup_rcc(struct body *b);
static int _setup_rec(struct body *b);
static int _setup_trc(struct body *b);
static int _setup_ell(struct body *b);

static int _setup_yzp(struct body *b);
static int _setup_xzp(struct body *b);
static int _setup_xyp(struct body *b);
static int _setup_pla(struct body *b);

static int _setup_rot(struct body *b);
static int _setup_cpy(struct body *b);
static int _setup_mov(struct body *b);

static void _vertex_index_arb_fluka(double d, int *i);

/**
 * @brief Setup all bodies in a gemca object.
 *
 * @details This is the only function exposed outside of this file.
 *
 * @param[in] - a gemca object
 * @param[out] - a gemca object
 *
 * @returns 1
 *
 * @author Niels Bassler
 */
int osh_gemca_body_setup(struct gemca_workspace *g) {

    size_t i;

    for (i = 0; i < g->nbodies; i++) {
        setup_body(g->bodies[i]);
    }
    return 1;
}

/**
 * @brief Setup the parametric description of a body.
 *
 * @details This function will setup the body according to its type.
 *
 * @param[in] b - the body which will be setup
 * @param[out] b - the body which will be setup
 *
 * @returns 1
 *
 * @author Niels Bassler
 */
static int setup_body(struct body *b) {

    switch (b->type) {
    case OSH_GEMCA_BODY_SPH:
        _setup_sph(b);
        break;
    case OSH_GEMCA_BODY_WED:
        _setup_wed(b);
        break;
    case OSH_GEMCA_BODY_ARB:
        _setup_arb(b);
        break;

    case OSH_GEMCA_BODY_BOX:
        _setup_box(b);
        break;
    case OSH_GEMCA_BODY_VOX:
        _setup_vox(b);
        break;
    case OSH_GEMCA_BODY_RPP:
        _setup_rpp(b);
        break;

    case OSH_GEMCA_BODY_RCC:
        _setup_rcc(b);
        break;
    case OSH_GEMCA_BODY_REC:
        _setup_rec(b);
        break;
    case OSH_GEMCA_BODY_TRC:
        _setup_trc(b);
        break;
    case OSH_GEMCA_BODY_ELL:
        _setup_ell(b);
        break;

    case OSH_GEMCA_BODY_YZP:
        _setup_yzp(b);
        break;
    case OSH_GEMCA_BODY_XYP:
        _setup_xyp(b);
        break;
    case OSH_GEMCA_BODY_XZP:
        _setup_xzp(b);
        break;
    case OSH_GEMCA_BODY_PLA:
        _setup_pla(b);
        break;

    case OSH_GEMCA_BODY_ROT:
        _setup_rot(b);
        break;
    case OSH_GEMCA_BODY_CPY:
        _setup_cpy(b);
        break;
    case OSH_GEMCA_BODY_MOV:
        _setup_mov(b);
        break;

    default:
        break;
    }
    return 1;
}

/**
 * @brief Setup a sphere.
 *
 * @details
 *          SPH given user args b->a[*]: (all in OSH_COORD_UNIVERSE)
 *          0,1,2: center coordinates of the sphere in OSH_COORD_UNIVERSE.
 *          3: radius
 *
 * @param[in,out] b - the body which will be setup
 *
 * @returns 1
 *
 * @author Niels Bassler
 */
static int _setup_sph(struct body *b) {
    int const nsurfs = 1;
    struct surface *sf;

    /* ----------- Setup translation matrix */
    b->coord = OSH_COORD_BCALIGN;
    b->t[3] = -b->a[0];
    b->t[7] = -b->a[1];
    b->t[11] = -b->a[2];

    /* ----------- Setup surfaces */
    osh_gemca2_add_surfaces(b, nsurfs);

    sf = b->surfs[0];
    osh_gemca2_add_surf_pars(sf, OSH_GEMCA_SURF_SPHERE);
    sf->p[0] = b->a[3] * b->a[3]; /* radius squared for computation speed */

    return 1;
}

/**
 * @brief Setup a wedge.
 *
 * @details
 *          WED given user args b->a[*]: (all in OSH_COORD_UNIVERSE)
 *          0,1,2: r0 - origin point
 *          3,4,5: r1 - spanning vector1 (height)
 *          6,7,8: r2 - spanning vector2
 *          9,10,11: r3 - spanning vector3
 *          Note: r2 and r3 must be orthogonal to r1
 *
 * @param[in,out] b - the body which will be setup
 *
 * @returns 1
 *
 * @author Niels Bassler
 */
static int _setup_wed(struct body *b) {

    int const nsurfs = 5;
    struct surface *sf; /* temporary surface */

    double v[3]; /* temporary vector */
    double u[3]; /* temporary vector */
    double p[3]; /* temporary point */
    double c[3]; /* some point guaranteed to be inside the wedge */

    double *r0, *r1, *r2, *r3; /* alias pointers to make code more readable */
    int i;

    r0 = &(b->a[0]); /* origin point */
    r1 = &(b->a[3]); /* spanning vector1 (height) */
    r3 = &(b->a[9]); /* spanning vector2 */
    r2 = &(b->a[6]); /* spanning vector3 */

    /* ----------- Setup translation matrix */
    b->coord = OSH_COORD_UNIVERSE; /* no translation matrix needed */

    /* ----------- Setup surfaces */
    /* choose a point inside the wedge */
    /*
       (X,Y,Z) is coordinate of a reference point of WED which lies:
       at a 1/4 wedge height counting from the plane give by R2 and R3
       and a 1/4 length of R2+R3 from vector R1
     */
    for (i = 0; i < 3; i++) {
        c[i] = b->a[0] + 0.25 * (b->a[3 + i] + b->a[6 + i] + b->a[9 + i]);
    }

    /* XYZ planes are always setup so the normal vector points out of the volume of interest */

    /* bottom of the wedge: point R0 with -R1 as normal vector */
    osh_gemca2_add_surfaces(b, nsurfs);
    osh_vect_reverse(r1, v); /* v points out of the volume */
    sf = b->surfs[0];
    osh_gemca2_add_surf_pars(sf, OSH_GEMCA_SURF_PLANE);
    osh_vect_eqpln(r0, v, sf->p); /* surface parameters, i.e. Ax + By + Cz + D = 0; */

    /* find a point in the opposite plane which is parallel to the base of the wedge */
    sf = b->surfs[1];
    osh_gemca2_add_surf_pars(sf, OSH_GEMCA_SURF_PLANE);
    osh_vect_add(r0, r1, v);
    osh_vect_eqpln(r0, v, sf->p);

    /* r1 x r2 gives a normal vector to the plane spanned by those two vectors */
    sf = b->surfs[2];
    osh_gemca2_add_surf_pars(sf, OSH_GEMCA_SURF_PLANE);
    osh_vect_cross(r1, r2, v); /* v is now a normal vector to the plane of r1 and r2 */
    osh_vect_sub(r0, c, u);    /* u is now a vector from a point in the plane to a point in the center */
    if (osh_vect_dot(v, u) > 0) {
        /* normal vector is pointing into the volume */
        osh_vect_reverse(v, u); /* let normal vector point out of the volume */
        osh_vect_eqpln(r0, v, sf->p);
    } else {
        osh_vect_eqpln(r0, u, sf->p);
    }

    /* r1 x r3 gives a normal vector to the plane spanned by those two vectors */
    sf = b->surfs[3];
    osh_gemca2_add_surf_pars(sf, OSH_GEMCA_SURF_PLANE);
    osh_vect_cross(r3, r1, v);
    osh_vect_sub(r0, c, u); /* u is now a vector from a point in the plane to a point in the center */
    if (osh_vect_dot(v, u) > 0) {
        osh_vect_reverse(v, u); /* let normal vector point out of the volume */
        osh_vect_eqpln(r0, v, sf->p);
    } else {
        osh_vect_eqpln(r0, u, sf->p);
    }

    sf = b->surfs[4];
    osh_gemca2_add_surf_pars(sf, OSH_GEMCA_SURF_PLANE);
    osh_vect_add(r0, r2, p);  /* find a point in the last plane */
    osh_vect_sub(r1, r2, v);  /* find a vector inside the last plane */
    osh_vect_cross(r0, v, u); /* u is now a normal vector to the last plane */
    osh_vect_eqpln(p, u, sf->p);
    osh_vect_sub(p, c, v); /* v is now a vector from a point in the plane to a point in the center */
    if (osh_vect_dot(v, u) > 0) {
        osh_vect_reverse(u, v); /* let normal vector point out of the volume */
        osh_vect_eqpln(p, u, sf->p);
    } else {
        osh_vect_eqpln(p, v, sf->p);
    }
    return 1;
}

/**
 * @brief Setup a arbitrary convex polyhedron.
 *
 * @details ARB is poorly defined, old SHIELD-HIT seems to be different from FLUKA:
 *          This solution tries to encompass both formats.
 *          Notice also that the SHIELD-HIT solution does not allow for vertices < 8.
 *          SHIELD-HIT: Points are specified in a certain order as given in the figure in the maunal.
 *          FLUKA: An extra card specifying the surfaces (called "faces" in the FLUKA manual)
 *                 Notice also that FLUKA allows 4 vertices to specify a "face", even if 3 vertices are enough.
 *
 *          ARB given user args b->a[*]: (all in OSH_COORD_UNIVERSE)
 *          0,1,2: p[0] - starting point of box
 *          3,4,5: p[1] - next point
 *          ...
 *          21,22,23: p[7] - last point
 *          24,25,26:      - optional vertice index (FLUKA FORMAT only)
 *
 * @param[in,out] b - the body which will be setup
 *
 * @returns 1
 *
 * @author Niels Bassler
 */

static int _setup_arb(struct body *b) {

    int const nsurfs = 6;
    struct surface *sf;      /* temporary surface */
    double *p[8];            /* 8 points describing the arb */
    double u[3], v[3], w[3]; /* temporary vectors */
    int i, j;
    int k[3];

    /* setup alias for the vertices of ARB for easier to read code */
    for (i = 0; i < 8; i++) {
        p[i] = &(b->a[i * 3]); /* list of vertices */
    }

    /* ----------- Setup translation matrix */
    b->coord = OSH_COORD_UNIVERSE; /* no translation matrix needed */

    /* ----------- Setup surfaces */
    for (i = 0; i < nsurfs; i++) {
        sf = b->surfs[i];
        osh_gemca2_add_surf_pars(sf, OSH_GEMCA_SURF_PLANE);
    }

    /* A normal arb is specified by 8 points, 8 * 3 = 24 arguemnts.
       FLUKA format has a 25th argument for the faces of the surfaces */
    if (b->na == 25) {
        /* this is FLUKA format ARB. */
        /* face descrpitor is stored in argument 25, i.e. b->a[24] */
        _vertex_index_arb_fluka(b->a[24], k); /* get the first two indices from face-descriptor */
        for (j = 24; j < 30; j++) {           /* loop over each face-descriptor */
            sf = b->surfs[0];
            osh_vect_sub(p[k[0]], p[k[1]], u); /* build vector 1 */
            osh_vect_sub(p[k[0]], p[k[2]], v); /* build vector 2 */
            osh_vect_cross(v, u, w);           /* cross product to get the normal vector */
            // TODO: check if w points outside or into the body, both is allowed in FLUKA
            osh_vect_eqpln(p[k[0]], w, sf->p); /* build the equation which describes this surface */
        }
    } else {
        /* this is SH12A/MORSE format ARB. */
        sf = b->surfs[0];
        osh_vect_sub(p[0], p[1], u);
        osh_vect_sub(p[1], p[2], v);
        osh_vect_cross(v, u, w);
        osh_vect_eqpln(p[0], w, sf->p);

        sf = b->surfs[1];
        osh_vect_sub(p[1], p[2], u);
        osh_vect_sub(p[1], p[5], v);
        osh_vect_cross(v, u, w);
        osh_vect_eqpln(p[1], w, sf->p);

        sf = b->surfs[2];
        osh_vect_sub(p[2], p[3], u);
        osh_vect_sub(p[6], p[2], v);
        osh_vect_cross(v, u, w);
        osh_vect_eqpln(p[2], w, sf->p);

        sf = b->surfs[3];
        osh_vect_sub(p[0], p[3], u);
        osh_vect_sub(p[0], p[4], v);
        osh_vect_cross(v, u, w);
        osh_vect_eqpln(p[0], w, sf->p);

        sf = b->surfs[4];
        osh_vect_sub(p[0], p[4], u);
        osh_vect_sub(p[0], p[1], v);
        osh_vect_cross(v, u, w);
        osh_vect_eqpln(p[0], w, sf->p);

        sf = b->surfs[5];
        osh_vect_sub(p[4], p[7], u);
        osh_vect_sub(p[4], p[5], v);
        osh_vect_cross(v, u, w);
        osh_vect_eqpln(p[4], w, sf->p);
    }
    return 1;
}

/**
 * @brief Setup a box geometry with any orientation.
 *
 * @details
 *          BOX given user args a[*]: (all in OSH_COORD_UNIVERSE)
 *          0,1,2: p - starting point of box
 *          3,4,5: r - spanning vector 1
 *          6,7,8: s - spanning vector 2
 *          9,10,11: t - spanning vector 3
 *          Note: r,s,t must be orthorgonal.
 *
 * @param[in,out] b - the body which will be setup
 *
 * @returns 1
 *
 * @author Niels Bassler
 */
static int _setup_box(struct body *b) {

    int const nsurfs = 6;
    struct surface *sf; /* temporary surface */

    double *p, *r, *s, *t;
    double u[3], v[3], w[3];
    int i;

    p = &(b->a[0]); /* starting point of box */
    r = &(b->a[3]); /* spanning vector 1 */
    s = &(b->a[6]); /* spanning vector 2 */
    t = &(b->a[9]); /* spanning vector 3 */

    /* ----------- Setup translation matrix */
    b->coord = OSH_COORD_UNIVERSE; /* no translation matrix needed */

    /* ----------- Setup surfaces */
    for (i = 0; i < nsurfs; i++) {
        sf = b->surfs[i];
        osh_gemca2_add_surf_pars(sf, OSH_GEMCA_SURF_PLANE);
    }

    /* Setting up a BOX is simple, since the vectors are normal vectors pointing out of the box seen
       from the given vertex. The other three planes are spanned from the opposite side. */
    sf = b->surfs[0];
    osh_vect_reverse(r, u); /* reverse the vectors, so they point into the box */
    osh_vect_eqpln(p, u, sf->p);

    sf = b->surfs[1];
    osh_vect_reverse(s, u); /* reverse the vectors, so they point into the box */
    osh_vect_eqpln(p, u, sf->p);

    sf = b->surfs[2];
    osh_vect_reverse(t, u); /* reverse the vectors, so they point into the box */
    osh_vect_eqpln(p, u, sf->p);

    /* Find opposite vertex to cube new point */
    osh_vect_add(p, r, u);
    osh_vect_add(u, s, v);
    osh_vect_add(v, t, w); /* w is now the oppoiste to point p of the box */

    sf = b->surfs[3];
    osh_vect_eqpln(w, u, sf->p);

    sf = b->surfs[4];
    osh_vect_eqpln(w, u, sf->p);

    sf = b->surfs[5];
    osh_vect_eqpln(w, u, sf->p);

    return 1;
}

/**
 * @brief Setup a CT voxel geometry.
 *
 * @details
 *          VOX given user args b->a[*]: (all in OSH_COORD_UNIVERSE)
 *          0,1,2: isocenter x,y,z in [cm]
 *          3: couch [degrees]
 *          4: gantry [degrees]
 *          5: target dose in [Gy]
 *          remaining parameters will be read from the .hed file.
 *
 * @param[in,out] b - the body which will be setup
 *
 * @returns 1
 *
 * @author Niels Bassler
 */
static int _setup_vox(struct body *b) {

    int const nsurfs = 6;
    struct surface *sf; /* temporary surface */
    /* define unit vectors along the axes */
    double tb[3][3] = {
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1},
    };

    double couch_angle = 0.0;
    double gantry_angle = 0.0;

    int i, j;

    // TODO:
    // open b->filename_vox
    // get parameters from it

    /* ----------- Setup translation matrix */
    /* translation and rotation needed lowest corner is at 0,0,0 */
    b->coord = OSH_COORD_BZALIGN;
    couch_angle = (b->a[3] / 180.0) * OSH_M_PI;
    gantry_angle = (b->a[4] / 180.0) * OSH_M_PI;

    /* rotate the base tb according to gantry / couch angles */
    for (i = 0; i < 3; i++) {
        osh_vect_rot_y(couch_angle, tb[i]);
        osh_vect_rot_z(gantry_angle, tb[i]);
    }

    /* copy transposed matrix into translation matrix */
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            b->t[j * 4 + i] = tb[j][i]; /* remember; t is a 4x4 matrix */
        }
    }

    /* isocenter shift goes into translation part of transformation matrix */
    b->t[3] = -b->a[0];
    b->t[7] = -b->a[1];
    b->t[11] = -b->a[2];

    /* ----------- Setup surfaces */
    osh_gemca2_add_surfaces(b, nsurfs);

    // TODO: get b->a[] parameters set properly from voxelload()

    sf = b->surfs[0];
    osh_gemca2_add_surf_pars(sf, OSH_GEMCA_SURF_PLANEX);
    sf->p[0] = -1;      /* = A */
    sf->p[1] = b->a[0]; /* = D = -Ax   =>   Ax - D = 0 */

    sf = b->surfs[1];
    osh_gemca2_add_surf_pars(sf, OSH_GEMCA_SURF_PLANEX);
    sf->p[0] = 1;        /* normal vector inverted, so it points into the zone */
    sf->p[1] = -b->a[1]; /* Ax - D = 0 */

    sf = b->surfs[2];
    osh_gemca2_add_surf_pars(sf, OSH_GEMCA_SURF_PLANEY);
    sf->p[0] = -1;
    sf->p[1] = b->a[2];

    sf = b->surfs[3];
    osh_gemca2_add_surf_pars(sf, OSH_GEMCA_SURF_PLANEY);
    sf->p[0] = 1;
    sf->p[1] = -b->a[3];

    sf = b->surfs[4];
    osh_gemca2_add_surf_pars(sf, OSH_GEMCA_SURF_PLANEZ);
    sf->p[0] = -1;
    sf->p[1] = b->a[4];

    sf = b->surfs[5];
    osh_gemca2_add_surf_pars(sf, OSH_GEMCA_SURF_PLANEZ);
    sf->p[0] = 1;
    sf->p[1] = -b->a[5];

    return 1;
}

/**
 * @brief Setup a box which is parallel to X,Y,Z axes.
 *
 * @details
 *          RPP given user args b->a[*]: (all in OSH_COORD_UNIVERSE)
 *          0,1: x-min x-max
 *          2,3: y-min y-max
 *          4,5: z-min z-max
 *
 * @param[in,out] b - the body which will be setup
 *
 * @returns 1
 *
 * @author Niels Bassler
 */
static int _setup_rpp(struct body *b) {

    int const nsurfs = 6;
    struct surface *sf; /* temporary surface */

    /* ----------- Setup translation matrix */
    b->coord = OSH_COORD_UNIVERSE; /* no translation matrix neededs */

    /* ----------- Setup surfaces */
    osh_gemca2_add_surfaces(b, nsurfs);

    sf = b->surfs[0];
    osh_gemca2_add_surf_pars(sf, OSH_GEMCA_SURF_PLANEX);
    sf->p[0] = -1;      /* = A */
    sf->p[1] = b->a[0]; /* = D = -Ax   =>   Ax - D = 0 */

    sf = b->surfs[1];
    osh_gemca2_add_surf_pars(sf, OSH_GEMCA_SURF_PLANEX);
    sf->p[0] = 1;        /* normal vector inverted, so it points into the zone */
    sf->p[1] = -b->a[1]; /* Ax - D = 0 */

    sf = b->surfs[2];
    osh_gemca2_add_surf_pars(sf, OSH_GEMCA_SURF_PLANEY);
    sf->p[0] = -1;
    sf->p[1] = b->a[2];

    sf = b->surfs[3];
    osh_gemca2_add_surf_pars(sf, OSH_GEMCA_SURF_PLANEY);
    sf->p[0] = 1;
    sf->p[1] = -b->a[3];

    sf = b->surfs[4];
    osh_gemca2_add_surf_pars(sf, OSH_GEMCA_SURF_PLANEZ);
    sf->p[0] = -1;
    sf->p[1] = b->a[4];

    sf = b->surfs[5];
    osh_gemca2_add_surf_pars(sf, OSH_GEMCA_SURF_PLANEZ);
    sf->p[0] = 1;
    sf->p[1] = -b->a[5];

    return 1;
}

/**
 * @brief Setup a cylinder.
 *
 * @details
 *          RCC given user args b->a[*]: (all in OSH_COORD_UNIVERSE)
 *          0,1,2: p - center of cylinder base
 *          3,4,5: r - vector spanning the cylinder from the cylinder base
 *          6: cylinder radius
 *
 * @param[in,out] b - the body which will be setup
 *
 * @returns 1
 *
 * @author Niels Bassler
 */
static int _setup_rcc(struct body *b) {

    int const nsurfs = 3;
    struct surface *sf;

    double *p, *r;

    /* Cylinder has a base point p, and is spanned by a user given vector t. */
    p = &(b->a[0]); /* center of cylinder base in OSH_COORD_UNIVERSE */
    r = &(b->a[3]); /* vector spanning the cylinder from the cylinder base in OSH_COORD_UNIVERSE. */

    /* ----------- Setup translation matrix */
    b->coord = OSH_COORD_BZALIGN;
    osh_vect_setup_tmatrix_bzalign(p, r, b->t);

    /* ----------- Setup surfaces */
    osh_gemca2_add_surfaces(b, nsurfs);

    /* base plane in BZALIGN is always at point 0,0,0  */
    sf = b->surfs[0];
    osh_gemca2_add_surf_pars(sf, OSH_GEMCA_SURF_PLANEZ);
    sf->p[0] = -1; /* C : normal vectors of surfaces usually point outwards of volumes */
    sf->p[1] = 0;  /* D = -C * z0 */

    /* top plane */
    sf = b->surfs[1];
    osh_gemca2_add_surf_pars(sf, OSH_GEMCA_SURF_PLANEZ);
    sf->p[0] = 1;                       /* C */
    sf->p[1] = -sqrt(osh_vect_len2(r)); /* D = -C * z0 . In BZALIGN system, this is simply the height of the cylinder */

    sf = b->surfs[2];
    osh_gemca2_add_surf_pars(sf, OSH_GEMCA_SURF_CYLZ);
    sf->p[0] = b->a[6] * b->a[6]; /* radius squared */

    return 1;
}

/**
 * @brief Setup a elliptical cylinder, freely located and oriented in space.
 *
 * @details
 *          REC given user args b->a[*]: (all in OSH_COORD_UNIVERSE)
 *          0,1,2: p - center of cylinder base
 *          3,4,5: r - height
 *          6,7,8: s - minor axis vector
 *          9,10,11: t - major axis vector
 *
 * @param[in,out] b - the body which will be setup
 *
 * @returns 1
 *
 * @author Niels Bassler
 */
static int _setup_rec(struct body *b) {
    int const nsurfs = 3;
    struct surface *sf; /* temporary surface */

    double *p, *r, *s, *t;
    double wt[3];

    p = &(b->a[0]); /* center of cylinder base */
    r = &(b->a[3]); /* height */
    s = &(b->a[6]); /* minor axis vector */
    t = &(b->a[9]); /* major axis vector */

    b->coord = OSH_COORD_BZALIGN;

    /* ----------- Setup translation matrix */
    osh_vect_setup_tmatrix_bzalign(p, r, b->t);

    /* ----------- Setup surfaces */
    osh_gemca2_add_surfaces(b, nsurfs);

    /* base plane */
    sf = b->surfs[0];
    osh_gemca2_add_surf_pars(sf, OSH_GEMCA_SURF_PLANEZ);
    sf->p[0] = -1;    /* C */
    sf->p[1] = -p[2]; /* D = -C * z0 */

    /* top plane */
    sf = b->surfs[1];
    osh_gemca2_add_surf_pars(sf, OSH_GEMCA_SURF_PLANEZ);
    osh_vect_add(p, r, wt); /* w is now the opposite point to p */
    sf->p[0] = 1;           /* C, but opposite so normal vector points out of the cone */
    sf->p[1] = -wt[2];      /* D = -C * z0 */

    sf = b->surfs[2];
    osh_gemca2_add_surf_pars(sf, OSH_GEMCA_SURF_ELLZ);
    sf->p[0] = osh_vect_len2(s); /* length^2 of minor axis */
    sf->p[1] = osh_vect_len2(t); /* length^2 of major axis */

    return 1;
}

/**
 * @brief Setup a truncated code.
 *
 * @details
 *          TRC given user args b->a[*]: (all in OSH_COORD_UNIVERSE)
 *          0,1,2: p - center point of cone base
 *          3,4,5: r - vector along cone center axis
 *          6: r1 - radius at cone base
 *          7: r2 - radius at cone top
 *          Note: r1 > r2.
 *
 * @param[in,out] b - the body which will be setup
 *
 * @returns 1
 *
 * @author Niels Bassler
 */
static int _setup_trc(struct body *b) {

    int const nsurfs = 3;
    struct surface *sf; /* temporary surface */

    double *p, *r;

    double wt[3];
    double h, r1, r2; /* scalar cylinder height, minor axis, major axis */

    p = &(b->a[0]); /* center point of cone base */
    r = &(b->a[3]); /* vector along cone center axis */
    r1 = (b->a[6]); /* radius at cone base */
    r2 = (b->a[7]); /* radius at cone top */

    /* ----------- Setup translation matrix */
    b->coord = OSH_COORD_BZALIGN;
    osh_vect_setup_tmatrix_bzalign(p, r, b->t);

    /* ----------- Setup surfaces */
    osh_gemca2_add_surfaces(b, nsurfs);

    /* base plane */
    sf = b->surfs[0];
    osh_gemca2_add_surf_pars(sf, OSH_GEMCA_SURF_PLANEZ);
    sf->p[0] = -1;    /* C */
    sf->p[1] = -p[2]; /* D = -C * z0 */

    /* top plane */
    sf = b->surfs[1];
    osh_gemca2_add_surf_pars(sf, OSH_GEMCA_SURF_PLANEZ);
    osh_vect_add(p, r, wt); /* w is now the opposite point to p */
    sf->p[0] = 1;           /* C, but opposite so normal vector points out of the cone */
    sf->p[1] = -wt[2];      /* D = -C * z0 */

    /* cone surface */
    sf = b->surfs[2];
    osh_gemca2_add_surf_pars(sf, OSH_GEMCA_SURF_CONE);
    h = sqrt(osh_vect_len2(r));                 /* height of the cone */
    sf->p[0] = h / (1.0 - r2 / r1);             /* slope factor */
    sf->p[1] = sf->p[0] * sf->p[0] / (r1 * r1); /* shape parameter */

    return 1;
}

/**
 * @brief Setup an ellipsoid.
 *
 * @details
 *          RPP given user args b->a[*]: (all in OSH_COORD_UNIVERSE)
 *          0,1,2: p - center point of ellipsoid
 *          3,4,5: r - vector1 from p
 *          6,7,8: s - vector2 from p
 *          9,10,11: t - vector3 from p
 *          r, s and t must be orthorgonal.
 *
 * @param[in,out] b - the body which will be setup
 *
 * @returns 1
 *
 * @author Niels Bassler
 */
static int _setup_ell(struct body *b) {

    int const nsurfs = 1;
    struct surface *sf; /* temporary surface */

    double *p, *r, *s, *t;

    p = &(b->a[0]); /* center point of ellipsoid */
    r = &(b->a[3]); /* vector1 from p */
    s = &(b->a[6]); /* vector2 from p */
    t = &(b->a[9]); /* vector3 from p */

    /* ----------- Setup translation matrix */
    b->coord = OSH_COORD_BZALIGN;
    osh_vect_setup_tmatrix_bzalign(p, r, b->t);

    /* ----------- Setup surfaces */
    osh_gemca2_add_surfaces(b, nsurfs);
    sf = b->surfs[0];
    osh_gemca2_add_surf_pars(sf, OSH_GEMCA_SURF_ELLZ);

    /* lengths of vectors are stored, r^2 = vx^2 + vy^2 + vz^2 */
    sf->p[0] = osh_vect_len2(r);
    sf->p[1] = osh_vect_len2(s);
    sf->p[2] = osh_vect_len2(t);

    return 1;
}

/**
 * @brief Setup a YZ plane for a given X.
 *
 * @details
 *          YZP given user args b->a[*]: (in OSH_COORD_UNIVERSE)
 *          0: position of X-plane
 *          Normal vector will point in positive X-direction.
 *
 * @param[in,out] b - the body which will be setup
 *
 * @returns 1
 *
 * @author Niels Bassler
 */
static int _setup_yzp(struct body *b) {
    /* YZP given user args a[*]:
       0: position of X-plane
     */

    int const nsurfs = 1;
    struct surface *sf; /* temporary surface */

    /* ----------- Setup translation matrix */
    b->coord = OSH_COORD_UNIVERSE; /* no translation matrix needed */

    /* ----------- Setup surfaces */
    osh_gemca2_add_surfaces(b, nsurfs);

    sf = b->surfs[0];
    osh_gemca2_add_surf_pars(sf, OSH_GEMCA_SURF_PLANEX);
    sf->p[0] = 1;        /* A */
    sf->p[1] = -b->a[0]; /* D = -A * x0 */

    return 1;
}

/**
 * @brief Setup a XZ plane for a given Y.
 *
 * @details
 *          XZP given user args b->a[*]: (in OSH_COORD_UNIVERSE)
 *          0: position of Y-plane
 *          Normal vector will point in positive Y-direction.
 *
 * @param[in,out] b - the body which will be setup
 *
 * @returns 1
 *
 * @author Niels Bassler
 */
static int _setup_xzp(struct body *b) {
    /* XZP given user args a[*]:
       0: position of Y-plane
     */

    int const nsurfs = 1;
    struct surface *sf; /* temporary surface */

    /* ----------- Setup translation matrix */
    b->coord = OSH_COORD_UNIVERSE; /* no translation matrix needed */

    /* ----------- Setup surfaces */
    osh_gemca2_add_surfaces(b, nsurfs);

    sf = b->surfs[0];
    osh_gemca2_add_surf_pars(sf, OSH_GEMCA_SURF_PLANEY);
    sf->p[0] = 1;        /* B */
    sf->p[1] = -b->a[0]; /* D = -B * y0 */

    return 1;
}

/**
 * @brief Setup a XY plane for a given Z.
 *
 * @details
 *          XYP given user args b->a[*]: (in OSH_COORD_UNIVERSE)
 *          0: position of Z-plane
 *          Normal vector will point in positive Z-direction.
 *
 * @param[in,out] b - the body which will be setup
 *
 * @returns 1
 *
 * @author Niels Bassler
 */
static int _setup_xyp(struct body *b) {
    /* XYP given user args a[*]:
       0: position of Z-plane
     */

    int const nsurfs = 1;
    struct surface *sf; /* temporary surface */

    /* ----------- Setup translation matrix */
    b->coord = OSH_COORD_UNIVERSE; /* no translation matrix needed */

    /* ----------- Setup surfaces */
    osh_gemca2_add_surfaces(b, nsurfs);

    sf = b->surfs[0];
    osh_gemca2_add_surf_pars(sf, OSH_GEMCA_SURF_PLANEZ);
    sf->p[0] = 1;        /* C */
    sf->p[1] = -b->a[0]; /* D = -C * z0 */

    return 1;
}

/**
 * @brief Setup a plane at an arbitrary position and oritentation.
 *
 * @details Normal vector will point into the volume defining the body.
 *
 * @param[in,out] b - the body which will be setup
 *
 * @returns 1
 *
 * @author Niels Bassler
 */
static int _setup_pla(struct body *b) {
    /* PLA given user args a[*]:
       0,1,2: normal vector of plane, pointing *away* from the body volume (this is FLUKA notation)
       3,4,5: a point in the plane
     */

    int const nsurfs = 1;
    struct surface *sf; /* temporary surface */

    double *p, *r;

    p = &(b->a[0]);
    r = &(b->a[3]);

    /* ----------- Setup translation matrix */
    b->coord = OSH_COORD_UNIVERSE; /* no translation matrix needed */

    /* ----------- Setup surfaces */
    osh_gemca2_add_surfaces(b, nsurfs);

    sf = b->surfs[0];
    osh_gemca2_add_surf_pars(sf, OSH_GEMCA_SURF_PLANE);
    osh_vect_eqpln(p, r, sf->p);

    return 1;
}

/**
 * @brief Rotate a body.
 *
 * @details
 *
 * @param[in,out] b - the body which will be rotated
 *
 * @returns 1
 *
 * @author Niels Bassler
 */
static int _setup_rot(struct body *b) {
    /* ROT given user args a[*]:
       0: body number affected
       1,2,3: rotation over some axes
     */
    (void) *b;

    // TODO
    return 1;
}

/**
 * @brief Copy a body.
 *
 * @details Copy a body to a new position.
 *
 * @param[in,out] b - the body which will be copied
 *
 * @returns 1
 *
 * @author Niels Bassler
 */
static int _setup_cpy(struct body *b) {
    /* REC given user args a[*]:
       0: body number to be copied
       1:---23: any parameter affectes ?
     */

    (void) *b;
    // TODO
    return 1;
}

/**
 * @brief Move a body.
 *
 * @details
 *
 * @param[in,out] b - the body which will be moved
 *
 * @returns 1
 *
 * @author Niels Bassler
 */
static int _setup_mov(struct body *b) {
    /* REC given user args a[*]:
       0: body number affected
       1:---23: any parameter affectes
     */
    (void) *b;
    // TODO
    return 1;
}

/**
 * @brief returns the first three vertices for ARB body in FLUKA format.
 *
 * @details The 4th vertex is redundant and not needed.
 *
 * @param[in] d - a faces vertex-index descriptor as double, see ARB manual in FLUKA. e.g. 1243.0 for vertex 1,2,4
 * and 3.
 * @param[out] i - indices -1 if error
 *
 * @returns
 *
 * @author Niels Bassler
 */
static void _vertex_index_arb_fluka(double d, int *i) {
    char s[12];
    int j, k;

    snprintf(s, 11, "%.0f", d);

    for (j = 0; j < 3; j++) {
        k = s[j] - '0';
        if ((k > 0) && (k < 9)) { /* than may max be 8 indices */
            i[j] = k - 1;         /* convert FORTRAN indices to C indices by subtracting one */
        } else {
            i[j] = -1;
        }
    }
}
