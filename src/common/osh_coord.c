#include "common/osh_coord.h"

#include <math.h>
#include <stdio.h>

/* sphetical coordinate angles are always ISO 80000-2:2019 convention */

/**
 * @brief for given spherical coordinate cosines calculate a unit vector pointing same direction
 *
 * @param[in] c - <cos(theta), sin(phi), cos(phi)>
 * @param[out] v - a unit vector pointing along given direction
 *
 * @returns 1
 *
 * @author Niels Bassler
 */
int osh_coord_c2v(double const *c, double *v) {

    double cost, sint, sinf, cosf;

    cost = c[0];
    sint = sqrt(1.0 - cost * cost);
    sinf = c[1];
    cosf = c[2];

    v[0] = sint * cosf;
    v[1] = sint * sinf;
    v[2] = cost;

    return 1;
}

/* TODO: check with similar function in osh_rot.c */
/**
 * @brief for a given vector calculate spherical coordinate cosines pointing same direction
 *
 * @param[in] v - a vector with length larger than 0.0
 * @param[out] c - direction in spherical coordinate cosines, <cos(theta), sin(phi), cos(phi)>
 *
 * @returns 1
 *
 * @author Niels Bassler
 */
int osh_coord_v2c(double const *v, double *c) {

    int i;
    double r = 0;
    double phi;

    for (i = 0; i < 3; i++) {
        r += v[i] * v[i];
    }

    if (r > 0.0) {

        phi = atan2(v[1], v[0]);
        c[0] = v[2] / r; /* cos(theta) */
        c[1] = sin(phi);
        c[2] = cos(phi);
    }
    return 1;
}

/**
 * @brief get spherical coordinate angles for a position vector
 *
 * @param[in] v - a vector
 * @param[out] theta - polar angle to z-axis
 * @param[out] phi - azimuthal angle φ (phi) (angle of rotation from the initial meridian plane)
 *
 * @returns 1 if |v| > 0.0, 0 otherwise.
 *
 * @author Niels Bassler
 */
int osh_coord_point2sph(double const *v, double *theta, double *phi) {
    int i;
    double r = 0;

    for (i = 0; i < 3; i++) {
        r += v[i] * v[i];
    }

    if (r > 0.0) {
        r = sqrt(r);
        *theta = acos(v[2] / r);
        *phi = atan2(v[1], v[0]);
        return 1;
    } else {
        *phi = 0.0;
        *theta = 0.0;
    }
    return 0;
}

/*
    C-specific: in a double loop, col must increase faster than row, to use caching.
    t[row][col] ...
    int a[3][4] = {
       {0, 1, 2, 3} ,          //  initializers for row indexed by 0
       {4, 5, 6, 7},           //  initializers for row indexed by 1
       {8, 9, 10, 11}          //  initializers for row indexed by 2
 */

/*
   transformation matrix:
   a0, b0, c0, d0
   a1, b1, c1, d1
   a2, b2, c2, d2
   a3, b3, c3, d3

   t[4][4] = {{a0, b0, c0, d0},  {a1, b1, c1, d1}, ... }
 */

/**
 * @brief transposes a point p[3] to point pt[3] using 4x4 affine transformation matrix t.
 *
 * @param[in]  *p - x,y,z in original coordinate system
 * @param[out] *pt - x,y,z, in new coordinate system
 * @param[in] t - 4x4 affine transformation matrix to transpose p to pt.
 *
 * @returns 1
 *
 * @see https://www.uio.no/studier/emner/matnat/ifi/nedlagte-emner/INF3320/h03/undervisningsmateriale/lecture3.pdf
 *
 * @author Niels Bassler
 */
int osh_coord_trans_point(double const p[3], double pt[3], double const t[16]) {
    int i;
    int j;
    /* this expands the given coodinates to homogeonous coords, and applies affine transformation */
    for (i = 0; i < 3; i++) {
        j = i * 4;
        pt[i] = p[0] * t[j] + p[1] * t[j + 1] + p[2] * t[j + 2] - t[j + 3];
        // TODO:
        // subtracting the translation t[j+3] is non-strandard. It should be added, an instead the
        // cooresponding t values should be stored negative.
        // Keeping this for now, in order to be compatible to gemca.f
    }
    return 1;
}

/**
 * @brief transposes a homogeonous coordinate point p[4] to h.c. point pt[4] using
 * 4x4 affine transformation matrix t.
 *
 * @details this function is prepared for future use for possible shear transformations and similar.
 *
 * @param[in]  *p - x,y,z,i in original coordinate system, i = -1 for point, = 0 for a vector
 * @param[out] *pt - x,y,z,i in new coordinate system, i = -1 for point, = 0 for a vector
 * @param[in] t - 4x4 affine transformation matrix to transpose p to pt.
 *
 * @returns 1
 *
 * @see https://www.uio.no/studier/emner/matnat/ifi/nedlagte-emner/INF3320/h03/undervisningsmateriale/lecture3.pdf
 *
 * @author Niels Bassler
 */
int osh_coord_trans_point_hc(double const p[4], double pt[4], double const t[16]) {

    int i;
    int j;

    for (i = 0; i < 4; i++) {
        j = i * 4;
        pt[i] = p[0] * t[j] + p[1] * t[j + 1] + p[2] * t[j + 2] + p[3] * t[j + 3];
    }
    return 1;
}

/**
 * @brief transposes a postion p to position pt using 4x4 affine transformation matrix t.
 *
 * @param[in]  *p - struct position input
 * @param[out] *pt - struct position output
 * @param[in] t - 4x4 affine transformation matrix to transpose p to pt.
 *
 * @returns 1
 *
 * @author Niels Bassler
 */
int osh_coord_trans_pos(struct position const *p, struct position *pt, double const t[16]) {
    int i;
    int j;

    for (i = 0; i < 3; i++) {
        j = i * 4;
        pt->p[i] = p->p[0] * t[j] + p->p[1] * t[j + 1] + p->p[2] * t[j + 2] - t[j + 3]; /* position */
        pt->v[i] = p->v[0] * t[j] + p->v[1] * t[j + 1] + p->v[2] * t[j + 2];            /* vector, often not needed */
    }
    return 1;
}

/**
 * @brief transforms a ray (position,direction vector) r to ray rt using 4x4 affine transformation matrix t.
 *
 * @param[in]  *r - input ray holding position and direction vector
 * @param[out] *rt - output ray holding position and direction vector
 * @param[in] t - 4x4 affine transformation matrix to transform r to rt.
 *
 * @returns 1
 *
 * @author Niels Bassler
 */
int osh_coord_trans_ray(struct ray_v const *r, struct ray_v *rt, double const t[16]) {
    int i;
    int j;

    for (i = 0; i < 3; i++) {
        j = i * 4;
        rt->p[i] = r->p[0] * t[j] + r->p[1] * t[j + 1] + r->p[2] * t[j + 2] - t[j + 3]; /* position */
        rt->v[i] = r->v[0] * t[j] + r->v[1] * t[j + 1] + r->v[2] * t[j + 2];            /* vector, often not needed */
    }
    return 1;
}

// as before, but with struct ray instead of struct ray_v. Temporary solution, until struct ray_v has been renamed to
// struct ray.
int osh_coord_trans_ray_r(struct ray const *r, struct ray *rt, double const t[16]) {
    int i;
    int j;

    for (i = 0; i < 3; i++) {
        j = i * 4;
        rt->p[i] = r->p[0] * t[j] + r->p[1] * t[j + 1] + r->p[2] * t[j + 2] - t[j + 3]; /* position */
        rt->cp[i] = r->cp[0] * t[j] + r->cp[1] * t[j + 1] + r->cp[2] * t[j + 2];        /* vector, often not needed */
    }
    return 1;
}

// TODO:
// int osh_coord_transinv_step(struct step st, struct step stt, double const t[4][4]) {
//     return 1;
// }

/**
 * @brief Inverts a 4x4 matrix stored in serially C-style.
 *
 * @param[in]  *m . 4x4 double precision matrix entries
 * @param[out] *im . 4x4 double precision inverse matrix
 *
 * @returns 1 or 0 if determinant is 0
 *
 * @śee https://stackoverflow.com/questions/1148309/inverting-a-4x4-matrix
 *
 * @author MESA library, adaptions by Niels Bassler
 *
 */
int osh_invert_matrix(const double m[16], double im[16]) {
    double det;
    int i;

    im[0] = m[5] * m[10] * m[15] - m[5] * m[11] * m[14] - m[9] * m[6] * m[15] + m[9] * m[7] * m[14]
            + m[13] * m[6] * m[11] - m[13] * m[7] * m[10];

    im[4] = -m[4] * m[10] * m[15] + m[4] * m[11] * m[14] + m[8] * m[6] * m[15] - m[8] * m[7] * m[14]
            - m[12] * m[6] * m[11] + m[12] * m[7] * m[10];

    im[8] = m[4] * m[9] * m[15] - m[4] * m[11] * m[13] - m[8] * m[5] * m[15] + m[8] * m[7] * m[13]
            + m[12] * m[5] * m[11] - m[12] * m[7] * m[9];

    im[12] = -m[4] * m[9] * m[14] + m[4] * m[10] * m[13] + m[8] * m[5] * m[14] - m[8] * m[6] * m[13]
             - m[12] * m[5] * m[10] + m[12] * m[6] * m[9];

    im[1] = -m[1] * m[10] * m[15] + m[1] * m[11] * m[14] + m[9] * m[2] * m[15] - m[9] * m[3] * m[14]
            - m[13] * m[2] * m[11] + m[13] * m[3] * m[10];

    im[5] = m[0] * m[10] * m[15] - m[0] * m[11] * m[14] - m[8] * m[2] * m[15] + m[8] * m[3] * m[14]
            + m[12] * m[2] * m[11] - m[12] * m[3] * m[10];

    im[9] = -m[0] * m[9] * m[15] + m[0] * m[11] * m[13] + m[8] * m[1] * m[15] - m[8] * m[3] * m[13]
            - m[12] * m[1] * m[11] + m[12] * m[3] * m[9];

    im[13] = m[0] * m[9] * m[14] - m[0] * m[10] * m[13] - m[8] * m[1] * m[14] + m[8] * m[2] * m[13]
             + m[12] * m[1] * m[10] - m[12] * m[2] * m[9];

    im[2] = m[1] * m[6] * m[15] - m[1] * m[7] * m[14] - m[5] * m[2] * m[15] + m[5] * m[3] * m[14] + m[13] * m[2] * m[7]
            - m[13] * m[3] * m[6];

    im[6] = -m[0] * m[6] * m[15] + m[0] * m[7] * m[14] + m[4] * m[2] * m[15] - m[4] * m[3] * m[14] - m[12] * m[2] * m[7]
            + m[12] * m[3] * m[6];

    im[10] = m[0] * m[5] * m[15] - m[0] * m[7] * m[13] - m[4] * m[1] * m[15] + m[4] * m[3] * m[13] + m[12] * m[1] * m[7]
             - m[12] * m[3] * m[5];

    im[14] = -m[0] * m[5] * m[14] + m[0] * m[6] * m[13] + m[4] * m[1] * m[14] - m[4] * m[2] * m[13]
             - m[12] * m[1] * m[6] + m[12] * m[2] * m[5];

    im[3] = -m[1] * m[6] * m[11] + m[1] * m[7] * m[10] + m[5] * m[2] * m[11] - m[5] * m[3] * m[10] - m[9] * m[2] * m[7]
            + m[9] * m[3] * m[6];

    im[7] = m[0] * m[6] * m[11] - m[0] * m[7] * m[10] - m[4] * m[2] * m[11] + m[4] * m[3] * m[10] + m[8] * m[2] * m[7]
            - m[8] * m[3] * m[6];

    im[11] = -m[0] * m[5] * m[11] + m[0] * m[7] * m[9] + m[4] * m[1] * m[11] - m[4] * m[3] * m[9] - m[8] * m[1] * m[7]
             + m[8] * m[3] * m[5];

    im[15] = m[0] * m[5] * m[10] - m[0] * m[6] * m[9] - m[4] * m[1] * m[10] + m[4] * m[2] * m[9] + m[8] * m[1] * m[6]
             - m[8] * m[2] * m[5];

    det = m[0] * im[0] + m[1] * im[4] + m[2] * im[8] + m[3] * im[12];

    if (det == 0)
        return 0;

    det = 1.0 / det;

    for (i = 0; i < 16; i++)
        im[i] = im[i] * det;

    return 1;
}
