#include "osh_vect.h"
#include "osh_vect_hd.h"

#include <math.h>
#include <stdio.h>

void osh_vect_add(double const *p, double const *q, double *u) {
    int i;
    for (i = 0; i < OSH_VECT_DIM; i++) {
        u[i] = p[i] + q[i];
    }
}

void osh_vect_addmul(double const *p, double const *q, double d, double *u) {
    int i;
    for (i = 0; i < OSH_VECT_DIM; i++) {
        u[i] = p[i] + q[i] * d;
    }
}

void osh_vect_sub(double const *p, double const *q, double *u) {
    int i;
    for (i = 0; i < OSH_VECT_DIM; i++) {
        u[i] = p[i] - q[i];
    }
}

void osh_vect_copy(double const *u, double *v) {
    int i;
    for (i = 0; i < OSH_VECT_DIM; i++) {
        v[i] = u[i];
    }
}

void osh_vect_reverse(double const *u, double *v) {
    int i;
    for (i = 0; i < OSH_VECT_DIM; i++) {
        v[i] = -u[i];
    }
}

double osh_vect_len2(double const *u) {
    return _osh_vect_len2_hd(u);
}

void osh_vect_cross(double const *u, double const *v, double *w) {
    _osh_vect_cross_hd(u, v, w);
}

double osh_vect_sproj(double const *u, double const *v) {

    double w[OSH_VECT_DIM];
    osh_vect_norm2(v, w); /* w is the unit vector of v */
    return osh_vect_dot(u, w);
}

void osh_vect_norm(double *u) {
    _osh_vect_norm_hd(u);
}

void osh_vect_norm2(double const *u, double *v) {
    _osh_vect_norm2_hd(u, v);
}

/* TODO: new code should use osh_vect_orthogonal_basis_norm instead. */
void osh_vect_orthogonal_basis(double const *w, double *u, double *v) {
    _osh_vect_orthogonal_basis_hd(w, u, v);
}

void osh_vect_orthogonal_basis_norm(double const *w, double *u, double *v) {
    _osh_vect_orthogonal_basis_norm_hd(w, u, v);
}

void osh_vect_eqpln(double const *p, double const *u, double *pp) {
    int i;

    for (i = 0; i < OSH_VECT_DIM; i++) {
        pp[i] = u[i];
    }

    pp[3] = -osh_vect_dot(u, p);

    return;
}

void osh_vect_rot_y(double alpha, double *u) {
    double tx;
    double tz;

    tx = u[0];
    tz = u[2];

    u[0] = cos(alpha) * tx - sin(alpha) * tz;
    u[2] = sin(alpha) * tx + cos(alpha) * tz;

    return;
}

void osh_vect_rot_z(double alpha, double *u) {
    double tx;
    double ty;

    tx = u[0];
    ty = u[1];

    u[0] = cos(alpha) * tx + sin(alpha) * ty;
    u[1] = -sin(alpha) * tx + cos(alpha) * ty;

    return;
}

void osh_vect_print(double const *v) {
    int i;
    for (i = 0; i < OSH_VECT_DIM; i++) {
        printf("%.3f\n", v[i]);
    }
}

void osh_vect_matrix4_print(double const *tm) {
    int i;

    printf("Transformation matrix:\n");
    printf("    ");
    for (i = 0; i < 16; i++) {
        printf("% .3f ", tm[i]); /* fixed width printout */
        if (!((i + 1) % 4)) {
            printf("\n    ");
        }
    }
    printf("\n");
}

void osh_vect_setup_tmatrix_bzalign(double *p, double *r, double *tm) {

    double s[OSH_VECT_DIM];
    double t[OSH_VECT_DIM];
    double r_norm[OSH_VECT_DIM];

    /* find two vectors S and T such that S,T,R forms orthogonal set
       after normalisation S/|S| T/|T| R/|R| form right-handed basis */

    /* calculate unit vector r_norm = R/|R| */
    osh_vect_norm2(r, r_norm);

    osh_vect_orthogonal_basis(r_norm, s, t);

    /* lets ensure S and T are normalised */
    osh_vect_norm(s);
    osh_vect_norm(t);

    /* osh_vect_print(r); printf("\n"); */
    /* osh_vect_print(s); printf("\n"); */
    /* osh_vect_print(t); printf("\n"); */

    /* RST will be positive oriented basis */
    /* here we need to use R as last basis vector, therefore to maintain correct
     * orientation */
    /* we choose STR order */

    /**
     * Lets then build a 3x3 rotation matrix M which maps
     *    e1  --->  S
     *    e2  --->  T
     *    e3  --->  R
     * such matrix will have columns with S,T,R vectors respectively
     *       Sx  Tx  Rx
     *   M = Sy  Ty  Ry
     *       Sz  Tz  Rz
     * quick crosscheck how M works
     *
     *  M e1 = M * [1,0,0]^T = [Sx*1 + Tx*0 + Rx*0,Sy*1 + Ty*0 + Ry*0,Sz*1 + Tz*0 + Rz*0]=[Sx,Sy,Sz] = S
     *  M e2 = M * [0,1,0]^T = [Sx*0 + Tx*1 + Rx*0,Sy*0 + Ty*1 + Ry*0,Sz*0 + Tz*1 + Rz*0]=[Tx,Ty,Tz] = T
     *  M e3 = M * [0,0,1]^T = [Sx*0 + Tx*0 + Rx*1,Sy*0 + Ty*0 + Ry*1,Sz*0 + Tz*0 + Rz*1]=[Rx,Ry,Rz] = R
     */

    /**
     * Affine transformation matrix has following indices
     *     C-specific: in a double loop, col must increase faster than row, to
     * use caching.
     *   t[row][col] ...
     *   int a[3][4] = {
     *       {0, 1, 2, 3} ,          initializers for row indexed by 0
     *       {4, 5, 6, 7},           initializers for row indexed by 1
     *       {8, 9, 10, 11}          initializers for row indexed by 2
     * */
    /* First row: x coordindates and translation vector */
    tm[0] = s[0];
    tm[1] = t[0];
    tm[2] = r_norm[0];
    tm[3] = osh_vect_dot(p, s); /* length of projection of P on S, <S,P> */

    /* Second row: Y coordindates and translation vector */
    tm[4] = s[1];
    tm[5] = t[1];
    tm[6] = r_norm[1];
    tm[7] = osh_vect_dot(p, t); /* length of projection of P on T, <T,P> */

    /* Third row: Y coordindates and translation vector */
    tm[8] = s[2];
    tm[9] = t[2];
    tm[10] = r_norm[2];
    tm[11] = osh_vect_dot(p, r_norm); /* length of projection of P on R, <R/|R|,P> */

    /* Last row */
    tm[12] = 0;
    tm[13] = 0;
    tm[14] = 0;
    tm[15] = 1;
}

void osh_vect_setup_tmatrix_bzalign_affine(double const *p_local, double const *r_world, double *tm) {
    double s[OSH_VECT_DIM];
    double t[OSH_VECT_DIM];
    double r_norm[OSH_VECT_DIM];
    double p_world[OSH_VECT_DIM];
    int i;

    osh_vect_norm2(r_world, r_norm);
    osh_vect_orthogonal_basis(r_norm, s, t);
    osh_vect_norm(s);
    osh_vect_norm(t);

    tm[0] = s[0];
    tm[1] = t[0];
    tm[2] = r_norm[0];

    tm[4] = s[1];
    tm[5] = t[1];
    tm[6] = r_norm[1];

    tm[8] = s[2];
    tm[9] = t[2];
    tm[10] = r_norm[2];

    for (i = 0; i < OSH_VECT_DIM; i++) {
        int j = i * 4;
        p_world[i] = tm[j] * p_local[0] + tm[j + 1] * p_local[1] + tm[j + 2] * p_local[2];
        tm[j + 3] = p_world[i];
    }

    tm[12] = 0.0;
    tm[13] = 0.0;
    tm[14] = 0.0;
    tm[15] = 1.0;
}

void osh_vect_trans_point_affine(double const *p, double *pt, double const *tm) {
    pt[0] = p[0] * tm[0] + p[1] * tm[1] + p[2] * tm[2] + tm[3];
    pt[1] = p[0] * tm[4] + p[1] * tm[5] + p[2] * tm[6] + tm[7];
    pt[2] = p[0] * tm[8] + p[1] * tm[9] + p[2] * tm[10] + tm[11];
}

void osh_vect_trans_vector_affine(double const *v, double *vt, double const *tm) {
    vt[0] = v[0] * tm[0] + v[1] * tm[1] + v[2] * tm[2];
    vt[1] = v[0] * tm[4] + v[1] * tm[5] + v[2] * tm[6];
    vt[2] = v[0] * tm[8] + v[1] * tm[9] + v[2] * tm[10];
}

double osh_vect_dot(double const *u, double const *v) {
    return _osh_vect_dot_hd(u, v);
}
