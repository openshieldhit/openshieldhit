#include "osh_vect.h"
#include "osh_logger.h"
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

    int i;
    double l = 0;

    for (i = 0; i < OSH_VECT_DIM; i++) {
        l += u[i] * u[i];
    }
    return l;
}

double osh_vect_dot(double const *u, double const *v) {
    int i;
    double dot = 0;
    for (i = 0; i < OSH_VECT_DIM; i++) {
        dot += u[i] * v[i];
    }
    return dot;
}

void osh_vect_cross(double const *u, double const *v, double *w) {
    w[0] = u[1] * v[2] - u[2] * v[1];
    w[1] = u[2] * v[0] - u[0] * v[2];
    w[2] = u[0] * v[1] - u[1] * v[0];
}

double osh_vect_sproj(double const *u, double const *v) {

    double w[OSH_VECT_DIM];
    osh_vect_norm2(v, w); /* w is the unit vector of v */
    return osh_vect_dot(u, w);
}

void osh_vect_norm(double *u) {
    double d = 0;
    int i;

    for (i = 0; i < OSH_VECT_DIM; i++) {
        d += u[i] * u[i];
    }
    if (d > 0.0) {
        d = 1.0 / sqrt(d);
    } else {
        osh_error(EX_SOFTWARE, "osh_vect_norm() division by zero.\n");
    }

    for (i = 0; i < OSH_VECT_DIM; i++) {
        u[i] *= d;
    }
    return;
}

void osh_vect_norm2(double const *u, double *v) {
    double d = 0;
    int i;

    for (i = 0; i < OSH_VECT_DIM; i++) {
        d += u[i] * u[i];
    }

    if (d > 0.0) {
        d = 1.0 / sqrt(d);
    } else {
        osh_error(EX_SOFTWARE, "osh_vect_norm2() division by zero.\n");
    }

    for (i = 0; i < OSH_VECT_DIM; i++) {
        v[i] = u[i] * d;
    }

    return;
}

void osh_vect_orthogonal_basis(double const *w, double *u, double *v) {

    int i;
    /* define unit vectors e1, e2 and e3 forming standard basis */
    double const im[3][3] = {
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1},
    };
    double sign = 1.;
    /*
       We need to maintain certain orientation of the basis, namely
       u x v = w
       v x w = u
       w x u = v
     */

    /*
       First lets check if w is colinear with w_1 (if it lies along X axis)
       u = w x e3
     */
    osh_vect_cross(w, im[2], u);

    /*
       If |u| == 0 then u = 0 which means that w is colinear with e3.
       If |u| > 0 then u and e3 span a plane in 3D space and v is well defined
       normal vector. It is safe here to compare with zero, as this is in fact
       check for non-zero vector afterwards only multiplications are used, so
       even if `u` has a small norm, the numerical error won't propagate.
     */
    if (osh_vect_len2(u) > 0) {
        /* v = w x u */
        osh_vect_cross(w, u, v);
    } else {
        if (u[2] < 0) {
            sign = -1.;
        }
        /*
           If u is oriented in the same direction as e1 (u==e1, as both are unit
           vectors) then we create a standard basis (u=e1, v=e2, w=e3). If u is
           oriented in opposite direction as e1(u == -e1, as both are unit
           vectors) then we choose another basis which has the same orientation
           as standard basis (u=-e1, v=-e2, w=e3).
         */
        for (i = 0; i < OSH_VECT_DIM; i++) {
            u[i] = sign * im[0][i];
            v[i] = im[1][i];
        }
    }
}

void osh_vect_eqpln(double const *p, double const *u, double *pp) {
    int i;

    for (i = 0; i < OSH_VECT_DIM; i++)
        pp[i] = u[i];

    pp[3] = -osh_vect_dot(u, p);

    return;
}

void osh_vect_rot_y(double alpha, double *u) {
    double tx, tz;

    tx = u[0];
    tz = u[2];

    u[0] = cos(alpha) * tx - sin(alpha) * tz;
    u[2] = sin(alpha) * tx + cos(alpha) * tz;

    return;
}

void osh_vect_rot_z(double alpha, double *u) {
    double tx, ty;

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
        printf("% .3f ", tm[i]); // fixed width printout
        if (!((i + 1) % 4))
            printf("\n    ");
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
     *  M e1 = M * [1,0,0]^T = [Sx*1 + Tx*0 + Rx*0,Sy*1 + Ty*0 + Ry*0,Sz*1 +
     Tz*0 + Rz*0]=[Sx,Sy,Sz] = S
     *  M e2 = M * [0,1,0]^T = [Sx*0 + Tx*1 + Rx*0,Sy*0 + Ty*1 + Ry*0,Sz*0 +
     Tz*1 + Rz*0]=[Tx,Ty,Tz] = T
     *  M e3 = M * [0,0,1]^T = [Sx*0 + Tx*0 + Rx*1,Sy*0 + Ty*0 + Ry*1,Sz*0 +
     Tz*0 + Rz*1]=[Rx,Ry,Rz] = R
     */

    /**
     * Affine transformation matrix has following indices
     *     C-specific: in a double loop, col must increase faster than row, to
     use caching.
     *   t[row][col] ...
     *   int a[3][4] = {
     *       {0, 1, 2, 3} ,          initializers for row indexed by 0
     *       {4, 5, 6, 7},           initializers for row indexed by 1
     *       {8, 9, 10, 11}          initializers for row indexed by 2
     */

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
