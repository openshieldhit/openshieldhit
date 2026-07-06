#ifndef OSH_VECT_HD_H
#define OSH_VECT_HD_H

/*
 * osh_vect_hd.h — device-compilable vector helpers used by the GPU transport
 * kernel.
 *
 * Bodies are marked OSH_HD static inline so they compile both as host
 * functions (via plain C compilation) and as device functions (via nvcc
 * with __host__ __device__).  The original .c file includes this header
 * and re-exports each with its unchanged public signature.
 */

#include "common/osh_hd.h"
#include "osh_vect.h"

#include <math.h>

OSH_HD static inline double _osh_vect_dot_hd(double const *u, double const *v) {
    int i;
    double dot = 0;
    for (i = 0; i < OSH_VECT_DIM; i++) {
        dot += u[i] * v[i];
    }
    return dot;
}

OSH_HD static inline double _osh_vect_len2_hd(double const *u) {
    int i;
    double l = 0;

    for (i = 0; i < OSH_VECT_DIM; i++) {
        l += u[i] * u[i];
    }
    return l;
}

OSH_HD static inline void _osh_vect_cross_hd(double const *u, double const *v, double *w) {
    w[0] = u[1] * v[2] - u[2] * v[1];
    w[1] = u[2] * v[0] - u[0] * v[2];
    w[2] = u[0] * v[1] - u[1] * v[0];
}

OSH_HD static inline void _osh_vect_norm_hd(double *u) {
    double d = 0;
    int i;

    for (i = 0; i < OSH_VECT_DIM; i++) {
        d += u[i] * u[i];
    }
    if (d > 0.0) {
        d = 1.0 / sqrt(d);
    }

    for (i = 0; i < OSH_VECT_DIM; i++) {
        u[i] *= d;
    }
}

OSH_HD static inline void _osh_vect_norm2_hd(double const *u, double *v) {
    double d = 0;
    int i;

    for (i = 0; i < OSH_VECT_DIM; i++) {
        d += u[i] * u[i];
    }

    if (d > 0.0) {
        d = 1.0 / sqrt(d);
    }

    for (i = 0; i < OSH_VECT_DIM; i++) {
        v[i] = u[i] * d;
    }
}

OSH_HD static inline void _osh_vect_orthogonal_basis_hd(double const *w, double *u, double *v) {
    int i;
    double const im[3][3] = {
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1},
    };
    double sign = 1.;

    _osh_vect_cross_hd(w, im[2], u);

    if (_osh_vect_len2_hd(u) > 0) {
        _osh_vect_cross_hd(w, u, v);
    } else {
        if (u[2] < 0) {
            sign = -1.;
        }
        for (i = 0; i < OSH_VECT_DIM; i++) {
            u[i] = sign * im[0][i];
            v[i] = im[1][i];
        }
    }
}

OSH_HD static inline void _osh_vect_orthogonal_basis_norm_hd(double const *w, double *u, double *v) {
    double ax = fabs(w[0]);
    double ay = fabs(w[1]);
    double az = fabs(w[2]);
    double norm_inv;

    if (ax <= ay && ax <= az) {
        u[0] = 1.0 - w[0] * w[0];
        u[1] = 0.0 - w[1] * w[0];
        u[2] = 0.0 - w[2] * w[0];
    } else if (ay <= az) {
        u[0] = 0.0 - w[0] * w[1];
        u[1] = 1.0 - w[1] * w[1];
        u[2] = 0.0 - w[2] * w[1];
    } else {
        u[0] = 0.0 - w[0] * w[2];
        u[1] = 0.0 - w[1] * w[2];
        u[2] = 1.0 - w[2] * w[2];
    }
    norm_inv = 1.0 / sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
    u[0] *= norm_inv;
    u[1] *= norm_inv;
    u[2] *= norm_inv;

    v[0] = w[1] * u[2] - w[2] * u[1];
    v[1] = w[2] * u[0] - w[0] * u[2];
    v[2] = w[0] * u[1] - w[1] * u[0];
}

#endif /* OSH_VECT_HD_H */
