#ifndef OSH_RAY_HD_H
#define OSH_RAY_HD_H

/*
 * osh_ray_hd.h — device-compilable ray transform.
 *
 * The body is marked OSH_HD static inline so it compiles both as a host
 * function (via plain C compilation) and as a device function (via nvcc
 * compilation with __host__ __device__).  osh_ray.c includes this header
 * and re-exports the function with its unchanged public signature.
 */

#include "common/osh_hd.h"
#include "common/osh_ray.h"

/*
 * _osh_ray_transform_hd() — body of osh_ray_transform() (see osh_ray.c).
 * Affine universe → body-local transform with a 4x4 row-major matrix; the
 * translation column is subtracted (GEMCA/SHIELD-HIT sign convention).
 */
OSH_HD static inline int _osh_ray_transform_hd(struct ray const *r, struct ray *rt, double const t[16]) {
    int i;
    int j;

    for (i = 0; i < 3; i++) {
        j = i * 4;
        rt->p[i] = r->p[0] * t[j] + r->p[1] * t[j + 1] + r->p[2] * t[j + 2] - t[j + 3];
        rt->cp[i] = r->cp[0] * t[j] + r->cp[1] * t[j + 1] + r->cp[2] * t[j + 2];
    }
    rt->system = r->system; /* caller should update to the target system after the call */
    return 1;
}

#endif /* OSH_RAY_HD_H */
