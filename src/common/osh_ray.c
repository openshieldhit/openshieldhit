#include "common/osh_ray.h"

#include <math.h>
#include <stdio.h>

#include "common/osh_const.h"
#include "common/osh_ray_hd.h"

void osh_ray_move(struct ray *r, double d) {
    int i;

    for (i = 0; i < 3; i++) {
        r->p[i] += r->cp[i] * d;
    }
}

void osh_ray_print(struct ray const *r) {
    printf(" x: %.9f [cm]\n", r->p[0]);
    printf(" y: %.9f [cm]\n", r->p[1]);
    printf(" z: %.9f [cm]\n", r->p[2]);
    printf(" cx: %.9f\n", r->cp[0]);
    printf(" cy: %.9f\n", r->cp[1]);
    printf(" cz: %.9f\n", r->cp[2]);
    printf(" c.system: %i\n", r->system);
}

void osh_ray_c_print(struct ray_c r) {
    printf(" x,y,z      : % .4f % .4f % .4f [cm]\n", r.p[0], r.p[1], r.p[2]);
    printf(" ct,sf,cf   : % .4f % .4f % .4f\n", r.c[0], r.c[1], r.c[2]);
    printf(
        " theta, phi : %.4f %.4f deg\n", acos(r.c[0]) * OSH_M_1_PI * 180.0, atan2(r.c[1], r.c[2]) * OSH_M_1_PI * 180.0);
    printf(" c.system   : %i \n", r.system);
}

void osh_ray_c_clear(struct ray_c *r) {
    int i;

    for (i = 0; i < 4; i++) {
        r->p[i] = 0.0; /* p[3] = energy [MeV]; zero-initialise to avoid leaking garbage */
    }
    for (i = 0; i < 3; i++) {
        r->c[i] = 1.0;
    }
    r->c[1] = 0.0; /* travel along +Z: c = (cos θ=1, sin φ=0, cos φ=1) */
    r->system = OSH_COORD_PZALIGN;
}

int osh_ray_v_transform(struct ray_v const *r, struct ray_v *rt, double const t[16]) {
    int i;
    int j;

    for (i = 0; i < 3; i++) {
        j = i * 4;
        /* Translation is subtracted — GEMCA/SHIELD-HIT sign convention. */
        rt->p[i] = r->p[0] * t[j] + r->p[1] * t[j + 1] + r->p[2] * t[j + 2] - t[j + 3];
        rt->v[i] = r->v[0] * t[j] + r->v[1] * t[j + 1] + r->v[2] * t[j + 2];
    }
    rt->p[3] = r->p[3];     /* energy is invariant under coordinate transform */
    rt->system = r->system; /* caller should update to the target system after the call */
    return 1;
}

/* Temporary — remove once struct ray is retired in favour of struct ray_v. */
/* Body lives in osh_ray_hd.h so device kernels compile the same lines. */
int osh_ray_transform(struct ray const *r, struct ray *rt, double const t[16]) {
    return _osh_ray_transform_hd(r, rt, t);
}
