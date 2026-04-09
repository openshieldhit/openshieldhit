#include "transport/osh_transport.h"

#include <math.h>
#include <stdio.h>

#include "common/osh_const.h"
#include "common/osh_coord.h"

/**
 * @brief Move a ray along its path by a distance d.
 *
 * @param[in,out] r - a ray to be moved.
 * @param[in] d - the distance to move the ray.
 *
 * @author Niels Bassler
 */
void osh_transport_move_ray(struct ray *r, double d) {
    int i;

    for (i = 0; i < 3; i++) {
        r->p[i] += r->cp[i] * d;
    }
}

/**
 * @brief Print a struct ray
 *
 * @param[in]  ray_c - a ray as used in gemca
 *
 * @author Niels Bassler
 *
 */
void osh_transport_print_ray(struct ray const *r) {
    printf(" x: %.9f [cm]\n", r->p[0]);
    printf(" y: %.9f [cm]\n", r->p[1]);
    printf(" z: %.9f [cm]\n", r->p[2]);
    printf(" cx: %.9f\n", r->cp[0]);
    printf(" cy: %.9f\n", r->cp[1]);
    printf(" cz: %.9f\n", r->cp[2]);
    printf(" c.system: %i\n", r->system);
}

/**
 * @brief Print a ray_c
 *
 * @param[in]  ray_c - a ray in spherical coordinates
 *
 * @author Niels Bassler
 */
void print_ray_c(struct ray_c r) {
    printf(" x,y,z      : % .4f % .4f % .4f [cm]\n", r.p[0], r.p[1], r.p[2]);
    printf(" ct,sf,cf   : % .4f % .4f % .4f\n", r.c[0], r.c[1], r.c[2]);
    printf(
        " theta, phi : %.4f %.4f deg\n", acos(r.c[0]) * OSH_M_1_PI * 180.0, atan2(r.c[1], r.c[2]) * OSH_M_1_PI * 180.0);
    printf(" c.system   : %i \n", r.system);
}

/**
 * @brief Initialize a ray for travel along the Z-axis.
 *
 * @param[in]  ray_c - a ray in spherical coordinates
 *
 * @author Niels Bassler
 */
void osh_clear_ray_c(struct ray_c *r) {
    int i;

    for (i = 0; i < 3; i++) {
        r->p[i] = 0;
        r->c[i] = 1;
    }
    r->c[1] = 0; /* initialized travel along z: ray->c = (1,0,1) */
    r->system = OSH_COORD_PZALIGN;
}
