#include "transport/osh_transport.h"

#include <math.h>
#include <stdio.h>

#include "common/osh_const.h"
#include "common/osh_coord.h"

/**
 * @brief Print all data of a struct step
 *
 * @param[in] st - a struct step to be printed to stdout.
 *
 * @author Niels Bassler
 */
void print_step(struct step st) {
    printf(" x: %.9f             %.9f [cm]\n", st.p[0], st.q[0]);
    printf(" y: %.9f   ----->    %.9f [cm]\n", st.p[1], st.q[1]);
    printf(" z: %.9f             %.9f [cm]\n", st.p[2], st.q[2]);
    printf(" E: %.9f             %.9f [MeV]\n", st.p[3], st.q[3]);
    printf(" length  : %.9f [cm]\n", st.ds);
    printf(" delta E : %.9f [MeV]\n", st.de);
    printf(" rho       : %.9f [g/cm^3]\n", st.rho);
    printf(" medium    : %i \n", st.medium);
    printf(" zone      : %i \n", st.zone);
    printf(" c.system  : %i \n", st.system);
}

/**
 * @brief Deep copy a step struct.
 *
 * @param[in]  src - a struct step.
 * @param[out]  dest - a struct step copy.
 *
 * @returns 1
 *
 * @author Niels Bassler
 */
int copy_step(struct step *dest, struct step *src) {
    int i;

    for (i = 0; i < 4; i++) {
        dest->p[i] = src->p[i];
        dest->q[i] = src->q[i];
        if (i != 3) {
            dest->v[i] = src->v[i];
        }
    }
    dest->ds = src->ds;
    dest->de = src->de;
    dest->rho = src->rho;
    dest->medium = src->medium;
    dest->zone = src->zone;
    dest->system = src->system;

    return 1;
}

/**
 * @brief Print a position struct
 *
 * @param[in]  p - a struct position.
 *
 * @author Niels Bassler
 */
void print_pos(struct position p) {
    printf(" x: %.9f [cm]\n", p.p[0]);
    printf(" y: %.9f [cm]\n", p.p[1]);
    printf(" z: %.9f [cm]\n", p.p[2]);
    printf(" E: %.9f [MeV]\n", p.p[3]);
    printf(" rho       : %.9f [g/cm^3]\n", p.rho);
    printf(" medium    : %i \n", p.medium);
    printf(" zone      : %i \n", p.zone);
    printf(" c.system  : %i \n", p.system);
}

/**
 * @brief Deep copy a position struct.
 *
 * @param[in]  src - a struct position.
 * @param[out]  dest - a struct position copy.
 *
 * @returns 1
 *
 * @author Niels Bassler
 */
int copy_pos(struct position *dest, struct position *src) {
    int i;

    for (i = 0; i < 4; i++) {
        dest->p[i] = src->p[i];
        if (i != 3) {
            dest->v[i] = src->v[i];
        }
    }
    dest->rho = src->rho;
    dest->medium = src->medium;
    dest->zone = src->zone;
    dest->system = src->system;

    return 1;
}

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
