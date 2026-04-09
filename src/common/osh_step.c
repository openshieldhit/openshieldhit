#include "common/osh_step.h"

#include <stdio.h>

/**
 * @brief Print all fields of a struct step to stdout.
 *
 * @param[in] st  Step to print.
 */
void osh_step_print(struct step const *st) {
    printf(" x: %.9f             %.9f [cm]\n", st->p[0], st->q[0]);
    printf(" y: %.9f   ----->    %.9f [cm]\n", st->p[1], st->q[1]);
    printf(" z: %.9f             %.9f [cm]\n", st->p[2], st->q[2]);
    printf(" E: %.9f             %.9f [MeV]\n", st->p[3], st->q[3]);
    printf(" v: %.9f %.9f %.9f\n", st->v[0], st->v[1], st->v[2]);
    printf(" w: %.9f %.9f %.9f\n", st->w[0], st->w[1], st->w[2]);
    printf(" ds:     %.9f [cm]\n", st->ds);
    printf(" de:     %.9f [MeV]\n", st->de);
    printf(" rho:    %.9f [g/cm3]\n", st->rho);
    printf(" medium: %i\n", st->medium);
    printf(" zone:   %i\n", st->zone);
    printf(" system: %i\n", st->system);
}

/**
 * @brief Copy a struct step.
 *
 * @param[out] dst  Destination.
 * @param[in]  src  Source.
 */
void osh_step_copy(struct step *dst, struct step const *src) {
    *dst = *src;
}

/**
 * @brief Print all fields of a struct position to stdout.
 *
 * @param[in] pos  Position to print.
 */
void osh_position_print(struct position const *pos) {
    printf(" x: %.9f [cm]\n", pos->p[0]);
    printf(" y: %.9f [cm]\n", pos->p[1]);
    printf(" z: %.9f [cm]\n", pos->p[2]);
    printf(" E: %.9f [MeV]\n", pos->p[3]);
    printf(" v: %.9f %.9f %.9f\n", pos->v[0], pos->v[1], pos->v[2]);
    printf(" rho:    %.9f [g/cm3]\n", pos->rho);
    printf(" medium: %i\n", pos->medium);
    printf(" zone:   %i\n", pos->zone);
    printf(" system: %i\n", pos->system);
}

/**
 * @brief Copy a struct position.
 *
 * @param[out] dst  Destination.
 * @param[in]  src  Source.
 */
void osh_position_copy(struct position *dst, struct position const *src) {
    *dst = *src;
}
