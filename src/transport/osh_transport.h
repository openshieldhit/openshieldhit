#ifndef OSH_TRANSPORT_H
#define OSH_TRANSPORT_H

#include "common/osh_coord.h"
#include "common/osh_step.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * struct step and struct position are defined in common/osh_step.h so that
 * domain modules (scoring, gemca) can receive them without depending on
 * transport/.  This header re-exports them via the include above.
 */

void osh_transport_move_ray(struct ray *r, double d);
void osh_transport_print_ray(struct ray const *r);

void print_ray_c(struct ray_c r);
void osh_clear_ray_c(struct ray_c *r);

#ifdef __cplusplus
}
#endif

#endif /* OSH_TRANSPORT_H */
