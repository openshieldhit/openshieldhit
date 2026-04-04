#ifndef _OSH_TRANSPORT
#define _OSH_TRANSPORT

#include "common/osh_coord.h"

/* position is a point with additional data : direction +  meta data on the medium where it is */
/* TODO: this is actually a ray */
struct position {
    double p[4]; /* x,y,z,E; E is total kinetic energy [MeV] (not per nucleon or per amu).*/
    double v[3]; /* unit vector pointing where particle is traveling (like CX,CY,CZ in gdatap)*/
    double rho;  /* CT-corrected density at this point [g/cm3] */
    int medium;  /* medium ID at this point, -1 if unknown */
    int zone;    /* zone number at this point, -1 if unknown */
    int system;  /* optional marker for saying what coordinate system we are in. 0 = unknown, 1 = universe ... */
};

struct step {
    double p[4]; /* p_i : start i = x,y,z,E; E is total kinetic energy [MeV] (not per nucleon or amu).*/
    double q[4]; /* q_i : stop  i = x,y,z,E; E is total kinetic energy [MeV] (not per nucleon or amu).*/
    double v[3]; /* unit vector pointing where particle is traveling */
    double ds;   /* track length of this step [cm] */
    // double ds; /* track length of this step including any detours on its way from p to q [cm] */
    double de; /* energy loss of this step (calculated positive. Energy gain is calculated negative) [MeV] */

    double rho; /* CT-corrected density at this point [g/cm3] */
    int medium; /* medium ID at this point, -1 if unknown */
    int zone;   /* zone number at this point, -1 if unknown */

    int system; /* optional marker for saying what coordinate system we are in. 0 = unknown, 1 = universe ... */
    /* TODO: consider refactoring step to hold two rays (in/out) instead of
     * p/q/v scalars. That would preserve both the incoming and outgoing
     * direction at boundaries while keeping transport metadata (ds, de, zone,
     * medium) in one record. */
};

void print_step(struct step st);
int copy_step(struct step *dest, struct step *src);
void osh_transport_move_ray(struct ray *r, double d);
void osh_transport_print_ray(struct ray const *r);

void print_pos(struct position pos);
int copy_pos(struct position *dest, struct position *src);

void print_ray_c(struct ray_c r);
void osh_clear_ray_c(struct ray_c *r);

#endif /* !_OSH_TRANSPORT */
