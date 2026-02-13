#ifndef _OSH_TRANSPORT
#define _OSH_TRANSPORT

struct point {
    double p[4]; /* x,y,z,E; E is total kinetic energy [MeV] (not per nucleon or per amu).*/
    int system;  /* optional marker for saying what coordinate system we are in. 0 = unknown, 1 = universe ... */
};

struct ray {
    double p[3];
    double cp[3]; /* direction vector */ // TODO: refactor to v
    int system;                          /* coordinate system */
};

/* position is a point with additional data : direction +  meta data on the medium where it is */
/* Todo: this is actually a ray */
struct position {
    double p[4]; /* x,y,z,E; E is total kinetic energy [MeV] (not per nucleon or per amu).*/
    double v[3]; /* unit vector pointing where particle is traveling (like CX,CY,CZ in gdatap)*/
    // double c[3];        /* spherical coordinate direction cosines (cos(theta), sin(phi), cos(phi) */
    double rho; /* CT-corrected density at this point [g/cm3] */
    int medium; /* medium ID at this point, -1 if unknown */
    int zone;   /* zone number at this point, -1 if unknown */
    int system; /* optional marker for saying what coordinate system we are in. 0 = unknown, 1 = universe ... */
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
};

/* a ray is a primitive struct which merily contains position, optional energy, and a direction vector */
/* a ray has no length */
struct ray_v {   // TODO: rename to ray
    double p[4]; /* x,y,z,E; E is total kinetic energy [MeV] (not per nucleon or per amu).*/
    double v[3]; /* unit vector pointing where particle is traveling (like CX,CY,CZ in gdatap)*/
    unsigned char
        system; /* optional marker for saying what coordinate system we are in. 0 = unknown, 1 = universe ... */
};

/* alternatively a ray can also be defined by cosines of spherical coordinates */
struct ray_c {
    double p[4]; /* x,y,z,E; E is total kinetic energy [MeV] (not per nucleon or per amu).*/
    double c[3]; /* spherical coordinate direction cosines (cos(theta), sin(phi), cos(phi) */
    unsigned char
        system; /* optional marker for saying what coordinate system we are in. 0 = unknown, 1 = universe ... */
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
