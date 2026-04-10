#ifndef OSH_COORD_H
#define OSH_COORD_H

/* list of coordinate system identifiers */
#define OSH_COORD_UNKNOWN 0  /* Unknown or not set */
#define OSH_COORD_UNIVERSE 1 /* Simulation universe, as specified by user */
#define OSH_COORD_PZALIGN 2  /* Particle Z-ALIGNed systen: Particle initially travels along z-axis. Used by STRAGL */
#define OSH_COORD_VOXELCT 3  /* CT system, where lowest X,Y,Z corner is at (0,0,0) cm and slices along Z-axis */
#define OSH_COORD_BZALIGN 4  /* System which aligned so one body corner is at (0,0,0) cm and turned along z-axis */
#define OSH_COORD_BCALIGN 5  /* System which aligned so the body center is at (0,0,0) cm */

/*
 * struct point: a location in space with optional energy and coordinate tag.
 * Ray types (struct ray, struct ray_v, struct ray_c) are in osh_ray.h.
 */
struct point {
    double p[4]; /* x,y,z [cm]; p[3] = total kinetic energy [MeV] when used as particle state */
    int system;  /* coordinate system (OSH_COORD_*); 0 = unknown */
};

struct position; /* defined in osh_step.h */

int osh_coord_c2v(double const *c, double *v);
int osh_coord_v2c(double const *v, double *c);
int osh_coord_point2sph(double const *v, double *theta, double *phi);
int osh_coord_trans_point(double const p[3], double pt[3], double const t[16]);
int osh_coord_trans_point_hc(double const p[4], double pt[4], double const t[16]);
int osh_coord_trans_pos(struct position const *p, struct position *pt, double const t[16]);

int osh_invert_matrix(double const m[16], double im[16]);

#endif /* OSH_COORD_H */
