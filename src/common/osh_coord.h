#ifndef _OSH_COORD
#define _OSH_COORD

/* list of coordinate system identifiers */
#define OSH_COORD_UNKNOWN 0  /* Unknown or not set */
#define OSH_COORD_UNIVERSE 1 /* Simulation universe, as specified by user */
#define OSH_COORD_PZALIGN 2  /* Particle Z-ALIGNed systen: Particle initially travels along z-axis. Used by STRAGL */
#define OSH_COORD_VOXELCT 3  /* CT system, where lowest X,Y,Z corner is at (0,0,0) cm and slices along Z-axis */
#define OSH_COORD_BZALIGN 4  /* System which aligned so one body corner is at (0,0,0) cm and turned along z-axis */
#define OSH_COORD_BCALIGN 5  /* System which aligned so the body center is at (0,0,0) cm */

#include "transport/osh_transport.h"

int osh_coord_c2v(double const *c, double *v);
int osh_coord_v2c(double const *v, double *c);
int osh_coord_point2sph(double const *v, double *theta, double *phi);
int osh_coord_trans_point(double const p[3], double pt[3], double const t[16]);
int osh_coord_trans_point_hc(double const p[4], double pt[4], double const t[16]);
int osh_coord_trans_pos(struct position const *p, struct position *pt, double const t[16]);
int osh_coord_trans_ray(struct ray_v const *r, struct ray_v *rt, double const t[16]);
int osh_coord_trans_ray_r(
    struct ray const *r, struct ray *rt,
    double const t[16]); // TODO temporary solution, until struct ray_v has been renamed to struct ray

// int osh_coord_trans_step(struct step st, struct step stt, double const t[4][4]);
// int osh_coord_transinv_step(struct step st, struct step stt, double const t[4][4]);

int osh_invert_matrix(const double m[16], double im[16]);

#endif /* !_OSH_COORD */
