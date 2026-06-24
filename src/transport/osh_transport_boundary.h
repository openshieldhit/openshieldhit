#ifndef OSH_TRANSPORT_BOUNDARY_H
#define OSH_TRANSPORT_BOUNDARY_H

/*
 * Shared transport boundary handling.
 *
 * The epsilon is a spatial nudge [cm] used after a particle reaches, or is
 * found sitting on, a geometry boundary.  Keep this transport-owned rather
 * than particle-specific: the floating-point issue is common to all transport
 * families that query GEMCA boundaries.
 */
#define OSH_TRANSPORT_BOUNDARY_EPS 1.0e-8

static inline void osh_transport_nudge_boundary(double *x,
                                                double *y,
                                                double *z,
                                                double ux,
                                                double uy,
                                                double uz) {
    *x += ux * OSH_TRANSPORT_BOUNDARY_EPS;
    *y += uy * OSH_TRANSPORT_BOUNDARY_EPS;
    *z += uz * OSH_TRANSPORT_BOUNDARY_EPS;
}

static inline void osh_transport_advance_to_boundary(double *x,
                                                     double *y,
                                                     double *z,
                                                     double ux,
                                                     double uy,
                                                     double uz,
                                                     double dist) {
    double const d = dist + OSH_TRANSPORT_BOUNDARY_EPS;
    *x += ux * d;
    *y += uy * d;
    *z += uz * d;
}

#endif /* OSH_TRANSPORT_BOUNDARY_H */
