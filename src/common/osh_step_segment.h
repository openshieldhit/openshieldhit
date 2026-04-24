#ifndef OSH_STEP_SEGMENT_H
#define OSH_STEP_SEGMENT_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * A step segment is one constant-density piece of a transport step.
 *
 * A transport step advances a particle from its current position to a new
 * position without crossing the current medium boundary.  In the current M5
 * transport path, analytic zones and CT voxel zones both supply exactly one
 * segment: analytic zones use the zone boundary distance, while CT voxel zones
 * use the current voxel-exit distance.
 *
 * The step function (osh_transport_ion_step) receives the segment list
 * produced by the geometry layer before the physics is evaluated.  It uses:
 *   ds  — to accumulate the geometric step length
 *   rho — to accumulate the areal density rho × ds for energy loss,
 *           straggling, and scoring
 *
 * The segment list is built in the wavefront outer loop from GEMCA's boundary
 * distance and material/runtime's density lookup.
 */
struct osh_step_segment {
    double ds;  /* path length through this constant-density piece [cm] */
    double rho; /* actual material density for this piece [g/cm³] */
};

#ifdef __cplusplus
}
#endif

#endif /* OSH_STEP_SEGMENT_H */
