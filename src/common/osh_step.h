#ifndef OSH_STEP_H
#define OSH_STEP_H

#include "common/osh_coord.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Simulation interface types: position and step.
 *
 * These structs are the handshake between the transport engine and domain
 * modules (scoring, gemca, beam). They live in common/ so any module can
 * receive them without depending on transport/.
 *
 * Field ordering follows the project convention: largest/most-aligned fields
 * first, reducing struct padding.
 */

/**
 * @brief A point in space with direction, energy, and transport context.
 *
 * @details
 * Represents a particle's state at a single location: where it is, where it
 * is going, how much kinetic energy it has, and the transport-engine context
 * (medium, zone, density) at that point.
 *
 * The coordinate system is identified by the @ref system field using the
 * OSH_COORD_* constants from osh_coord.h.
 *
 * TODO: this struct is currently named position for historical reasons but is
 * functionally a ray_v with transport metadata attached.  Consider renaming
 * once the full transport refactor is settled.
 */
struct position {
    double p[4]; /* x,y,z,E — E is total kinetic energy [MeV] */
    double v[3]; /* unit direction vector */
    double rho;  /* CT-corrected density at this point [g/cm3] */
    int medium;  /* medium ID at this point, -1 if unknown */
    int zone;    /* zone number at this point, -1 if unknown */
    int system;  /* coordinate system identifier (OSH_COORD_*) */
};

/**
 * @brief One transport step: start, stop, energy loss, and transport context.
 *
 * @details
 * Records a single particle step from p to q along direction v, with scalar
 * track length ds and energy deposit de.  The transport context (rho, medium,
 * zone) describes the medium the step was taken in; the step must lie entirely
 * within one zone.
 *
 * Sign convention: de is positive for energy loss, negative for energy gain.
 *
 * TODO: consider refactoring to hold two struct ray_v (entry and exit rays)
 * instead of p/q/v scalars, which would preserve both incoming and outgoing
 * directions at zone boundaries while keeping all transport metadata together.
 */
struct step {
    double p[4]; /* start: x,y,z,E — E is total kinetic energy [MeV] */
    double q[4]; /* stop:  x,y,z,E */
    double v[3]; /* entry direction (unit vector at p) */
    double w[3]; /* exit direction  (unit vector at q); equals v for straight steps */
    double ds;   /* track length of this step [cm] */
    double de;   /* energy loss (positive) or gain (negative) [MeV] */
    double rho;  /* CT-corrected density [g/cm3] */
    int medium;  /* medium ID, -1 if unknown */
    int zone;    /* zone number, -1 if unknown */
    int system;  /* coordinate system identifier (OSH_COORD_*) */
};

void osh_step_print(struct step const *st);
void osh_step_copy(struct step *dst, struct step const *src);

void osh_position_print(struct position const *pos);
void osh_position_copy(struct position *dst, struct position const *src);

#ifdef __cplusplus
}
#endif

#endif /* OSH_STEP_H */
