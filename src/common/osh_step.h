#ifndef OSH_STEP_H
#define OSH_STEP_H

#include <stddef.h>
#include <stdint.h>

#include "common/osh_coord.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Simulation interface types: position and step.
 *
 * These structs are the handshake between the transport engine and domain
 * modules (scoring, gemca, beam).  They live in common/ so any module can
 * include them without depending on transport/.
 *
 * Field ordering follows the project convention: largest/most-aligned fields
 * first, minimising struct padding.
 */

/**
 * @brief A point in space with direction, energy, and transport context.
 *
 * @details
 * Represents a particle's instantaneous state: where it is, where it is
 * going, its kinetic energy, and the local transport context (material,
 * zone, density).  Used as an intermediate result inside the transport loop
 * and as input to geometry queries.
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
    double rho;  /* CT-corrected density at this point [g/cm³] */
    int medium;  /* material index at this point; -1 if unknown */
    int zone;    /* geometry zone index at this point; -1 if unknown */
    int system;  /* coordinate system (OSH_COORD_*) */
};

/**
 * @brief One completed transport step, carrying geometry, physics, and
 *        history context for scoring.
 *
 * @details
 * Produced by the transport engine for every step taken by a particle.
 * Passed to osh_scoring_score_step() and any other consumer that needs
 * per-step information.
 *
 * The struct is intentionally self-contained: a consumer needs no other
 * argument to answer "where did this step happen, how much energy was
 * deposited, and which history does it belong to?".
 *
 * Geometry: the step runs from point p to point q along direction v inside
 * a single zone.  For straight-line (CSDA) transport v == w; once multiple
 * Coulomb scattering is added, w will diverge from v at zone boundaries.
 *
 * Physics: de is the energy deposited in the medium over track length ds.
 * Sign convention: de > 0 for energy loss (the normal case), de < 0 for
 * energy gain (e.g. delta-electron backscatter, currently not implemented).
 *
 * History context: wt, gen, and prim_idx are per-history metadata supplied
 * by the transport engine from the particle pool.  They are NOT species
 * properties — two protons in the same beam can have different gen or wt.
 * Scoring filters (GEN, NPRIM) read these fields rather than struct particle
 * so that struct particle can remain a cold, shared, const species descriptor.
 *
 * wt (statistical weight): 1.0 for unweighted transport.  Set by beam spot
 * selection (SOBP weights) and modified by variance-reduction techniques.
 *
 * gen (generation): 0 for beam primaries, 1 for their direct secondaries
 * (e.g. nuclear-reaction products), 2 for the next generation, and so on.
 * Replaces the TREE branch-depth concept from the legacy SHIELD-HIT code
 * without requiring an explicit cascade tree structure.
 *
 * prim_idx (primary ancestor index): 0-based index into the beam batch
 * identifying which original primary particle spawned this history.
 * Enables per-primary dose tallies and correlated variance estimation.
 * Filled by the transport engine from the particle pool, which in turn
 * receives it from osh_beam_runtime_fill_pool().
 *
 * TODO: consider refactoring to hold two struct ray_v (entry and exit rays)
 * instead of p/q/v/w scalars to preserve both directions at zone boundaries.
 */
struct step {
    double p[4];       /* start point: x,y,z [cm] and total KE [MeV] */
    double q[4];       /* stop  point: x,y,z [cm] and total KE [MeV] */
    double v[3];       /* entry direction unit vector at p */
    double w[3];       /* exit  direction unit vector at q; equals v for straight steps */
    double ds;         /* track length [cm] */
    double de;         /* energy deposit [MeV]; positive = loss, negative = gain */
    double rho;        /* local material density [g/cm³] */
    double wt;         /* statistical weight of this history; 1.0 = unweighted */
    size_t voxel_idx;  /* flat CT voxel index; valid when zone is inside a voxel body */
    int medium;        /* material index; -1 if unknown */
    int zone;          /* geometry zone index; -1 if unknown */
    int system;        /* coordinate system (OSH_COORD_*) */
    uint32_t prim_idx; /* 0-based index of the beam-primary ancestor in the batch */
    uint8_t gen;       /* generation: 0 = beam primary, 1 = first secondary, … */
    uint8_t has_voxel; /* non-zero when voxel_idx is valid */
};

void osh_step_print(struct step const *st);
void osh_step_copy(struct step *dst, struct step const *src);

void osh_position_print(struct position const *pos);
void osh_position_copy(struct position *dst, struct position const *src);

#ifdef __cplusplus
}
#endif

#endif /* OSH_STEP_H */
