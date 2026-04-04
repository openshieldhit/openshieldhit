#ifndef OSH_BEAM_MODEL_H
#define OSH_BEAM_MODEL_H

#include "beam/osh_beam.h"
#include "particle/osh_particle.h"
#include "transport/osh_transport.h"

/* ---- Primary particle sampling -------------------------------------------
 *
 * osh_beam_new_primary() is the single entry point called once per history.
 * It samples a complete primary ray in the UNIVERSE coordinate system.
 *
 * Spot position convention
 * ------------------------
 * All spot positions (spot->p[x,y,z]) are beam-local offsets relative to
 * isocenter BEFORE gantry/table rotation. In other words, BEAMPOS lives in
 * the beam-local PZALIGN frame. The post-parse step folds this offset into the
 * precomputed spot->_tm matrix as a standard affine transform:
 *
 *   p_universe = R * p_local + R * spot->p
 *
 * where R is the beam direction rotation derived from shared.theta/phi.
 *
 * Sampling pipeline
 * -----------------
 *   1. select_spot      — pick the active beam_spot (spot[0] for single-spot;
 *                         weighted draw from the spotlist for SOBP mode).
 *
 *   2. sample_energy    — draw total kinetic energy T around spot->t0 with
 *                         spread spot->tsigma.
 *
 *   3. sample_phasespace — sample the transverse phase space in the beam-local
 *                         PZALIGN frame (beam along +Z).
 *
 *   4. apply_sad        — if shared.use_sad: apply source-axis-distance
 *                         correction (scanning-magnet fan-out) in PZALIGN.
 *
 *   5. apply_transform  — apply the precomputed standard affine matrix
 *                         spot->_tm[16] to map PZALIGN -> UNIVERSE.
 *
 * @param wb       Fully initialised beam workspace.
 * @param part_out Receives a pointer to the particle species (owned by wb).
 * @param ray_out  Receives the sampled ray in OSH_COORD_UNIVERSE.
 *                 ray_out->p[3] holds total kinetic energy [MeV].
 * @return OSH_OK on success, negative OSH_E* on error. */
int osh_beam_new_primary(struct beam_workspace const *wb, struct particle **part_out, struct ray_v *ray_out);

#endif /* OSH_BEAM_MODEL_H */
