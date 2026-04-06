#ifndef OSH_BEAM_MODEL_H
#define OSH_BEAM_MODEL_H

#include "beam/osh_beam.h"
#include "particle/osh_particle.h"
#include "random/osh_rng.h"
#include "transport/osh_transport.h"

/* ---- Primary particle sampling -------------------------------------------
 *
 * osh_beam_new_primaries() is the batched primary-source entry point. It
 * samples n complete primary rays in the UNIVERSE coordinate system.
 *
 * osh_beam_new_primary() is a thin convenience wrapper around the batched API
 * with n=1.
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
 * Coordinate-language convention
 * ------------------------------
 * Beam-model calculations use local upstream/downstream language, not global
 * sign conventions such as "negative z". In the local beam frame:
 *   - +z always means downstream
 *   - the virtual source / scanning-magnet deflection point is always upstream
 *   - SAD is therefore always a positive upstream distance from isocenter
 *
 * This avoids ambiguity when users rotate the beam in UNIVERSE coordinates.
 * All SAD and fan-out calculations are performed in the local beam frame
 * before the final affine transform to UNIVERSE.
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
 *                         SAD is a positive source-to-isocenter distance to an
 *                         upstream virtual source, not a signed z coordinate.
 *                         It is defined relative to isocenter/focal geometry,
 *                         not relative to the beam start plane. BEAMPOS and
 *                         sampled local offsets are therefore included in the
 *                         fan-out angle calculation.
 *
 *   5. apply_transform  — apply the precomputed standard affine matrix
 *                         spot->_tm[16] to map PZALIGN -> UNIVERSE.
 *
 * @param wb       Fully initialised beam workspace.
 * @param rng      Random-number generator state. Required so the caller can
 *                 control reproducibility and thread-local stream ownership.
 *                 Current placeholder sampling may not use it yet, but the
 *                 API is designed around explicit RNG flow.
 * @param n        Number of primaries to sample.
 * @param part_out Receives n particle-species pointers (owned by wb).
 * @param ray_out  Receives n sampled rays in OSH_COORD_UNIVERSE.
 *                 ray_out[i].p[3] holds total kinetic energy [MeV].
 * @return OSH_OK on success, OSH_E* on error. */
int osh_beam_new_primaries(
    struct beam_workspace const *wb, struct osh_rng *rng, size_t n, struct particle **part_out, struct ray_v *ray_out);

int osh_beam_new_primary(struct beam_workspace const *wb,
                         struct osh_rng *rng,
                         struct particle **part_out,
                         struct ray_v *ray_out);

#endif /* OSH_BEAM_MODEL_H */
