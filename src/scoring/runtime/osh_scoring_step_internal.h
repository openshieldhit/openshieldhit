#ifndef OSH_SCORING_STEP_INTERNAL_H
#define OSH_SCORING_STEP_INTERNAL_H

/*
 * Module-internal deposit helpers shared between the step-based scorer
 * (osh_scoring_step.c) and the point-based scorer (osh_scoring_point.c).
 *
 * NOT part of the public scoring API — do not include outside src/scoring.
 * These carry the filter, differential-axis and dose-override logic so the
 * point scorer can reuse it instead of duplicating it.
 */

#include "common/osh_step.h"
#include "common/raytrace/osh_raytrace.h"
#include "openshieldhit/status.h"
#include "particle/osh_particle.h"
#include "scoring/runtime/osh_scoring_accumulator.h"
#include "scoring/runtime/osh_scoring_runtime.h"

/** Build a raytrace grid + 1/voxel-volume from a Mesh scoring geometry. */
enum osh_status mesh_geometry_to_grid(struct osh_scoring_geometry_runtime const *geo,
                                      struct osh_raytrace_grid *grid,
                                      double *voxel_volume_inv_out);

/**
 * Deposit st->de into the ENERGY pages of one score group over @p crossings.
 * Per crossing the contribution is st->de * (crossing.path_len / score_len).
 */
enum osh_status score_group_energy(struct osh_scoring_runtime const *rt,
                                   struct osh_scoring_accumulator *acc_set,
                                   struct osh_scoring_geometry_score_group const *group,
                                   struct osh_voxel_crossing const *crossings,
                                   size_t ncross,
                                   struct particle const *part,
                                   struct step const *st,
                                   double score_len);

/**
 * Deposit st->de into the DOSE pages of one score group over @p crossings.
 * Per crossing the contribution is crossing.path_len * crossing.vol_inv *
 * st->de / (score_len * st->rho) [MeV/g], with per-page stopping-power-ratio
 * overrides (dose-to-medium) applied.
 */
enum osh_status score_group_dose(struct osh_scoring_runtime const *rt,
                                 struct osh_scoring_accumulator *acc_set,
                                 struct osh_scoring_geometry_score_group const *group,
                                 struct osh_voxel_crossing const *crossings,
                                 size_t ncross,
                                 struct particle const *part,
                                 struct step const *st,
                                 double score_len);

#endif /* OSH_SCORING_STEP_INTERNAL_H */
