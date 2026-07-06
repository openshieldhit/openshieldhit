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
 * Map a transport step's zone id to the dense Zone-scorer bin index.  Linear scan
 * over the (small) selected-zone list; returns 1 and sets @p *idx_out on a hit, 0
 * when @p zone is negative or is not one of this geometry's selected zones.
 */
int zone_bin_index(struct osh_scoring_geometry_runtime const *geo, int zone, size_t *idx_out);

/**
 * Build a raytrace grid from a Cyl (R,Z) scoring geometry.  Field convention:
 * origin/spacing/n[0] = r_min/dr/nr; origin/spacing/n[2] = z_min/dz/nz; index [1]
 * is unused (n[1] = 1).  Per-voxel 1/V is geo->cyl_vol_inv[r_bin] with
 * r_bin = flat_idx % nr and flat_idx = z_bin * nr + r_bin.
 */
enum osh_status cyl_geometry_to_grid(struct osh_scoring_geometry_runtime const *geo, struct osh_raytrace_grid *grid);

/**
 * Deposit st->de into the ENERGY pages of one score group over @p crossings.
 * Per crossing the contribution is st->de * (crossing.path_len / score_len).
 */
enum osh_status score_step_energy(struct osh_scoring_runtime const *rt,
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
enum osh_status score_step_dose(struct osh_scoring_runtime const *rt,
                                struct osh_scoring_accumulator *acc_set,
                                struct osh_scoring_geometry_score_group const *group,
                                struct osh_voxel_crossing const *crossings,
                                size_t ncross,
                                struct particle const *part,
                                struct step const *st,
                                double score_len);

/* Shared parameter list of every per-estimator deposit handler (step and point),
 * so a new estimator declares/defines its handler with one token. */
#define OSH_SCORING_DEPOSIT_PARAMS                                                                                     \
    struct osh_scoring_runtime const *rt, struct osh_scoring_accumulator *acc_set,                                     \
        struct osh_scoring_geometry_score_group const *group, struct osh_voxel_crossing const *crossings,              \
        size_t ncross, struct particle const *part, struct step const *st, double score_len

/* The remaining per-estimator deposit handlers share the same signature and are
 * registered in osh_scoring_estimator.c.  FLUENCE deposits track length; the LET
 * and Qeff handlers fill a two-pass (data, data2) accumulator finalised as a ratio
 * in postprocess.  score_point_dose is the point counterpart of score_step_dose
 * (skips neutral particles, which deposit energy but no local dose). */
enum osh_status score_step_fluence(OSH_SCORING_DEPOSIT_PARAMS);
enum osh_status score_step_dlet(OSH_SCORING_DEPOSIT_PARAMS);
enum osh_status score_step_tlet(OSH_SCORING_DEPOSIT_PARAMS);
enum osh_status score_step_dqeff(OSH_SCORING_DEPOSIT_PARAMS);
enum osh_status score_step_tqeff(OSH_SCORING_DEPOSIT_PARAMS);
enum osh_status score_point_dose(OSH_SCORING_DEPOSIT_PARAMS);

#endif /* OSH_SCORING_STEP_INTERNAL_H */
