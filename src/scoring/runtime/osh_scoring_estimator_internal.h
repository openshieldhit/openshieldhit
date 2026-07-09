#ifndef OSH_SCORING_ESTIMATOR_INTERNAL_H
#define OSH_SCORING_ESTIMATOR_INTERNAL_H

/*
 * Module-internal estimator handler declarations.
 *
 * The registry in osh_scoring_estimator.c maps score kinds to these handlers.
 * Geometry drivers locate crossings or bins first; estimator handlers only book
 * values into the already-selected accumulator pages.
 */

#include "common/osh_step.h"
#include "common/raytrace/osh_raytrace.h"
#include "openshieldhit/status.h"
#include "particle/osh_particle.h"
#include "scoring/runtime/osh_scoring_accumulator.h"
#include "scoring/runtime/osh_scoring_runtime.h"

enum osh_status osh_scoring_estimator_step_energy(struct osh_scoring_runtime const *rt,
                                                  struct osh_scoring_accumulator *acc_set,
                                                  struct osh_scoring_geometry_score_group const *group,
                                                  struct osh_voxel_crossing const *crossings,
                                                  size_t ncross,
                                                  struct particle const *part,
                                                  struct step const *st,
                                                  double score_len);

enum osh_status osh_scoring_estimator_step_fluence(struct osh_scoring_runtime const *rt,
                                                   struct osh_scoring_accumulator *acc_set,
                                                   struct osh_scoring_geometry_score_group const *group,
                                                   struct osh_voxel_crossing const *crossings,
                                                   size_t ncross,
                                                   struct particle const *part,
                                                   struct step const *st,
                                                   double score_len);

enum osh_status osh_scoring_estimator_step_dose(struct osh_scoring_runtime const *rt,
                                                struct osh_scoring_accumulator *acc_set,
                                                struct osh_scoring_geometry_score_group const *group,
                                                struct osh_voxel_crossing const *crossings,
                                                size_t ncross,
                                                struct particle const *part,
                                                struct step const *st,
                                                double score_len);

enum osh_status osh_scoring_estimator_step_dirtydose(struct osh_scoring_runtime const *rt,
                                                     struct osh_scoring_accumulator *acc_set,
                                                     struct osh_scoring_geometry_score_group const *group,
                                                     struct osh_voxel_crossing const *crossings,
                                                     size_t ncross,
                                                     struct particle const *part,
                                                     struct step const *st,
                                                     double score_len);

enum osh_status osh_scoring_estimator_step_dlet(struct osh_scoring_runtime const *rt,
                                                struct osh_scoring_accumulator *acc_set,
                                                struct osh_scoring_geometry_score_group const *group,
                                                struct osh_voxel_crossing const *crossings,
                                                size_t ncross,
                                                struct particle const *part,
                                                struct step const *st,
                                                double score_len);

enum osh_status osh_scoring_estimator_step_tlet(struct osh_scoring_runtime const *rt,
                                                struct osh_scoring_accumulator *acc_set,
                                                struct osh_scoring_geometry_score_group const *group,
                                                struct osh_voxel_crossing const *crossings,
                                                size_t ncross,
                                                struct particle const *part,
                                                struct step const *st,
                                                double score_len);

enum osh_status osh_scoring_estimator_step_dqeff(struct osh_scoring_runtime const *rt,
                                                 struct osh_scoring_accumulator *acc_set,
                                                 struct osh_scoring_geometry_score_group const *group,
                                                 struct osh_voxel_crossing const *crossings,
                                                 size_t ncross,
                                                 struct particle const *part,
                                                 struct step const *st,
                                                 double score_len);

enum osh_status osh_scoring_estimator_step_tqeff(struct osh_scoring_runtime const *rt,
                                                 struct osh_scoring_accumulator *acc_set,
                                                 struct osh_scoring_geometry_score_group const *group,
                                                 struct osh_voxel_crossing const *crossings,
                                                 size_t ncross,
                                                 struct particle const *part,
                                                 struct step const *st,
                                                 double score_len);

enum osh_status osh_scoring_estimator_point_energy(struct osh_scoring_runtime const *rt,
                                                   struct osh_scoring_accumulator *acc_set,
                                                   struct osh_scoring_geometry_score_group const *group,
                                                   size_t spatial_idx,
                                                   struct particle const *part,
                                                   struct step const *st);

enum osh_status osh_scoring_estimator_point_dose(struct osh_scoring_runtime const *rt,
                                                 struct osh_scoring_accumulator *acc_set,
                                                 struct osh_scoring_geometry_score_group const *group,
                                                 size_t spatial_idx,
                                                 struct particle const *part,
                                                 struct step const *st);

enum osh_status osh_scoring_estimator_point_dirtydose(struct osh_scoring_runtime const *rt,
                                                      struct osh_scoring_accumulator *acc_set,
                                                      struct osh_scoring_geometry_score_group const *group,
                                                      size_t spatial_idx,
                                                      struct particle const *part,
                                                      struct step const *st);

#endif /* OSH_SCORING_ESTIMATOR_INTERNAL_H */
