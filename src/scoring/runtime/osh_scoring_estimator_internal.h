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

/* Shared parameter list of every per-estimator deposit handler.  Point handlers
 * receive one located crossing and ignore score_len; this keeps the registry
 * compact until step and point handler types are split. */
#define OSH_SCORING_DEPOSIT_PARAMS                                                                                     \
    struct osh_scoring_runtime const *rt, struct osh_scoring_accumulator *acc_set,                                     \
        struct osh_scoring_geometry_score_group const *group, struct osh_voxel_crossing const *crossings,              \
        size_t ncross, struct particle const *part, struct step const *st, double score_len

enum osh_status score_step_energy(OSH_SCORING_DEPOSIT_PARAMS);
enum osh_status score_step_fluence(OSH_SCORING_DEPOSIT_PARAMS);
enum osh_status score_step_dose(OSH_SCORING_DEPOSIT_PARAMS);
enum osh_status score_step_dlet(OSH_SCORING_DEPOSIT_PARAMS);
enum osh_status score_step_tlet(OSH_SCORING_DEPOSIT_PARAMS);
enum osh_status score_step_dqeff(OSH_SCORING_DEPOSIT_PARAMS);
enum osh_status score_step_tqeff(OSH_SCORING_DEPOSIT_PARAMS);
enum osh_status score_point_energy(OSH_SCORING_DEPOSIT_PARAMS);
enum osh_status score_point_dose(OSH_SCORING_DEPOSIT_PARAMS);

#endif /* OSH_SCORING_ESTIMATOR_INTERNAL_H */
