#include "scoring/runtime/osh_scoring_estimator.h"

#include "scoring/runtime/osh_scoring_postprocess.h"   /* postprocess_* handlers   */
#include "scoring/runtime/osh_scoring_step_internal.h" /* score_step_* / score_point_* */

/*
 * The estimator table.  One entry per Quantity; adding support for a new one is a
 * single row here plus its handler(s).  Kept as a switch (rather than a
 * kind-indexed array) for portability and because score kinds are sparse.
 *
 * See docs/dev/scoring.md for the tabulated score_step / score_point / postprocess
 * contract that this table implements.
 */
static struct osh_scoring_estimator const k_energy = {score_step_energy, score_step_energy, NULL};
static struct osh_scoring_estimator const k_fluence = {score_step_fluence, NULL, postprocess_volume};
static struct osh_scoring_estimator const k_dose = {score_step_dose, score_point_dose, postprocess_volume};
static struct osh_scoring_estimator const k_dosegy = {score_step_dose, score_point_dose, postprocess_dosegy};
static struct osh_scoring_estimator const k_dlet = {score_step_dlet, NULL, postprocess_ratio};
static struct osh_scoring_estimator const k_tlet = {score_step_tlet, NULL, postprocess_ratio};
static struct osh_scoring_estimator const k_dqeff = {score_step_dqeff, NULL, postprocess_ratio};
static struct osh_scoring_estimator const k_tqeff = {score_step_tqeff, NULL, postprocess_ratio};
static struct osh_scoring_estimator const k_nkerma = {NULL, NULL, postprocess_volume};

struct osh_scoring_estimator const *osh_scoring_estimator_for(enum osh_scoring_score_kind kind) {
    switch (kind) {
    case OSH_SCORING_SCORE_ENERGY:
        return &k_energy;
    case OSH_SCORING_SCORE_FLUENCE:
        return &k_fluence;
    case OSH_SCORING_SCORE_DOSE:
        return &k_dose;
    case OSH_SCORING_SCORE_DOSEGY:
        return &k_dosegy;
    case OSH_SCORING_SCORE_DLET:
        return &k_dlet;
    case OSH_SCORING_SCORE_TLET:
        return &k_tlet;
    case OSH_SCORING_SCORE_DQEFF:
        return &k_dqeff;
    case OSH_SCORING_SCORE_TQEFF:
        return &k_tqeff;
    case OSH_SCORING_SCORE_NKERMA:
        return &k_nkerma;
    default:
        return NULL;
    }
}
