#include "scoring/runtime/osh_scoring_estimator.h"

#include "scoring/runtime/osh_scoring_estimator_internal.h" /* estimator handler declarations */
#include "scoring/runtime/osh_scoring_postprocess.h"        /* postprocess_* handlers   */

/*
 * The estimator table.  One entry per Quantity; adding support for a new one is a
 * single row here plus its handler(s).  Kept as a switch (rather than a
 * kind-indexed array) for portability and because score kinds are sparse.
 *
 * See docs/dev/scoring.md for the tabulated score_step / score_point / postprocess
 * contract that this table implements.
 */
static struct osh_scoring_estimator const k_energy = {
    osh_scoring_estimator_step_energy, osh_scoring_estimator_point_energy, NULL};
static struct osh_scoring_estimator const k_fluence = {osh_scoring_estimator_step_fluence, NULL, postprocess_volume};
static struct osh_scoring_estimator const k_dose = {
    osh_scoring_estimator_step_dose, osh_scoring_estimator_point_dose, postprocess_volume};
static struct osh_scoring_estimator const k_dosegy = {
    osh_scoring_estimator_step_dose, osh_scoring_estimator_point_dose, postprocess_dosegy};
static struct osh_scoring_estimator const k_dirtydose = {
    osh_scoring_estimator_step_dirtydose, osh_scoring_estimator_point_dirtydose, postprocess_volume};
static struct osh_scoring_estimator const k_dirtydosegy = {
    osh_scoring_estimator_step_dirtydose, osh_scoring_estimator_point_dirtydose, postprocess_dosegy};
static struct osh_scoring_estimator const k_dlet = {
    osh_scoring_estimator_step_dlet, osh_scoring_estimator_point_dlet, postprocess_ratio};
static struct osh_scoring_estimator const k_tlet = {
    osh_scoring_estimator_step_tlet, osh_scoring_estimator_point_tlet, postprocess_ratio};
static struct osh_scoring_estimator const k_dqeff = {osh_scoring_estimator_step_dqeff, NULL, postprocess_ratio};
static struct osh_scoring_estimator const k_tqeff = {osh_scoring_estimator_step_tqeff, NULL, postprocess_ratio};
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
    case OSH_SCORING_SCORE_DIRTYDOSE:
        return &k_dirtydose;
    case OSH_SCORING_SCORE_DIRTYDOSEGY:
        return &k_dirtydosegy;
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
