#ifndef OSH_SCORING_STEP_H
#define OSH_SCORING_STEP_H

#include "common/osh_step.h"
#include "openshieldhit/status.h"
#include "particle/osh_particle.h"
#include "scoring/runtime/osh_scoring_runtime.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Score one transport step into the compiled scoring runtime.
 *
 * @details
 * This is the hot-path entry for step-based scorers. The current
 * implementation supports unrotated Cartesian mesh geometries and the
 * `ENERGY` / `FLUENCE` page kinds.
 */
enum osh_status
osh_scoring_score_step(struct osh_scoring_runtime *rt, struct particle const *part, struct step const *st);

/**
 * @brief Score one point event into the compiled scoring runtime.
 *
 * @details
 * Point scoring is not implemented yet.
 */
enum osh_status
osh_scoring_score_point(struct osh_scoring_runtime *rt, struct particle const *part, struct position const *pos);

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_STEP_H */
