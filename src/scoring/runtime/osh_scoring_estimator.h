#ifndef OSH_SCORING_ESTIMATOR_H
#define OSH_SCORING_ESTIMATOR_H

#include <stddef.h>

#include "openshieldhit/status.h"
#include "scoring/runtime/osh_scoring_defs.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @file osh_scoring_estimator.h
 * @brief Per-"Quantity" estimator registry: maps a score kind to its
 *        score_step / score_point / postprocess handlers.
 *
 * @details
 * Geometry and estimator are orthogonal.  The geometry decides *which* bins a
 * transport step touches (raytraced crossings for Mesh/Cyl, one bin for Zone/a
 * point).  The estimator decides *what* to book into those bins and *how* to
 * finalise it.  score_step / score_point deposit the raw (extensive) quantity;
 * postprocess turns it into presentation form (÷volume, ratio, unit conversion),
 * once per bin.  Adding a new Quantity = add one handler trio + one registry row.
 * The contract is tabulated in docs/dev/scoring.md.
 */

struct osh_scoring_runtime;
struct osh_scoring_accumulator;
struct osh_scoring_geometry_score_group;
struct osh_scoring_geometry_runtime;
struct osh_scoring_page_runtime;
struct osh_voxel_crossing;
struct particle;
struct step;

/** Deposit handler: books this estimator's raw quantity over @p ncross bin
 *  crossings.
 *
 * This single function type is an implementation convenience shared by the step
 * path and point path.  Step handlers receive raytraced crossings and a physical
 * score length.  Point handlers receive one already-located "crossing" whose
 * idx is the destination bin; point-specific kernels must not infer track
 * semantics from score_len.  If the point path grows more independent, split
 * this into separate step and point function pointer types. */
typedef enum osh_status (*osh_scoring_deposit_fn)(struct osh_scoring_runtime const *rt,
                                                  struct osh_scoring_accumulator *acc_set,
                                                  struct osh_scoring_geometry_score_group const *group,
                                                  struct osh_voxel_crossing const *crossings,
                                                  size_t ncross,
                                                  struct particle const *part,
                                                  struct step const *st,
                                                  double score_len);

/** Postprocess handler: finalise one page's accumulator into presentation form
 *  (@p geo carries the per-bin volume for volume-normalised estimators). */
typedef enum osh_status (*osh_scoring_postprocess_fn)(struct osh_scoring_page_runtime *dst,
                                                      struct osh_scoring_page_runtime const *src,
                                                      struct osh_scoring_geometry_runtime const *geo);

/**
 * One estimator's handler trio.  A NULL handler means the estimator is not active
 * on that path:
 *   score_step  == NULL : not scored on the transport-step path
 *   score_point == NULL : no point meaning (needs a track length, e.g. FLUENCE, LET)
 *   postprocess == NULL : accumulator is already final (e.g. ENERGY, COUNT)
 */
struct osh_scoring_estimator {
    osh_scoring_deposit_fn score_step;
    osh_scoring_deposit_fn score_point;
    osh_scoring_postprocess_fn postprocess;
};

/** Registry lookup by Quantity / score kind.  Returns NULL for an unknown kind. */
struct osh_scoring_estimator const *osh_scoring_estimator_for(enum osh_scoring_score_kind kind);

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_ESTIMATOR_H */
