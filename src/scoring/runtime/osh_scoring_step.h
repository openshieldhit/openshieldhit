#ifndef OSH_SCORING_STEP_H
#define OSH_SCORING_STEP_H

#include "common/osh_step.h"
#include "openshieldhit/status.h"
#include "particle/osh_particle.h"
#include "scoring/runtime/osh_scoring_accumulator.h"
#include "scoring/runtime/osh_scoring_runtime.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Score one transport step into the compiled scoring runtime.
 *
 * @details
 * This is the hot-path entry for step-based scorers.  The current
 * implementation supports Mesh (X,Y,Z) and Cyl (R,Z) geometries and the
 * step-based page kinds implemented in osh_scoring_step.c (ENERGY, FLUENCE,
 * DOSE, DOSEGY, DLET, TLET, DQEFF, TQEFF).
 *
 * @p rt is the read-only compiled descriptor (bin layout, filters, geometry,
 * group ranges); the deposit path never writes through it.  @p acc_set is the
 * mutable accumulator storage to deposit into, indexed in lockstep with
 * @c rt->pages (length @c rt->npages).  @p scratch is the caller-owned per-step
 * voxel-crossing scratch traversal writes into; it must be private to the worker
 * calling this function.  The single-worker serial driver passes the master views
 * (@ref osh_scoring_runtime_master_accumulators and
 * @ref osh_scoring_runtime_master_scratch) so deposits land straight in the shared
 * master pages; a parallel worker passes its own private accumulator set and scratch,
 * and the driver folds the accumulators into the master afterwards with
 * @ref osh_scoring_accumulator_merge.  Because both the deposit target and the
 * traversal scratch are caller-owned, multiple workers may share one read-only
 * @p rt without racing.
 */
enum osh_status osh_scoring_score_step(struct osh_scoring_runtime const *rt,
                                       struct osh_scoring_accumulator *acc_set,
                                       struct osh_scoring_scratch *scratch,
                                       struct particle const *part,
                                       struct step const *st);

/* Point (zero-track-length) energy deposits are scored via
 * osh_scoring_score_point() in osh_scoring_point.h. */

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_STEP_H */
