#ifndef OSH_SCORING_RUNTIME_H
#define OSH_SCORING_RUNTIME_H

#include <stddef.h>

#include "common/raytrace/osh_raytrace.h"
#include "material/runtime/osh_material_runtime.h"
#include "scoring/runtime/osh_scoring_filter_runtime.h"
#include "scoring/runtime/osh_scoring_geometry_runtime.h"
#include "scoring/runtime/osh_scoring_output_runtime.h"
#include "scoring/runtime/osh_scoring_settings_runtime.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---- Per-worker traversal scratch ---------------------------------------- */

/**
 * @brief Caller-owned voxel-crossing scratch for one @ref osh_scoring_score_step call.
 *
 * @details
 * Holds the reusable per-step buffer that raytrace / voxel-crossing traversal
 * writes into before the per-group deposits consume it.  It is mutated on every
 * step, so it must be private to the worker doing the traversal: the serial
 * single-worker driver hands the long-lived @ref osh_scoring_runtime.master_scratch
 * (built once at compile time, see @ref osh_scoring_runtime_master_scratch), while a
 * parallel worker passes its own instance.  Keeping this off the read-only compiled
 * descriptor is what lets multiple workers share one @c rt without racing on the
 * traversal scratch.
 *
 * @c crossing_buf is grown to fit by the owner (the serial scratch is pre-sized at
 * compile time to the largest geometry, so the serial hot path never reallocates).
 */
struct osh_scoring_scratch {
    struct osh_voxel_crossing *crossing_buf; /* reusable per-step voxel-crossing scratch buffer */
    size_t crossing_cap;                     /* allocated length of crossing_buf */
};

/* ---- Deposit target ------------------------------------------------------ */

struct osh_scoring_accumulator;

/**
 * @brief Where a transport worker deposits its scoring: an accumulator set plus
 *        the per-step traversal scratch, threaded together as one unit.
 *
 * @details
 * The pair @ref osh_scoring_score_step / @ref osh_scoring_score_point need in
 * lockstep: @p acc_set is the mutable accumulator storage (indexed alongside
 * @c rt->pages, length @c rt->npages) and @p scratch is the caller-owned
 * voxel-crossing scratch the traversal writes into.  Bundling them means the
 * transport call chain threads a single small struct instead of two parallel
 * pointers through several signatures — and it is the natural unit a worker owns.
 *
 * The serial single-worker driver points this at the shared master views
 * (@ref osh_scoring_runtime_master_accumulators / @ref osh_scoring_runtime_master_scratch);
 * a parallel or replica worker points it at its own private set, and the driver
 * folds the accumulators into the master afterwards with
 * @ref osh_scoring_accumulator_merge.  A NULL @p acc_set / @p scratch means "use
 * the master view" so a call site can fall back to today's exact behaviour.
 */
struct osh_score_target {
    struct osh_scoring_accumulator *acc_set; /* deposit target (npages-long); NULL ⇒ master */
    struct osh_scoring_scratch *scratch;     /* per-worker traversal scratch; NULL ⇒ master */
};

/* ---- Top-level compiled runtime ------------------------------------------ */

/**
 * @brief Scoring-owned compiled scorer workspace.
 *
 * @details
 * The prepare step resolves raw `detect.dat` names to dense indices and
 * allocates page-local accumulators. Runtime traversal is geometry-centric:
 * each geometry owns one contiguous page span so Siddon/Jacobs traversal can
 * be executed once per geometry, then all attached pages updated in one loop.
 *
 * Output files remain a separate cold-path grouping. Each output stores the
 * page indices that should later be written into one BDO file/page sequence.
 */
struct osh_scoring_runtime {
    struct osh_scoring_filter_runtime *filters;
    struct osh_scoring_settings_runtime *settings;
    struct osh_scoring_geometry_runtime *geometries;
    struct osh_scoring_page_runtime *pages;
    struct osh_scoring_output_runtime *outputs;
    size_t nfilters;
    size_t nsettings;
    size_t ngeometries;
    size_t npages;
    size_t noutputs;
    /* Serial single-worker traversal scratch: the long-lived per-step voxel-crossing
     * buffer the serial driver passes to osh_scoring_score_step() (see
     * osh_scoring_runtime_master_scratch).  Pre-sized once by osh_scoring_compile() to
     * the largest geometry, so the serial hot path never reallocates.  Not consulted
     * by the deposit path itself: score_step uses the caller-supplied scratch, so a
     * parallel worker passes its own and this stays the serial driver's private copy. */
    struct osh_scoring_scratch master_scratch;
    struct osh_material_runtime const *mat_tables; /* SP tables; NULL until wired by simulation */
    /* Partial-result completeness label for the save layer (issue #193 / #195).
     * NULL (the default after compile) means a physically complete, family-drained
     * result — the writers stamp "exact".  A mid-run dump of an ion-only or
     * otherwise family-incomplete snapshot sets this on its shadow view (never on
     * the live runtime) to e.g. "families_pending=neutron" so the honesty label
     * on disk (BDO OSHBDO_RT_COMPLETENESS / ASCII "# COMPLETENESS:") reflects it.
     * Borrowed: points at a static string, never owned/freed by the runtime. */
    char const *completeness;
    /* Master accumulator view: an npages-long array whose element i shallow-aliases
     * pages[i].acc (same data pointers).  Built once by osh_scoring_compile() so the
     * single-worker serial driver can hand osh_scoring_score_step() a deposit target
     * that is identical in shape to a parallel worker's private set, with no per-step
     * allocation.  NULL when npages == 0.  Owns only the array; the arrays it points
     * into are owned by the pages. */
    struct osh_scoring_accumulator *master_acc;
    /* Set once by osh_scoring_postprocess_into() after a successful in-place
     * (dst == src) finalisation — the path both osh_scoring_postprocess() and a
     * direct in-place _into(rt, rt) take; guards against a second destructive
     * in-place call (which would double-divide by volume / re-ratio data÷data2).
     * Zero after osh_scoring_compile() (which memsets the runtime). The repeatable,
     * non-destructive primitive for dump/checkpoint is osh_scoring_postprocess_into()
     * with a distinct dst, which never touches this flag. See docs/dev/scoring.md §6. */
    int postprocessed;
};

/**
 * @brief Master accumulator-set view (length npages) aliasing each page's storage.
 *
 * @details
 * The serial single-worker driver passes this to @ref osh_scoring_score_step so its
 * deposits land straight in the shared master pages, while a parallel worker passes
 * its own private set instead.  Built once by @ref osh_scoring_compile; never
 * allocated on the hot path.  Returns NULL when there are no pages.
 */
static inline struct osh_scoring_accumulator *osh_scoring_runtime_master_accumulators(struct osh_scoring_runtime *rt) {
    return rt ? rt->master_acc : NULL;
}

/**
 * @brief Serial single-worker traversal scratch, pre-sized for the largest geometry.
 *
 * @details
 * The serial driver passes this to @ref osh_scoring_score_step as its per-step
 * voxel-crossing scratch, mirroring @ref osh_scoring_runtime_master_accumulators for the
 * deposit target.  A parallel worker passes its own @ref osh_scoring_scratch instead, so
 * the compiled @c rt is never the source of the mutable traversal scratch.  Built once by
 * @ref osh_scoring_compile; never allocated on the hot path.
 */
static inline struct osh_scoring_scratch *osh_scoring_runtime_master_scratch(struct osh_scoring_runtime *rt) {
    return rt ? &rt->master_scratch : NULL;
}

/**
 * @brief Resolve a deposit target's accumulator set, falling back to the master.
 *
 * @details
 * Returns @p t->acc_set when a worker supplied a private set, else the shared
 * master view.  A NULL @p t (no target threaded) resolves to the master too, so
 * the serial path is bit-for-bit today's behaviour.  Single seam for the
 * "private-or-master" decision so no call site hard-codes the master again.
 */
static inline struct osh_scoring_accumulator *osh_score_target_accumulators(struct osh_score_target const *t,
                                                                            struct osh_scoring_runtime *rt) {
    return (t && t->acc_set) ? t->acc_set : osh_scoring_runtime_master_accumulators(rt);
}

/**
 * @brief Resolve a deposit target's traversal scratch, falling back to the master.
 *
 * @details
 * The scratch counterpart to @ref osh_score_target_accumulators: @p t->scratch
 * when supplied, else the master scratch (also the NULL-@p t case).
 */
static inline struct osh_scoring_scratch *osh_score_target_scratch(struct osh_score_target const *t,
                                                                   struct osh_scoring_runtime *rt) {
    return (t && t->scratch) ? t->scratch : osh_scoring_runtime_master_scratch(rt);
}

/* Batch-means observations used to derive a checkpoint cadence for a plain
 * single-threaded run when at least one page tracks variance (issue #209).  Enough
 * degrees of freedom for a stable error estimate without paying many merge
 * boundaries.  Kept internal: the run-wide batch count is deliberately not exposed
 * as a card or flag yet (see PR #247 discussion) — an explicit dump cadence or
 * --score-replicas already provides the batch split when the user needs control. */
#define OSH_SCORING_VARIANCE_DEFAULT_BATCHES 10

/**
 * @brief Does this runtime track per-bin Monte-Carlo variance (Welford M2) on any page?
 *
 * @details
 * True when at least one page has a "Variance On" Settings block, so @ref
 * osh_scoring_compile allocated the M2 arrays for it (issue #209).  Variance is
 * enabled per estimator, so this scans every page rather than assuming page 0 is
 * representative.  The transport driver uses this once at setup to decide whether to
 * score each checkpoint batch into a private set and fold it (with the
 * Schubert–Gertz cross-term) instead of depositing cumulatively into the master;
 * the finalize/save layers use it to emit the standard-error columns.  When no page
 * tracks variance, the accumulators and the hot path are exactly as before.
 */
static inline int osh_scoring_runtime_tracks_variance(struct osh_scoring_runtime const *rt) {
    size_t p;

    if (!rt || rt->npages == 0u || !rt->pages) {
        return 0;
    }
    for (p = 0u; p < rt->npages; ++p) {
        if (rt->pages[p].acc.data_var != NULL) {
            return 1;
        }
    }
    return 0;
}

/**
 * @brief Completeness label to stamp on a saved result — never NULL.
 *
 * @details
 * Resolves @ref osh_scoring_runtime::completeness to a printable token for the
 * save layer, collapsing the NULL default (a fully family-drained result) to the
 * literal "exact".  Used by the ASCII and BDO writers so the on-disk honesty
 * label is well-defined for every output, dump or final.
 */
static inline char const *osh_scoring_runtime_completeness_label(struct osh_scoring_runtime const *rt) {
    /* A NULL completeness (the default after compile) means a fully family-drained
     * result, which the writers stamp as the literal "exact". */
    if (rt && rt->completeness) {
        return rt->completeness;
    }
    return "exact";
}

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_RUNTIME_H */
