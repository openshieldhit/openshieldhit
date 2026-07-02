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
 * @brief Completeness label to stamp on a saved result — never NULL.
 *
 * @details
 * Resolves @ref osh_scoring_runtime::completeness to a printable token for the
 * save layer, collapsing the NULL default (a fully family-drained result) to the
 * literal "exact".  Used by the ASCII and BDO writers so the on-disk honesty
 * label is well-defined for every output, dump or final.
 */
static inline char const *osh_scoring_runtime_completeness_label(struct osh_scoring_runtime const *rt) {
    return (rt && rt->completeness) ? rt->completeness : "exact";
}

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_RUNTIME_H */
