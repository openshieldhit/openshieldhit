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
    struct osh_voxel_crossing *crossing_buf;       /* reusable per-step scratch buffer */
    size_t crossing_cap;                           /* allocated length of crossing_buf */
    struct osh_material_runtime const *mat_tables; /* SP tables; NULL until wired by simulation */
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

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_RUNTIME_H */
