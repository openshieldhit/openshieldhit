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
};

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_RUNTIME_H */
