#ifndef OSH_SCORING_SAVE_PLOT_H
#define OSH_SCORING_SAVE_PLOT_H

#include "openshieldhit/status.h"
#include "scoring/save/osh_scoring_save.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Save one compiled output as a native, dependency-free SVG line plot.
 *
 * @details
 * Prototype for issue #238 — a zero-dependency "quick-look" plot so a run can
 * drop a viewable Bragg curve / 1-D profile next to its `.bdo`/`.dat` on a
 * machine that only has the `openshieldhit` binary (no Python/matplotlib).  The
 * whole feature is compiled in only when the project is configured with
 * `-DOSH_ENABLE_PLOT=ON`; the stock binary carries no plotting code.
 *
 * Scope of this prototype is deliberately narrow — SVG only, 1-D spatial
 * profiles only:
 *
 *   - MESH with exactly one non-singleton spatial axis (the profile axis), e.g.
 *     a depth-dose / Bragg curve scored on a `1 x 1 x N` mesh.
 *   - CYL with exactly one non-singleton axis (R or Z).
 *
 * Every page of the output is drawn as its own polyline (auto-scaled to a shared
 * y-axis) with a small legend, so a single-quantity output yields one clean
 * curve.  The x-axis is the spatial coordinate (bin centres, in cm); values are
 * normalised per primary exactly as the ASCII writer does (NORM/SUM ÷ nstat,
 * AVER written as the physical mean).  The plot is an **additional** artifact:
 * BDO remains the source of truth and is never replaced.
 *
 * Anything outside that scope returns @c OSH_ENOTSUP for now (2-D / 3-D meshes,
 * differential spectra, ZONE geometry, rotated meshes, two-pass pages not yet
 * postprocessed).  These are follow-up items for the discussion on #238.
 *
 * @param[in] ws          Scoring workspace with output metadata and file paths.
 * @param[in] rt          Compiled scoring runtime with postprocessed accumulators.
 * @param[in] nstat       Actual number of primary particles simulated; must be > 0.
 * @param[in] output_idx  Zero-based index into ws->outputs / rt->outputs.
 *
 * @returns OSH_OK on success, OSH_ENOTSUP for an unsupported plot shape,
 *          OSH_EINVAL on a bad argument, or OSH_EIO on a write error.
 */
enum osh_status osh_scoring_save_plot_output(struct osh_scoring_workspace const *ws,
                                             struct osh_scoring_runtime const *rt,
                                             unsigned long long nstat,
                                             size_t output_idx);

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_SAVE_PLOT_H */
