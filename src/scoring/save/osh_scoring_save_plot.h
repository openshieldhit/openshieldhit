#ifndef OSH_SCORING_SAVE_PLOT_H
#define OSH_SCORING_SAVE_PLOT_H

#include "openshieldhit/status.h"
#include "scoring/save/osh_scoring_save.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Save one compiled output as native, dependency-free SVG line plots.
 *
 * @details
 * Issue #238 — a zero-dependency "quick-look" plot so a run can drop a viewable
 * Bragg curve / 1-D profile on a machine that only has the `openshieldhit`
 * binary (no Python/matplotlib).  Always compiled in and reached by
 * `FileFormat SVG` on an `Output` block; the drawing is delegated to the
 * generic renderer in @ref osh_scoring_svg.h.  `FileFormat SVG` is that
 * `Output`'s only writer — it does not also emit `.bdo`/`.dat` for the same
 * block; add a second `Output` with `FileFormat BDO`/`TEXT` for the numeric
 * data.
 *
 * Two 1-D shapes are chosen automatically from the page model:
 *
 *   - **Spatial profile** — MESH with exactly one non-singleton X/Y/Z axis (e.g.
 *     a depth-dose / Bragg curve on a `1 x 1 x N` mesh), or CYL with one
 *     non-singleton R/Z axis.  x = spatial coordinate (bin centres, cm).
 *   - **Spectrum** — a differential page (dPhi/dE, dose-vs-LET, ...) scored over
 *     a single spatial bin: a `1 x 1 x 1` mesh voxel or a one-`Zone` geometry
 *     (both detected uniformly as `diff_stride == 1`).  x = the Diff1 bin
 *     centres of the differential quantity, log-scaled when the binning is LOG.
 *
 * An output may mix pages that fit the chosen shape with pages that do not — for
 * example a 1-D depth mesh whose extra page adds a Diff1 axis (making that page
 * 2-D spatial x energy).  Only the pages that are truly 1-D for the shape are
 * plotted; the rest are skipped.  Each plotted page is written to its own file
 * as a single black curve: one plotted quantity keeps the given filename, while
 * several each get a `_p1`/`_p2`/... page suffix before the `.svg` extension so
 * none is overwritten.  Values are normalised per primary exactly as the ASCII
 * writer does (NORM/SUM ÷ nstat, AVER written as the physical mean; a
 * differential page's bin-width division is already applied by @ref
 * osh_scoring_postprocess).
 *
 * @c OSH_ENOTSUP is returned when no page fits either shape (2-D / 3-D meshes,
 * 2-D spectra with a second differential axis, spectra spanning more than one
 * spatial bin, categorical multi-zone profiles, rotated meshes, two-pass pages
 * not yet postprocessed).  In that case the `Output` produces no file at all;
 * use `FileFormat TEXT`/`BDO` for that geometry and an external tool such as
 * pymchelper to visualise it.  These are follow-up items for the discussion on
 * #238.
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
