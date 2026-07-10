#ifndef OSH_SCORING_PLOT_SPEC_H
#define OSH_SCORING_PLOT_SPEC_H

#include <stddef.h>

#include "openshieldhit/status.h"
#include "scoring/runtime/osh_scoring_geometry_runtime.h"
#include "scoring/runtime/osh_scoring_output_runtime.h"
#include "scoring/runtime/osh_scoring_runtime.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @file
 * @brief Turn one compiled Output into a 1-D plot description (issue #238).
 *
 * @details
 * The scoring-model half of the SVG writer: it classifies an Output, selects the
 * pages that form a 1-D curve, and produces the shared x-coordinates and axis
 * labels.  The drawing is done elsewhere (osh_scoring_svg.h) from these plain
 * arrays; the file writing and normalisation live in osh_scoring_save_plot.c.
 */

/**
 * @brief A 1-D plot description built from one Output.
 *
 * @details
 * @c pages lists the page indices (into @c rt->pages) that are truly 1-D for the
 * chosen shape; each is drawn as one series at the shared @c xs coordinates.
 * Owns its @c xs and @c pages arrays — release with @ref osh_plot_spec_free.
 */
struct osh_plot_spec {
    size_t npts;     /* number of data points (== used page data length) */
    double *xs;      /* [npts] x-coordinates in data space (owned) */
    int xlog;        /* 1 = log10 x-axis */
    char xlabel[32]; /* x-axis title */
    char ylabel[96]; /* y-axis title */
    size_t *pages;   /* [nsel] plotted page indices into rt->pages (owned) */
    size_t nsel;     /* number of plotted pages */
};

/**
 * @brief Build a plot spec from one Output, selecting its plottable pages.
 *
 * @details
 * Tries a spatial profile first (MESH/CYL with one non-singleton axis, drawing
 * the plain pages), then a single-bin spectrum (a Diff1 page over a 0-D voxel or
 * single Zone, drawing the pages that share that differential axis).  Pages that
 * do not fit the chosen shape are skipped.  On failure @p spec is left safe to
 * pass to @ref osh_plot_spec_free.
 *
 * @returns OSH_OK, OSH_ENOTSUP when no page fits either shape, or OSH_ENOMEM.
 */
enum osh_status osh_plot_spec_build(struct osh_scoring_runtime const *rt,
                                    struct osh_scoring_output_runtime const *out,
                                    struct osh_scoring_geometry_runtime const *geo,
                                    struct osh_plot_spec *spec);

/** @brief Free the owned arrays and zero the pointers. */
void osh_plot_spec_free(struct osh_plot_spec *spec);

/**
 * @brief Per-primary scaling for a page's values: @p inv_nstat for NORM/SUM
 *        quantities, 1.0 for AVER (DLET/TLET) which already hold a physical mean.
 */
double osh_plot_page_scale(struct osh_scoring_page_runtime const *page, double inv_nstat);

/**
 * @brief Write the legend label (the uppercased quantity keyword) for the
 *        plotted series at index @p k into @p buf.
 */
void osh_plot_spec_series_label(
    struct osh_scoring_runtime const *rt, struct osh_plot_spec const *spec, size_t k, char *buf, size_t cap);

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_PLOT_SPEC_H */
