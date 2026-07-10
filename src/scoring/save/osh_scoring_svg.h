#ifndef OSH_SCORING_SVG_H
#define OSH_SCORING_SVG_H

#include <stddef.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @file
 * @brief Generic, scoring-agnostic SVG line-plot renderer (issue #238).
 *
 * @details
 * Draws a framed 1-D line chart into an open @c FILE*: title/subtitle, a light
 * major + minor grid, "nice" axis ticks (linear or log10 x), one polyline per
 * series, and a right-hand legend.  Pure @c fprintf — no vendored library, no
 * font engine, no raster pipeline.
 *
 * This module knows nothing about the scoring model: callers pass plain arrays
 * of doubles and label strings.  The scoring-specific glue that turns pages and
 * geometry into those arrays lives in @c osh_scoring_save_plot.c.  Keeping the
 * two apart is what the #238 discussion anticipated for a future shared @c io/
 * renderer.
 */

/** @brief Fixed canvas size in SVG user units (pixels). */
#define OSH_SVG_W 760
#define OSH_SVG_H 480

/**
 * @brief Plot geometry: the pixel drawing area plus the "nice" data-space
 *        ranges that @ref osh_svg_plot_init derives from the raw data extents.
 *
 * Treat as opaque after init; the fields are public only so callers can
 * stack-allocate it.
 */
struct osh_svg_plot {
    double px0;   /* left pixel of the plot area */
    double px1;   /* right pixel */
    double py0;   /* top pixel */
    double py1;   /* bottom pixel */
    double xlo;   /* nice x data lower bound */
    double xhi;   /* nice x data upper bound */
    double ylo;   /* nice y data lower bound */
    double yhi;   /* nice y data upper bound */
    double xstep; /* linear x tick spacing (unused when xlog) */
    double ystep; /* linear y tick spacing */
    int xlog;     /* 0 = linear x, 1 = log10 x */
};

/**
 * @brief Derive the pixel plot area and nice axis bounds from raw data extents.
 *
 * @param[out] p     Plot to populate.
 * @param[in]  xmin  Smallest x value in the data.
 * @param[in]  xmax  Largest x value in the data.
 * @param[in]  xlog  Non-zero for a log10 x-axis (all x must be > 0).
 * @param[in]  ymin  Smallest y value to show (e.g. a 0 baseline).
 * @param[in]  ymax  Largest y value in the data.
 */
void osh_svg_plot_init(struct osh_svg_plot *p, double xmin, double xmax, int xlog, double ymin, double ymax);

/**
 * @brief Open the document and draw background, grid, ticks, border, and titles.
 *
 * Text arguments are XML-escaped internally; pass raw UTF-8.
 */
void osh_svg_begin(FILE *fp,
                   struct osh_svg_plot const *p,
                   char const *title,
                   char const *subtitle,
                   char const *xlabel,
                   char const *ylabel);

/**
 * @brief Draw one data series as a polyline (with point markers when @p npts is
 *        small).  @c ys[i] is plotted at @c xs[i]; the stroke colour is chosen
 *        from an internal palette by @p series_index.
 */
void osh_svg_series(
    FILE *fp, struct osh_svg_plot const *p, double const *xs, double const *ys, size_t npts, size_t series_index);

/**
 * @brief Draw one legend row (colour swatch + label) at slot @p series_index.
 *        @p label is XML-escaped internally.
 */
void osh_svg_legend_entry(FILE *fp, struct osh_svg_plot const *p, size_t series_index, char const *label);

/** @brief Close the document. */
void osh_svg_end(FILE *fp);

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_SVG_H */
