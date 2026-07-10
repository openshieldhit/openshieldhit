/*
 * Native SVG plot writer (issue #238) — orchestration.
 *
 * Turns one compiled Output into a viewable 1-D line plot next to its .bdo/.dat,
 * so a run can be eyeballed on a machine that has only the openshieldhit binary
 * (no Python/matplotlib).  This unit validates the request, asks
 * osh_scoring_plot_spec.c which pages form a 1-D curve, computes the data
 * ranges, and drives the generic renderer in osh_scoring_svg.c to a file.  It is
 * compiled unconditionally; `FileFormat SVG` on an Output block reaches it
 * through the save dispatcher.
 *
 * Normalisation matches the ASCII writer: NORM/SUM quantities are divided by
 * nstat, AVER quantities (DLET/TLET) are written as the physical mean already
 * produced by osh_scoring_postprocess() — which is also where a differential
 * page's bin-width division happens, so acc.data already holds dPhi/dE here.
 * The plot is an additional artifact — BDO stays the source of truth.  Save is a
 * cold path, so the small allocations here are fine (the §10 hot-path ban does
 * not apply).
 */

#include "scoring/save/osh_scoring_save_plot.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "scoring/runtime/osh_scoring_geometry_runtime.h"
#include "scoring/runtime/osh_scoring_output_runtime.h"
#include "scoring/runtime/osh_scoring_runtime.h"
#include "scoring/save/osh_scoring_plot_spec.h"
#include "scoring/save/osh_scoring_svg.h"

static enum osh_status validate_plot_output(struct osh_scoring_workspace const *ws,
                                            struct osh_scoring_runtime const *rt,
                                            size_t output_idx,
                                            struct osh_scoring_output_runtime const **out_out,
                                            struct osh_scoring_geometry_runtime const **geo_out);
static void compute_ranges(struct osh_plot_spec const *spec,
                           struct osh_scoring_runtime const *rt,
                           double inv_nstat,
                           double *xmin_out,
                           double *xmax_out,
                           double *ylo_out,
                           double *yhi_out);
static void render(FILE *fp,
                   struct osh_svg_plot const *plot,
                   struct osh_plot_spec const *spec,
                   struct osh_scoring_runtime const *rt,
                   double inv_nstat,
                   double *ys,
                   char const *title,
                   char const *subtitle);
static char *plot_output_path(char const *filename);
static char const *path_basename(char const *path);

enum osh_status osh_scoring_save_plot_output(struct osh_scoring_workspace const *ws,
                                             struct osh_scoring_runtime const *rt,
                                             unsigned long long nstat,
                                             size_t output_idx) {
    struct osh_scoring_output_runtime const *out;
    struct osh_scoring_geometry_runtime const *geo;
    struct osh_plot_spec spec;
    struct osh_svg_plot plot;
    FILE *fp;
    char *out_path;
    double *ys;
    double inv_nstat;
    double xmin;
    double xmax;
    double ylo;
    double yhi;
    char subtitle[128];
    enum osh_status rc;

    rc = validate_plot_output(ws, rt, output_idx, &out, &geo);
    if (rc != OSH_OK) {
        return rc;
    }
    if (nstat == 0ull) {
        return OSH_EINVAL;
    }

    rc = osh_plot_spec_build(rt, out, geo, &spec);
    if (rc != OSH_OK) {
        osh_plot_spec_free(&spec);
        return rc;
    }
    inv_nstat = 1.0 / (double) nstat;

    ys = (double *) malloc(spec.npts * sizeof(*ys));
    if (!ys) {
        osh_plot_spec_free(&spec);
        return OSH_ENOMEM;
    }
    compute_ranges(&spec, rt, inv_nstat, &xmin, &xmax, &ylo, &yhi);
    osh_svg_plot_init(&plot, xmin, xmax, spec.xlog, ylo, yhi);
    snprintf(
        subtitle, sizeof(subtitle), "geometry: %s  \xe2\x80\x94  %llu primaries", geo->name ? geo->name : "?", nstat);

    out_path = plot_output_path(out->filename);
    if (!out_path) {
        free(ys);
        osh_plot_spec_free(&spec);
        return OSH_ENOMEM;
    }
    fp = fopen(out_path, "w");
    if (!fp) {
        free(out_path);
        free(ys);
        osh_plot_spec_free(&spec);
        return OSH_EIO;
    }

    render(fp, &plot, &spec, rt, inv_nstat, ys, path_basename(out_path), subtitle);

    free(out_path);
    free(ys);
    osh_plot_spec_free(&spec);
    if (fclose(fp) != 0) {
        return OSH_EIO;
    }
    return OSH_OK;
}

static enum osh_status validate_plot_output(struct osh_scoring_workspace const *ws,
                                            struct osh_scoring_runtime const *rt,
                                            size_t output_idx,
                                            struct osh_scoring_output_runtime const **out_out,
                                            struct osh_scoring_geometry_runtime const **geo_out) {
    struct osh_scoring_output_runtime const *out;
    struct osh_scoring_geometry_runtime const *geo;
    size_t ip;

    if (!ws || !rt || !out_out || !geo_out) {
        return OSH_EINVAL;
    }
    if (output_idx >= rt->noutputs || output_idx >= ws->noutputs) {
        return OSH_EINVAL;
    }

    out = &rt->outputs[output_idx];
    if (out->geometry_idx >= rt->ngeometries) {
        return OSH_ESTATE;
    }
    geo = &rt->geometries[out->geometry_idx];
    if (out->npages == 0u || !out->page_indices) {
        return OSH_ENOTSUP; /* nothing to draw */
    }
    for (ip = 0; ip < out->npages; ++ip) {
        size_t page_idx = out->page_indices[ip];
        struct osh_scoring_page_runtime const *page;

        if (page_idx >= rt->npages) {
            return OSH_ESTATE;
        }
        page = &rt->pages[page_idx];
        if (page->geometry_idx != out->geometry_idx) {
            return OSH_ESTATE;
        }
        /* Two-pass pages must be finalised by postprocess before save. */
        if (!page->acc.data || page->has_data2 || page->divide) {
            return OSH_ENOTSUP;
        }
    }

    *out_out = out;
    *geo_out = geo;
    return OSH_OK;
}

/* x range from the shared coordinate array; y range across every plotted page.
 * Point i indexes acc.data[i] for both shapes (see osh_scoring_plot_spec.c).
 * Positive-only data gets a natural 0 baseline; data that dips negative keeps
 * its true minimum. */
static void compute_ranges(struct osh_plot_spec const *spec,
                           struct osh_scoring_runtime const *rt,
                           double inv_nstat,
                           double *xmin_out,
                           double *xmax_out,
                           double *ylo_out,
                           double *yhi_out) {
    double xmin = spec->xs[0];
    double xmax = spec->xs[0];
    double ymin = HUGE_VAL;
    double ymax = -HUGE_VAL;
    double ylo;
    size_t k;
    size_t i;

    for (i = 0; i < spec->npts; ++i) {
        if (spec->xs[i] < xmin) {
            xmin = spec->xs[i];
        }
        if (spec->xs[i] > xmax) {
            xmax = spec->xs[i];
        }
    }
    for (k = 0; k < spec->nsel; ++k) {
        struct osh_scoring_page_runtime const *page = &rt->pages[spec->pages[k]];
        double scale = osh_plot_page_scale(page, inv_nstat);
        for (i = 0; i < spec->npts; ++i) {
            double v = page->acc.data[i] * scale;
            if (v < ymin) {
                ymin = v;
            }
            if (v > ymax) {
                ymax = v;
            }
        }
    }
    if (ymin > ymax) {
        ymin = 0.0;
        ymax = 1.0;
    }
    ylo = (ymin >= 0.0) ? 0.0 : ymin;
    if (ymax <= ylo) {
        ymax = ylo + 1.0;
    }
    *xmin_out = xmin;
    *xmax_out = xmax;
    *ylo_out = ylo;
    *yhi_out = ymax;
}

/* Draw the document: frame, one series per plotted page (into the caller's @p ys
 * scratch), then the legend. */
static void render(FILE *fp,
                   struct osh_svg_plot const *plot,
                   struct osh_plot_spec const *spec,
                   struct osh_scoring_runtime const *rt,
                   double inv_nstat,
                   double *ys,
                   char const *title,
                   char const *subtitle) {
    size_t k;
    size_t i;

    osh_svg_begin(fp, plot, title, subtitle, spec->xlabel, spec->ylabel);
    for (k = 0; k < spec->nsel; ++k) {
        struct osh_scoring_page_runtime const *page = &rt->pages[spec->pages[k]];
        double scale = osh_plot_page_scale(page, inv_nstat);
        for (i = 0; i < spec->npts; ++i) {
            ys[i] = page->acc.data[i] * scale;
        }
        osh_svg_series(fp, plot, spec->xs, ys, spec->npts, k);
    }
    for (k = 0; k < spec->nsel; ++k) {
        char label[64];
        osh_plot_spec_series_label(rt, spec, k, label, sizeof(label));
        osh_svg_legend_entry(fp, plot, k, label);
    }
    osh_svg_end(fp);
}

/* Return a heap copy of @p filename with ".svg" appended unless already present,
 * so the artifact opens in a browser by double-click.  Mirrors the RTDOSE
 * writer's ".dcm" handling.  Caller frees. */
static char *plot_output_path(char const *filename) {
    size_t len;
    char *result;

    if (!filename) {
        return NULL;
    }
    len = strlen(filename);
    if (len >= 4u && strcmp(filename + len - 4u, ".svg") == 0) {
        result = (char *) malloc(len + 1u);
        if (result) {
            memcpy(result, filename, len + 1u);
        }
        return result;
    }
    result = (char *) malloc(len + 5u);
    if (result) {
        memcpy(result, filename, len);
        memcpy(result + len, ".svg", 5u);
    }
    return result;
}

/* Pointer to the last path component of @p path (portable over '/' and '\\'). */
static char const *path_basename(char const *path) {
    char const *base;
    char const *c;

    base = path;
    for (c = path; *c; ++c) {
        if (*c == '/' || *c == '\\') {
            base = c + 1;
        }
    }
    return base;
}
