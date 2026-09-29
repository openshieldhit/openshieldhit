/*
 * Native SVG plot writer (issue #238) — orchestration.
 *
 * Turns one compiled Output into viewable 1-D line plots, so a run can be
 * eyeballed on a machine that has only the openshieldhit binary (no
 * Python/matplotlib).  `FileFormat SVG` is that Output's only writer — it does
 * not also emit `.bdo`/`.dat` for the same block; a separate Output with
 * `FileFormat BDO`/`TEXT` covers the numeric data.  This unit validates the
 * request, asks osh_scoring_plot_spec.c which pages form a 1-D curve, and
 * writes one SVG file per plotted page via the generic renderer in
 * osh_scoring_svg.c.  A single plotted quantity keeps the given filename;
 * several quantities each get their own file with a "_p1"/"_p2"/... page
 * suffix so none is overwritten.  It is compiled unconditionally; `FileFormat
 * SVG` on an Output block reaches it through the save dispatcher.
 *
 * Normalisation matches the ASCII writer: NORM/SUM quantities are divided by
 * nstat, AVER quantities (DLET/TLET) are written as the physical mean already
 * produced by osh_scoring_postprocess() — which is also where a differential
 * page's bin-width division happens, so acc.data already holds dPhi/dE here.
 * Save is a cold path, so the small allocations here are fine (the §10
 * hot-path ban does not apply).
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
static void compute_x_range(struct osh_plot_spec const *spec, double *xmin_out, double *xmax_out);
static void compute_page_y_range(
    struct osh_scoring_page_runtime const *page, size_t npts, double scale, double *ylo_out, double *yhi_out);
static enum osh_status save_one_page(char const *filename,
                                     struct osh_plot_spec const *spec,
                                     struct osh_scoring_runtime const *rt,
                                     size_t k,
                                     double inv_nstat,
                                     double xmin,
                                     double xmax,
                                     char const *subtitle,
                                     double *ys);
static char *plot_output_path(char const *filename, size_t page_no, size_t npages);
static char const *path_basename(char const *path);

enum osh_status osh_scoring_save_plot_output(struct osh_scoring_workspace const *ws,
                                             struct osh_scoring_runtime const *rt,
                                             unsigned long long nstat,
                                             size_t output_idx) {
    struct osh_scoring_output_runtime const *out;
    struct osh_scoring_geometry_runtime const *geo;
    struct osh_plot_spec spec;
    double *ys;
    double inv_nstat;
    double xmin;
    double xmax;
    char subtitle[128];
    char const *geo_name;
    size_t k;
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
    compute_x_range(&spec, &xmin, &xmax);
    geo_name = geo->name;
    if (!geo_name) {
        geo_name = "?";
    }
    snprintf(subtitle, sizeof(subtitle), "geometry: %s  \xe2\x80\x94  %llu primaries", geo_name, nstat);

    /* One file per plotted page (see plot_output_path for the _pN naming). */
    rc = OSH_OK;
    for (k = 0; k < spec.nsel; ++k) {
        rc = save_one_page(out->filename, &spec, rt, k, inv_nstat, xmin, xmax, subtitle, ys);
        if (rc != OSH_OK) {
            break;
        }
    }

    free(ys);
    osh_plot_spec_free(&spec);
    return rc;
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
    /* Expanded multi-format targets live past ws->noutputs (issue #308); every
     * field this writer needs comes from rt->outputs[output_idx]. */
    if (output_idx >= rt->noutputs) {
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

/* x range from the shared coordinate array (identical for every plotted page,
 * see osh_scoring_plot_spec.c). */
static void compute_x_range(struct osh_plot_spec const *spec, double *xmin_out, double *xmax_out) {
    double xmin = spec->xs[0];
    double xmax = spec->xs[0];
    size_t i;

    for (i = 0; i < spec->npts; ++i) {
        if (spec->xs[i] < xmin) {
            xmin = spec->xs[i];
        }
        if (spec->xs[i] > xmax) {
            xmax = spec->xs[i];
        }
    }
    *xmin_out = xmin;
    *xmax_out = xmax;
}

/* y range for one page's scaled values.  Point i indexes page->acc.data[i].
 * Positive-only data gets a natural 0 baseline; data that dips negative keeps
 * its true minimum. */
static void compute_page_y_range(
    struct osh_scoring_page_runtime const *page, size_t npts, double scale, double *ylo_out, double *yhi_out) {
    double ymin = HUGE_VAL;
    double ymax = -HUGE_VAL;
    double ylo;
    size_t i;

    for (i = 0; i < npts; ++i) {
        double v = page->acc.data[i] * scale;
        if (v < ymin) {
            ymin = v;
        }
        if (v > ymax) {
            ymax = v;
        }
    }
    if (ymin > ymax) {
        ymin = 0.0;
        ymax = 1.0;
    }
    if (ymin >= 0.0) {
        ylo = 0.0;
    } else {
        ylo = ymin;
    }
    if (ymax <= ylo) {
        ymax = ylo + 1.0;
    }
    *ylo_out = ylo;
    *yhi_out = ymax;
}

/* Draw plotted page @p k to its own SVG file: derive its y-range and title,
 * open the file, and emit the frame plus the single (black) data curve into the
 * caller's @p ys scratch. */
static enum osh_status save_one_page(char const *filename,
                                     struct osh_plot_spec const *spec,
                                     struct osh_scoring_runtime const *rt,
                                     size_t k,
                                     double inv_nstat,
                                     double xmin,
                                     double xmax,
                                     char const *subtitle,
                                     double *ys) {
    struct osh_scoring_page_runtime const *page = &rt->pages[spec->pages[k]];
    double scale = osh_plot_page_scale(page, inv_nstat);
    struct osh_svg_plot plot;
    char ylabel[96];
    char *out_path;
    FILE *fp;
    double ylo;
    double yhi;
    size_t i;

    compute_page_y_range(page, spec->npts, scale, &ylo, &yhi);
    osh_svg_plot_init(&plot, xmin, xmax, spec->xlog, ylo, yhi);

    out_path = plot_output_path(filename, k, spec->nsel);
    if (!out_path) {
        return OSH_ENOMEM;
    }
    fp = fopen(out_path, "w");
    if (!fp) {
        free(out_path);
        return OSH_EIO;
    }

    osh_plot_spec_page_ylabel(rt, spec, k, ylabel, sizeof(ylabel));
    osh_svg_begin(fp, &plot, path_basename(out_path), subtitle, spec->xlabel, ylabel);
    for (i = 0; i < spec->npts; ++i) {
        ys[i] = page->acc.data[i] * scale;
    }
    osh_svg_series(fp, &plot, spec->xs, ys, spec->npts);
    osh_svg_end(fp);

    free(out_path);
    if (fclose(fp) != 0) {
        return OSH_EIO;
    }
    return OSH_OK;
}

/* Build the output path for plotted page @p page_no of @p npages.  ".svg" is
 * appended unless already present, and — only when several pages are plotted —
 * a "_p<n>" suffix is inserted before the extension so the files don't collide
 * (bragg.svg -> bragg_p1.svg, bragg_p2.svg, ...).  A lone page keeps the given
 * name.  Mirrors the RTDOSE writer's ".dcm" handling.  Caller frees. */
static char *plot_output_path(char const *filename, size_t page_no, size_t npages) {
    size_t stem_len;
    size_t cap;
    char *result;

    if (!filename) {
        return NULL;
    }
    stem_len = strlen(filename);
    if (stem_len >= 4u && strcmp(filename + stem_len - 4u, ".svg") == 0) {
        stem_len -= 4u; /* drop the extension; re-appended below */
    }
    /* stem + "_p<number>" + ".svg" + NUL; 32 bytes covers the suffix, the
     * extension, and any size_t rendered by %zu. */
    cap = stem_len + 32u;
    result = (char *) malloc(cap);
    if (!result) {
        return NULL;
    }
    memcpy(result, filename, stem_len);
    if (npages > 1u) {
        snprintf(result + stem_len, cap - stem_len, "_p%zu.svg", page_no + 1u);
    } else {
        snprintf(result + stem_len, cap - stem_len, ".svg");
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
