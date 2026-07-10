/*
 * Native SVG plot writer (prototype, issue #238).
 *
 * A zero-dependency "quick-look" renderer: pure fprintf into an SVG document, no
 * vendored library, no font engine, no raster pipeline — just a plot frame,
 * "nice" axis ticks, and one <polyline> per page.  It exists so a run can drop a
 * viewable 1-D profile (a depth-dose / Bragg curve) next to its .bdo/.dat on a
 * machine that has only the openshieldhit binary.  The whole translation unit is
 * compiled in only under -DOSH_ENABLE_PLOT=ON, so the stock binary is unchanged.
 *
 * Normalisation matches the ASCII writer: NORM/SUM quantities are divided by
 * nstat at write time, AVER quantities (DLET/TLET) are written as the physical
 * mean already produced by osh_scoring_postprocess().  The plot is an additional
 * artifact — BDO stays the source of truth.  Save is a cold path, so the §10
 * hot-path allocation ban does not apply (this writer allocates nothing anyway).
 *
 * Scope is deliberately narrow for the discussion prototype: 1-D spatial
 * profiles on MESH (one non-singleton X/Y/Z axis) or CYL (one non-singleton
 * R/Z axis).  Everything else returns OSH_ENOTSUP; 2-D heatmaps, differential
 * spectra, multipage PDF and PNG are follow-up items on #238.
 */

#include "scoring/save/osh_scoring_save_plot.h"

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "scoring/runtime/osh_scoring_geometry_runtime.h"
#include "scoring/runtime/osh_scoring_output_runtime.h"
#include "scoring/runtime/osh_scoring_runtime.h"

/* Canvas geometry, in SVG user units (pixels). */
#define PLOT_W 760
#define PLOT_H 480
#define PLOT_MARGIN_L 74
#define PLOT_MARGIN_R 166 /* right margin leaves room for the legend column */
#define PLOT_MARGIN_T 48
#define PLOT_MARGIN_B 58
#define PLOT_NTICK 6       /* target tick count per axis (Heckbert picks the exact count) */
#define PLOT_MARKER_MAX 40 /* draw per-point markers only up to this many bins */

/* Resolved pixel plot area plus the "nice" data-space ranges it maps. */
struct plot_layout {
    double px0; /* left pixel of the plot area */
    double px1; /* right pixel */
    double py0; /* top pixel */
    double py1; /* bottom pixel */
    double xlo; /* data-space x range (nice bounds) */
    double xhi;
    double ylo; /* data-space y range (nice bounds) */
    double yhi;
    double xstep; /* tick spacing in data space */
    double ystep;
};

static enum osh_status validate_plot_output(struct osh_scoring_workspace const *ws,
                                            struct osh_scoring_runtime const *rt,
                                            size_t output_idx,
                                            struct osh_scoring_output_runtime const **out_out,
                                            struct osh_scoring_geometry_runtime const **geo_out);
static enum osh_status pick_profile_axis(struct osh_scoring_geometry_runtime const *geo, size_t *axis_idx_out);
static double page_scale(struct osh_scoring_page_runtime const *page, double inv_nstat);
static char *plot_output_path(char const *filename);
static char const *path_basename(char const *path);
static char const *plot_data_unit(struct osh_scoring_page_runtime const *page);
static char const *series_color(size_t i);
static double nice_num(double range, int do_round);
static void nice_axis(double lo, double hi, int ntick, double *nlo, double *nhi, double *step);
static void upcase_copy(char *dst, size_t cap, char const *src);
static void svg_escape(FILE *fp, char const *s);
static void svg_write_frame(FILE *fp,
                            struct plot_layout const *lay,
                            char const *title,
                            char const *subtitle,
                            char const *xlabel,
                            char const *ylabel);
static void svg_write_series(FILE *fp,
                             struct plot_layout const *lay,
                             double x0,
                             double dx,
                             size_t nbins,
                             double const *data,
                             double scale,
                             char const *color);
static void svg_write_legend(FILE *fp,
                             struct plot_layout const *lay,
                             struct osh_scoring_runtime const *rt,
                             struct osh_scoring_output_runtime const *out);

static inline double _map_x(struct plot_layout const *lay, double x) {
    return lay->px0 + (x - lay->xlo) / (lay->xhi - lay->xlo) * (lay->px1 - lay->px0);
}

static inline double _map_y(struct plot_layout const *lay, double y) {
    /* Larger data values map to smaller pixel-y (top of the canvas). */
    return lay->py1 - (y - lay->ylo) / (lay->yhi - lay->ylo) * (lay->py1 - lay->py0);
}

enum osh_status osh_scoring_save_plot_output(struct osh_scoring_workspace const *ws,
                                             struct osh_scoring_runtime const *rt,
                                             unsigned long long nstat,
                                             size_t output_idx) {
    struct osh_scoring_output_runtime const *out;
    struct osh_scoring_geometry_runtime const *geo;
    struct osh_scoring_axis_runtime const *axis;
    struct osh_scoring_page_runtime const *p0;
    struct plot_layout lay;
    FILE *fp;
    char *out_path;
    double inv_nstat;
    double ymin;
    double ymax;
    double ylo_data;
    double x0;
    double dx;
    size_t axis_idx;
    size_t nbins;
    size_t ip;
    size_t i;
    char qty_upper[64];
    char ylabel[96];
    char xlabel[16];
    char subtitle[128];
    enum osh_status rc;

    rc = validate_plot_output(ws, rt, output_idx, &out, &geo);
    if (rc != OSH_OK) {
        return rc;
    }
    if (nstat == 0ull) {
        return OSH_EINVAL;
    }

    rc = pick_profile_axis(geo, &axis_idx);
    if (rc != OSH_OK) {
        return rc;
    }
    axis = &geo->axes[axis_idx];
    nbins = (size_t) axis->nbins;
    x0 = axis->lo;
    dx = (axis->hi - axis->lo) / (double) nbins;
    inv_nstat = 1.0 / (double) nstat;

    /* y-range across every page.  Since only one spatial axis is non-singleton,
     * the canonical flat index equals the profile bin index, so data[i] is the
     * value at profile bin i for every page (see pick_profile_axis). */
    ymin = HUGE_VAL;
    ymax = -HUGE_VAL;
    for (ip = 0; ip < out->npages; ++ip) {
        struct osh_scoring_page_runtime const *page = &rt->pages[out->page_indices[ip]];
        double scale = page_scale(page, inv_nstat);
        for (i = 0; i < nbins; ++i) {
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
    /* Positive-only data gets a natural 0 baseline (dose/fluence curves); data
     * that dips negative keeps its true minimum. */
    ylo_data = (ymin >= 0.0) ? 0.0 : ymin;
    if (ymax <= ylo_data) {
        ymax = ylo_data + 1.0;
    }

    lay.px0 = (double) PLOT_MARGIN_L;
    lay.px1 = (double) (PLOT_W - PLOT_MARGIN_R);
    lay.py0 = (double) PLOT_MARGIN_T;
    lay.py1 = (double) (PLOT_H - PLOT_MARGIN_B);
    nice_axis(axis->lo, axis->hi, PLOT_NTICK, &lay.xlo, &lay.xhi, &lay.xstep);
    nice_axis(ylo_data, ymax, PLOT_NTICK, &lay.ylo, &lay.yhi, &lay.ystep);

    p0 = &rt->pages[out->page_indices[0]];
    upcase_copy(qty_upper, sizeof(qty_upper), p0->quantity ? p0->quantity : "value");
    if (out->npages == 1u) {
        snprintf(ylabel, sizeof(ylabel), "%s [%s]", qty_upper, plot_data_unit(p0));
    } else {
        snprintf(ylabel, sizeof(ylabel), "value [%s]", plot_data_unit(p0));
    }
    snprintf(xlabel, sizeof(xlabel), "%s [cm]", axis->label);
    snprintf(
        subtitle, sizeof(subtitle), "geometry: %s  \xe2\x80\x94  %llu primaries", geo->name ? geo->name : "?", nstat);

    out_path = plot_output_path(out->filename);
    if (!out_path) {
        return OSH_ENOMEM;
    }
    fp = fopen(out_path, "w");
    if (!fp) {
        free(out_path);
        return OSH_EIO;
    }

    svg_write_frame(fp, &lay, path_basename(out_path), subtitle, xlabel, ylabel);
    for (ip = 0; ip < out->npages; ++ip) {
        struct osh_scoring_page_runtime const *page = &rt->pages[out->page_indices[ip]];
        svg_write_series(fp, &lay, x0, dx, nbins, page->acc.data, page_scale(page, inv_nstat), series_color(ip));
    }
    svg_write_legend(fp, &lay, rt, out);
    fputs("</svg>\n", fp);

    free(out_path);
    if (fclose(fp) != 0) {
        return OSH_EIO;
    }
    return OSH_OK;
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

/* NORM/SUM quantities are divided by nstat; AVER quantities (DLET/TLET) already
 * hold a physical mean after postprocess and must not be. */
static double page_scale(struct osh_scoring_page_runtime const *page, double inv_nstat) {
    return (page->postproc == OSH_SCORING_POSTPROC_AVER) ? 1.0 : inv_nstat;
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
    if (geo->geo_kind != OSH_SCORING_GEO_MESH && geo->geo_kind != OSH_SCORING_GEO_CYL) {
        return OSH_ENOTSUP; /* ZONE / voxel line plots are deferred */
    }
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
        /* Differential spectra (dPhi/dE, dose-vs-LET) are a separate 1-D case
         * that this prototype does not draw yet — see #238 / #215. */
        if (page->diff_nbins > 0u) {
            return OSH_ENOTSUP;
        }
    }

    *out_out = out;
    *geo_out = geo;
    return OSH_OK;
}

/* Find the single non-singleton spatial axis (the profile axis).  Returns
 * OSH_ENOTSUP unless exactly one axis has nbins > 1: zero means a single point
 * (nothing to plot as a curve), more than one means a 2-D/3-D map (deferred).
 * A rotated MESH is rejected because its X/Y/Z are local-frame and ambiguous,
 * mirroring the ASCII writer; CYL R/Z is always local-frame and is allowed. */
static enum osh_status pick_profile_axis(struct osh_scoring_geometry_runtime const *geo, size_t *axis_idx_out) {
    size_t i;
    size_t found;
    size_t count;

    if (geo->geo_kind == OSH_SCORING_GEO_MESH && geo->has_rotation) {
        return OSH_ENOTSUP;
    }
    found = 0u;
    count = 0u;
    for (i = 0; i < geo->naxes; ++i) {
        if (geo->axes[i].nbins > 1) {
            found = i;
            ++count;
        }
    }
    if (count != 1u) {
        return OSH_ENOTSUP;
    }
    *axis_idx_out = found;
    return OSH_OK;
}

/* Base unit of the scored quantity, mirroring page_data_unit() in the BDO
 * writer.  Differential axes are rejected upstream, so no denominator handling. */
static char const *plot_data_unit(struct osh_scoring_page_runtime const *page) {
    switch (page->score_kind) {
    case OSH_SCORING_SCORE_ENERGY:
        return "MeV";
    case OSH_SCORING_SCORE_FLUENCE:
        return "/cm^2";
    case OSH_SCORING_SCORE_DOSE:
    case OSH_SCORING_SCORE_DIRTYDOSE:
        return "MeV/g";
    case OSH_SCORING_SCORE_DOSEGY:
    case OSH_SCORING_SCORE_DIRTYDOSEGY:
        return "Gy";
    case OSH_SCORING_SCORE_DLET:
    case OSH_SCORING_SCORE_TLET:
        return "MeV/cm";
    case OSH_SCORING_SCORE_DQEFF:
    case OSH_SCORING_SCORE_TQEFF:
        return "dim.less";
    default:
        return "arb";
    }
}

/* A small qualitative palette (matplotlib "tab10" leaders), cycled per page. */
static char const *series_color(size_t i) {
    static char const *const palette[] = {
        "#1f77b4",
        "#d62728",
        "#2ca02c",
        "#ff7f0e",
        "#9467bd",
        "#8c564b",
        "#17becf",
        "#bcbd22",
    };
    return palette[i % (sizeof(palette) / sizeof(palette[0]))];
}

/* Heckbert's "nice numbers for graph labels": round @p range to 1/2/5 x 10^k.
 * @p do_round selects nearest (1) vs. ceiling (0). */
static double nice_num(double range, int do_round) {
    double exponent;
    double fraction;
    double nice;

    if (!(range > 0.0)) {
        return 1.0;
    }
    exponent = floor(log10(range));
    fraction = range / pow(10.0, exponent);
    if (do_round) {
        if (fraction < 1.5) {
            nice = 1.0;
        } else if (fraction < 3.0) {
            nice = 2.0;
        } else if (fraction < 7.0) {
            nice = 5.0;
        } else {
            nice = 10.0;
        }
    } else {
        if (fraction <= 1.0) {
            nice = 1.0;
        } else if (fraction <= 2.0) {
            nice = 2.0;
        } else if (fraction <= 5.0) {
            nice = 5.0;
        } else {
            nice = 10.0;
        }
    }
    return nice * pow(10.0, exponent);
}

/* Expand [lo, hi] to nice round bounds and pick a nice tick step. */
static void nice_axis(double lo, double hi, int ntick, double *nlo, double *nhi, double *step) {
    double range;
    double d;

    if (hi <= lo) {
        hi = lo + 1.0;
    }
    range = nice_num(hi - lo, 0);
    d = nice_num(range / (double) (ntick - 1), 1);
    *nlo = floor(lo / d) * d;
    *nhi = ceil(hi / d) * d;
    *step = d;
}

static void upcase_copy(char *dst, size_t cap, char const *src) {
    size_t i;

    if (cap == 0u) {
        return;
    }
    for (i = 0; i + 1u < cap && src[i]; ++i) {
        dst[i] = (char) toupper((unsigned char) src[i]);
    }
    dst[i] = '\0';
}

/* Emit @p s as SVG text content, escaping the XML metacharacters. */
static void svg_escape(FILE *fp, char const *s) {
    char const *c;

    if (!s) {
        return;
    }
    for (c = s; *c; ++c) {
        switch (*c) {
        case '&':
            fputs("&amp;", fp);
            break;
        case '<':
            fputs("&lt;", fp);
            break;
        case '>':
            fputs("&gt;", fp);
            break;
        default:
            fputc(*c, fp);
            break;
        }
    }
}

static void svg_write_frame(FILE *fp,
                            struct plot_layout const *lay,
                            char const *title,
                            char const *subtitle,
                            char const *xlabel,
                            char const *ylabel) {
    int nx;
    int ny;
    int k;

    fprintf(fp, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    fprintf(fp,
            "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%d\" height=\"%d\" "
            "viewBox=\"0 0 %d %d\" font-family=\"sans-serif\">\n",
            PLOT_W,
            PLOT_H,
            PLOT_W,
            PLOT_H);
    fprintf(fp, "<rect width=\"%d\" height=\"%d\" fill=\"#ffffff\"/>\n", PLOT_W, PLOT_H);

    /* Title + subtitle. */
    fprintf(fp,
            "<text x=\"%.1f\" y=\"22\" text-anchor=\"middle\" font-size=\"16\" font-weight=\"bold\">",
            0.5 * (lay->px0 + lay->px1));
    svg_escape(fp, title);
    fputs("</text>\n", fp);
    fprintf(fp,
            "<text x=\"%.1f\" y=\"39\" text-anchor=\"middle\" font-size=\"11\" fill=\"#555555\">",
            0.5 * (lay->px0 + lay->px1));
    svg_escape(fp, subtitle);
    fputs("</text>\n", fp);

    nx = (int) floor((lay->xhi - lay->xlo) / lay->xstep + 0.5) + 1;
    ny = (int) floor((lay->yhi - lay->ylo) / lay->ystep + 0.5) + 1;

    /* Gridlines + tick labels. */
    for (k = 0; k < nx; ++k) {
        double v = lay->xlo + (double) k * lay->xstep;
        double px = _map_x(lay, v);
        fprintf(fp,
                "<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"#e6e6e6\"/>\n",
                px,
                lay->py0,
                px,
                lay->py1);
        fprintf(fp,
                "<text x=\"%.1f\" y=\"%.1f\" text-anchor=\"middle\" font-size=\"11\">%g</text>\n",
                px,
                lay->py1 + 16.0,
                v);
    }
    for (k = 0; k < ny; ++k) {
        double v = lay->ylo + (double) k * lay->ystep;
        double py = _map_y(lay, v);
        fprintf(fp,
                "<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"#e6e6e6\"/>\n",
                lay->px0,
                py,
                lay->px1,
                py);
        fprintf(fp,
                "<text x=\"%.1f\" y=\"%.1f\" text-anchor=\"end\" font-size=\"11\">%g</text>\n",
                lay->px0 - 8.0,
                py + 4.0,
                v);
    }

    /* Plot border. */
    fprintf(fp,
            "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"none\" stroke=\"#333333\"/>\n",
            lay->px0,
            lay->py0,
            lay->px1 - lay->px0,
            lay->py1 - lay->py0);

    /* Axis titles. */
    fprintf(fp,
            "<text x=\"%.1f\" y=\"%d\" text-anchor=\"middle\" font-size=\"13\">",
            0.5 * (lay->px0 + lay->px1),
            PLOT_H - 16);
    svg_escape(fp, xlabel);
    fputs("</text>\n", fp);
    fprintf(fp,
            "<text x=\"18\" y=\"%.1f\" text-anchor=\"middle\" font-size=\"13\" "
            "transform=\"rotate(-90 18 %.1f)\">",
            0.5 * (lay->py0 + lay->py1),
            0.5 * (lay->py0 + lay->py1));
    svg_escape(fp, ylabel);
    fputs("</text>\n", fp);
}

static void svg_write_series(FILE *fp,
                             struct plot_layout const *lay,
                             double x0,
                             double dx,
                             size_t nbins,
                             double const *data,
                             double scale,
                             char const *color) {
    size_t i;

    fprintf(fp, "<polyline fill=\"none\" stroke=\"%s\" stroke-width=\"2\" points=\"", color);
    for (i = 0; i < nbins; ++i) {
        double xc = x0 + ((double) i + 0.5) * dx;
        double v = data[i] * scale;
        fprintf(fp, "%.2f,%.2f ", _map_x(lay, xc), _map_y(lay, v));
    }
    fputs("\"/>\n", fp);

    if (nbins <= PLOT_MARKER_MAX) {
        for (i = 0; i < nbins; ++i) {
            double xc = x0 + ((double) i + 0.5) * dx;
            double v = data[i] * scale;
            fprintf(fp,
                    "<circle cx=\"%.2f\" cy=\"%.2f\" r=\"2.5\" fill=\"%s\"/>\n",
                    _map_x(lay, xc),
                    _map_y(lay, v),
                    color);
        }
    }
}

static void svg_write_legend(FILE *fp,
                             struct plot_layout const *lay,
                             struct osh_scoring_runtime const *rt,
                             struct osh_scoring_output_runtime const *out) {
    double lx;
    double ly;
    size_t ip;

    lx = lay->px1 + 16.0;
    ly = lay->py0 + 6.0;
    for (ip = 0; ip < out->npages; ++ip) {
        struct osh_scoring_page_runtime const *page = &rt->pages[out->page_indices[ip]];
        char qty_upper[64];
        double row = ly + (double) ip * 20.0;

        upcase_copy(qty_upper, sizeof(qty_upper), page->quantity ? page->quantity : "?");
        fprintf(fp,
                "<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"%s\" stroke-width=\"3\"/>\n",
                lx,
                row,
                lx + 22.0,
                row,
                series_color(ip));
        fprintf(fp, "<text x=\"%.1f\" y=\"%.1f\" font-size=\"12\">", lx + 28.0, row + 4.0);
        svg_escape(fp, qty_upper);
        fputs("</text>\n", fp);
    }
}
