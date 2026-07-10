/*
 * Native SVG plot writer (prototype, issue #238).
 *
 * A zero-dependency "quick-look" renderer: pure fprintf into an SVG document, no
 * vendored library, no font engine, no raster pipeline — just a plot frame,
 * "nice" axis ticks, and one <polyline> per page.  It exists so a run can drop a
 * viewable 1-D curve next to its .bdo/.dat on a machine that has only the
 * openshieldhit binary.  The whole translation unit is compiled in only under
 * -DOSH_ENABLE_PLOT=ON, so the stock binary is unchanged.
 *
 * Two 1-D shapes are drawn, chosen automatically from the page model:
 *
 *   - Spatial profile — a depth-dose / Bragg curve or any 1-D profile: MESH with
 *     exactly one non-singleton X/Y/Z axis, or CYL with one non-singleton R/Z
 *     axis.  x = spatial coordinate (bin centres, cm).
 *   - Spectrum — a differential page (dPhi/dE, dose-vs-LET, ...) scored over a
 *     single spatial bin (a 0-D voxel or a single Zone): x = differential bin
 *     centres of the quantity's Diff1 axis, log-scaled when the binning is LOG.
 *     A single spatial bin is detected uniformly as diff_stride == 1, so a
 *     1x1x1 mesh and a one-zone geometry take the same path.
 *
 * Normalisation matches the ASCII writer: NORM/SUM quantities are divided by
 * nstat at write time, AVER quantities (DLET/TLET) are written as the physical
 * mean already produced by osh_scoring_postprocess() — which is also where a
 * differential page's bin-width division happens, so acc.data already holds
 * dPhi/dE here.  The plot is an additional artifact — BDO stays the source of
 * truth.  Save is a cold path, so the small malloc for the x-coordinate array is
 * fine (the §10 hot-path allocation ban does not apply).
 *
 * Out of scope for this prototype (all -> OSH_ENOTSUP): 2-D/3-D maps, 2-D
 * spectra (Diff1 x Diff2), spectra over more than one spatial bin, categorical
 * multi-zone profiles, rotated meshes.  PNG heatmaps and multipage PDF, plus
 * log-y and MC error bars, are follow-up items on #238.
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
    double xstep; /* linear x tick spacing (unused when xlog) */
    double ystep; /* linear y tick spacing */
    int xlog;     /* 0 = linear x, 1 = log10 x */
};

/* What to draw: the x-coordinate of each data point plus the two axis titles.
 * Built once (per mode) from the page model; owns its @c xs array. */
struct plot_spec {
    size_t npts;     /* number of data points (== page data length used) */
    double *xs;      /* [npts] x-coordinates in data space (owned) */
    int xlog;        /* 1 = log10 x-axis (log differential binning) */
    char xlabel[32]; /* x-axis title */
    char ylabel[96]; /* y-axis title */
};

static enum osh_status validate_plot_output(struct osh_scoring_workspace const *ws,
                                            struct osh_scoring_runtime const *rt,
                                            size_t output_idx,
                                            struct osh_scoring_output_runtime const **out_out,
                                            struct osh_scoring_geometry_runtime const **geo_out);
static enum osh_status build_profile_spec(struct osh_scoring_runtime const *rt,
                                          struct osh_scoring_output_runtime const *out,
                                          struct osh_scoring_geometry_runtime const *geo,
                                          struct plot_spec *spec);
static enum osh_status build_spectrum_spec(struct osh_scoring_runtime const *rt,
                                           struct osh_scoring_output_runtime const *out,
                                           struct plot_spec *spec);
static void plot_spec_free(struct plot_spec *spec);
static enum osh_status pick_profile_axis(struct osh_scoring_geometry_runtime const *geo, size_t *axis_idx_out);
static double page_scale(struct osh_scoring_page_runtime const *page, double inv_nstat);
static double diff_bin_center(struct osh_scoring_page_runtime const *page, size_t db);
static char *plot_output_path(char const *filename);
static char const *path_basename(char const *path);
static char const *plot_data_unit(struct osh_scoring_page_runtime const *page);
static char const *diff_kind_name(enum osh_scoring_diff_kind kind);
static char const *diff_kind_unit(enum osh_scoring_diff_kind kind);
static void page_value_unit_str(struct osh_scoring_page_runtime const *page, char *buf, size_t cap);
static void
build_ylabel(char *buf, size_t cap, struct osh_scoring_runtime const *rt, struct osh_scoring_output_runtime const *out);
static char const *series_color(size_t i);
static double nice_num(double range, int do_round);
static void nice_axis(double lo, double hi, int ntick, double *nlo, double *nhi, double *step);
static void nice_axis_log(double lo, double hi, double *nlo, double *nhi);
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
                             double const *xs,
                             size_t npts,
                             double const *data,
                             double scale,
                             char const *color);
static void svg_write_legend(FILE *fp,
                             struct plot_layout const *lay,
                             struct osh_scoring_runtime const *rt,
                             struct osh_scoring_output_runtime const *out);

static inline double _map_x(struct plot_layout const *lay, double x) {
    if (lay->xlog) {
        double a = log10(lay->xlo);
        double b = log10(lay->xhi);
        return lay->px0 + (log10(x) - a) / (b - a) * (lay->px1 - lay->px0);
    }
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
    struct osh_scoring_page_runtime const *p0;
    struct plot_spec spec;
    struct plot_layout lay;
    FILE *fp;
    char *out_path;
    double inv_nstat;
    double ymin;
    double ymax;
    double ylo_data;
    double xmin;
    double xmax;
    size_t ip;
    size_t i;
    char subtitle[128];
    enum osh_status rc;

    rc = validate_plot_output(ws, rt, output_idx, &out, &geo);
    if (rc != OSH_OK) {
        return rc;
    }
    if (nstat == 0ull) {
        return OSH_EINVAL;
    }

    /* A differential page 0 means a spectrum; otherwise a spatial profile. */
    memset(&spec, 0, sizeof(spec));
    p0 = &rt->pages[out->page_indices[0]];
    if (p0->diff_nbins > 0u) {
        rc = build_spectrum_spec(rt, out, &spec);
    } else {
        rc = build_profile_spec(rt, out, geo, &spec);
    }
    if (rc != OSH_OK) {
        plot_spec_free(&spec);
        return rc;
    }
    inv_nstat = 1.0 / (double) nstat;

    /* x data range from the coordinate array; y data range across every page.
     * Point i indexes acc.data[i] for both modes: a spatial profile has one
     * non-singleton axis, and a spectrum has diff_stride == 1, so in each case
     * the canonical flat index equals the point index. */
    xmin = spec.xs[0];
    xmax = spec.xs[0];
    for (i = 0; i < spec.npts; ++i) {
        if (spec.xs[i] < xmin) {
            xmin = spec.xs[i];
        }
        if (spec.xs[i] > xmax) {
            xmax = spec.xs[i];
        }
    }
    ymin = HUGE_VAL;
    ymax = -HUGE_VAL;
    for (ip = 0; ip < out->npages; ++ip) {
        struct osh_scoring_page_runtime const *page = &rt->pages[out->page_indices[ip]];
        double scale = page_scale(page, inv_nstat);
        for (i = 0; i < spec.npts; ++i) {
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
    lay.xlog = spec.xlog;
    if (spec.xlog) {
        nice_axis_log(xmin, xmax, &lay.xlo, &lay.xhi);
        lay.xstep = 0.0;
    } else {
        nice_axis(xmin, xmax, PLOT_NTICK, &lay.xlo, &lay.xhi, &lay.xstep);
    }
    nice_axis(ylo_data, ymax, PLOT_NTICK, &lay.ylo, &lay.yhi, &lay.ystep);

    snprintf(
        subtitle, sizeof(subtitle), "geometry: %s  \xe2\x80\x94  %llu primaries", geo->name ? geo->name : "?", nstat);

    out_path = plot_output_path(out->filename);
    if (!out_path) {
        plot_spec_free(&spec);
        return OSH_ENOMEM;
    }
    fp = fopen(out_path, "w");
    if (!fp) {
        free(out_path);
        plot_spec_free(&spec);
        return OSH_EIO;
    }

    svg_write_frame(fp, &lay, path_basename(out_path), subtitle, spec.xlabel, spec.ylabel);
    for (ip = 0; ip < out->npages; ++ip) {
        struct osh_scoring_page_runtime const *page = &rt->pages[out->page_indices[ip]];
        svg_write_series(fp, &lay, spec.xs, spec.npts, page->acc.data, page_scale(page, inv_nstat), series_color(ip));
    }
    svg_write_legend(fp, &lay, rt, out);
    fputs("</svg>\n", fp);

    free(out_path);
    plot_spec_free(&spec);
    if (fclose(fp) != 0) {
        return OSH_EIO;
    }
    return OSH_OK;
}

static void plot_spec_free(struct plot_spec *spec) {
    if (spec) {
        free(spec->xs);
        spec->xs = NULL;
    }
}

/* NORM/SUM quantities are divided by nstat; AVER quantities (DLET/TLET) already
 * hold a physical mean after postprocess and must not be. */
static double page_scale(struct osh_scoring_page_runtime const *page, double inv_nstat) {
    return (page->postproc == OSH_SCORING_POSTPROC_AVER) ? 1.0 : inv_nstat;
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

/* Spatial 1-D profile: MESH/CYL with exactly one non-singleton axis and no
 * differential pages.  Fills @p spec with the spatial bin centres (cm). */
static enum osh_status build_profile_spec(struct osh_scoring_runtime const *rt,
                                          struct osh_scoring_output_runtime const *out,
                                          struct osh_scoring_geometry_runtime const *geo,
                                          struct plot_spec *spec) {
    struct osh_scoring_axis_runtime const *axis;
    size_t axis_idx;
    size_t i;
    size_t ip;
    double lo;
    double dx;
    enum osh_status rc;

    if (geo->geo_kind != OSH_SCORING_GEO_MESH && geo->geo_kind != OSH_SCORING_GEO_CYL) {
        return OSH_ENOTSUP; /* ZONE / voxel spatial profiles are deferred */
    }
    for (ip = 0; ip < out->npages; ++ip) {
        if (rt->pages[out->page_indices[ip]].diff_nbins > 0u) {
            return OSH_ENOTSUP; /* mixed spatial + differential pages */
        }
    }
    rc = pick_profile_axis(geo, &axis_idx);
    if (rc != OSH_OK) {
        return rc;
    }
    axis = &geo->axes[axis_idx];
    spec->npts = (size_t) axis->nbins;
    spec->xs = (double *) malloc(spec->npts * sizeof(*spec->xs));
    if (!spec->xs) {
        return OSH_ENOMEM;
    }
    lo = axis->lo;
    dx = (axis->hi - axis->lo) / (double) spec->npts;
    for (i = 0; i < spec->npts; ++i) {
        spec->xs[i] = lo + ((double) i + 0.5) * dx;
    }
    spec->xlog = 0;
    snprintf(spec->xlabel, sizeof(spec->xlabel), "%s [cm]", axis->label);
    build_ylabel(spec->ylabel, sizeof(spec->ylabel), rt, out);
    return OSH_OK;
}

/* Differential spectrum over a single spatial bin (0-D voxel or single Zone).
 * Every page must share the same Diff1 layout, have no Diff2, and score into a
 * single spatial bin (diff_stride == 1).  Fills @p spec with the differential
 * bin centres, log-scaled when the binning is logarithmic. */
static enum osh_status build_spectrum_spec(struct osh_scoring_runtime const *rt,
                                           struct osh_scoring_output_runtime const *out,
                                           struct plot_spec *spec) {
    struct osh_scoring_page_runtime const *p0;
    char const *name;
    char const *unit;
    size_t i;
    size_t ip;

    p0 = &rt->pages[out->page_indices[0]];
    if (p0->diff_nbins == 0u || p0->diff2_nbins > 0u || p0->diff_stride != 1u) {
        return OSH_ENOTSUP; /* not a 1-D single-bin spectrum */
    }
    for (ip = 0; ip < out->npages; ++ip) {
        struct osh_scoring_page_runtime const *page = &rt->pages[out->page_indices[ip]];
        /* All pages must plot against the same differential axis. */
        if (page->diff_nbins != p0->diff_nbins || page->diff2_nbins != 0u || page->diff_stride != 1u
            || page->diff_kind != p0->diff_kind || page->diff_lo != p0->diff_lo || page->diff_hi != p0->diff_hi
            || page->diff_log != p0->diff_log) {
            return OSH_ENOTSUP;
        }
    }
    spec->npts = p0->diff_nbins;
    spec->xs = (double *) malloc(spec->npts * sizeof(*spec->xs));
    if (!spec->xs) {
        return OSH_ENOMEM;
    }
    for (i = 0; i < spec->npts; ++i) {
        spec->xs[i] = diff_bin_center(p0, i);
    }
    spec->xlog = p0->diff_log ? 1 : 0;
    name = diff_kind_name(p0->diff_kind);
    unit = diff_kind_unit(p0->diff_kind);
    if (unit[0]) {
        snprintf(spec->xlabel, sizeof(spec->xlabel), "%s [%s]", name, unit);
    } else {
        snprintf(spec->xlabel, sizeof(spec->xlabel), "%s", name);
    }
    build_ylabel(spec->ylabel, sizeof(spec->ylabel), rt, out);
    return OSH_OK;
}

/* Centre of differential bin @p db: arithmetic midpoint (linear binning) or
 * geometric midpoint (log binning).  Mirrors the ASCII writer. */
static double diff_bin_center(struct osh_scoring_page_runtime const *page, size_t db) {
    double t0 = (double) db / (double) page->diff_nbins;
    double t1 = (double) (db + 1u) / (double) page->diff_nbins;

    if (page->diff_log) {
        double ratio = page->diff_hi / page->diff_lo;
        return page->diff_lo * sqrt(pow(ratio, t0) * pow(ratio, t1));
    }
    return page->diff_lo + 0.5 * (t0 + t1) * (page->diff_hi - page->diff_lo);
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
 * writer.  The differential denominator (if any) is appended separately. */
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

/* Human-readable name of a differential axis quantity (for the x-axis title). */
static char const *diff_kind_name(enum osh_scoring_diff_kind kind) {
    switch (kind) {
    case OSH_SCORING_DIFF_EKIN:
        return "E";
    case OSH_SCORING_DIFF_ENUC:
        return "E/nucleon";
    case OSH_SCORING_DIFF_EAMU:
        return "E/amu";
    case OSH_SCORING_DIFF_LET:
        return "LET";
    case OSH_SCORING_DIFF_QEFF:
        return "(zeff/beta)^2";
    default:
        return "diff";
    }
}

/* Unit of a differential axis quantity, mirroring page_diff_unit() in the BDO
 * writer. */
static char const *diff_kind_unit(enum osh_scoring_diff_kind kind) {
    switch (kind) {
    case OSH_SCORING_DIFF_EKIN:
        return "MeV";
    case OSH_SCORING_DIFF_ENUC:
    case OSH_SCORING_DIFF_EAMU:
        return "MeV/u";
    case OSH_SCORING_DIFF_LET:
        return "MeV/cm";
    case OSH_SCORING_DIFF_QEFF:
        return "dim.less";
    default:
        return "";
    }
}

/* Compose a page's value unit into @p buf: the base quantity unit, plus the
 * differential denominator when the page is a spectrum (e.g. "/cm^2/MeV",
 * "/cm^2/(MeV/cm)").  Matches the BDO OSHBDO_PAG_DATA_UNIT string. */
static void page_value_unit_str(struct osh_scoring_page_runtime const *page, char *buf, size_t cap) {
    snprintf(buf, cap, "%s", plot_data_unit(page));
    if (page->diff_nbins > 0u) {
        char const *du = diff_kind_unit(page->diff_kind);
        size_t used = strlen(buf);
        if (du[0] && used < cap) {
            if (strchr(du, '/')) {
                snprintf(buf + used, cap - used, "/(%s)", du);
            } else {
                snprintf(buf + used, cap - used, "/%s", du);
            }
        }
    }
}

/* Build the y-axis title.  A single page reads "QUANTITY [unit]".  With several
 * pages the unit is shown only when every page agrees on it ("value [unit]"),
 * and omitted otherwise ("value") so the label never misrepresents a series. */
static void build_ylabel(char *buf,
                         size_t cap,
                         struct osh_scoring_runtime const *rt,
                         struct osh_scoring_output_runtime const *out) {
    struct osh_scoring_page_runtime const *p0 = &rt->pages[out->page_indices[0]];
    char unit0[48];
    char unitk[48];
    int shared;
    size_t ip;

    page_value_unit_str(p0, unit0, sizeof(unit0));
    shared = 1;
    for (ip = 1; ip < out->npages; ++ip) {
        page_value_unit_str(&rt->pages[out->page_indices[ip]], unitk, sizeof(unitk));
        if (strcmp(unit0, unitk) != 0) {
            shared = 0;
            break;
        }
    }

    if (out->npages == 1u) {
        char qty[64];
        upcase_copy(qty, sizeof(qty), p0->quantity ? p0->quantity : "value");
        snprintf(buf, cap, "%s [%s]", qty, unit0);
    } else if (shared) {
        snprintf(buf, cap, "value [%s]", unit0);
    } else {
        snprintf(buf, cap, "value");
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

/* Expand [lo, hi] outward to the enclosing powers of ten (log axis bounds).
 * @p lo must be positive; log binning guarantees it. */
static void nice_axis_log(double lo, double hi, double *nlo, double *nhi) {
    if (!(lo > 0.0)) {
        lo = 1.0;
    }
    if (hi <= lo) {
        hi = lo * 10.0;
    }
    *nlo = pow(10.0, floor(log10(lo)));
    *nhi = pow(10.0, ceil(log10(hi)));
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

    /* Vertical gridlines + x tick labels (log decades or linear nice ticks). */
    if (lay->xlog) {
        int e0 = (int) floor(log10(lay->xlo) + 0.5);
        int e1 = (int) floor(log10(lay->xhi) + 0.5);
        for (k = e0; k <= e1; ++k) {
            double v = pow(10.0, (double) k);
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
    } else {
        int nx = (int) floor((lay->xhi - lay->xlo) / lay->xstep + 0.5) + 1;
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
    }

    /* Horizontal gridlines + y tick labels. */
    ny = (int) floor((lay->yhi - lay->ylo) / lay->ystep + 0.5) + 1;
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
                             double const *xs,
                             size_t npts,
                             double const *data,
                             double scale,
                             char const *color) {
    size_t i;

    fprintf(fp, "<polyline fill=\"none\" stroke=\"%s\" stroke-width=\"2\" points=\"", color);
    for (i = 0; i < npts; ++i) {
        fprintf(fp, "%.2f,%.2f ", _map_x(lay, xs[i]), _map_y(lay, data[i] * scale));
    }
    fputs("\"/>\n", fp);

    if (npts <= PLOT_MARKER_MAX) {
        for (i = 0; i < npts; ++i) {
            fprintf(fp,
                    "<circle cx=\"%.2f\" cy=\"%.2f\" r=\"2.5\" fill=\"%s\"/>\n",
                    _map_x(lay, xs[i]),
                    _map_y(lay, data[i] * scale),
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
