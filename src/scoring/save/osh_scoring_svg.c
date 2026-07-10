/*
 * Generic SVG line-plot renderer (issue #238).
 *
 * Scoring-agnostic: everything here operates on plain doubles and strings, so
 * the same primitives can later back other line plots.  See osh_scoring_svg.h
 * for the contract; osh_scoring_save_plot.c is the scoring-aware caller.
 */

#include "scoring/save/osh_scoring_svg.h"

#include <math.h>
#include <string.h>

/* Plot-area margins inside the fixed OSH_SVG_W x OSH_SVG_H canvas. */
#define SVG_MARGIN_L 74
#define SVG_MARGIN_R 30
#define SVG_MARGIN_T 48
#define SVG_MARGIN_B 58
#define SVG_NTICK 6       /* target major tick count per axis */
#define SVG_MINOR_DIV 5   /* minor grid subdivisions per linear major interval */
#define SVG_MARKER_MAX 40 /* draw per-point markers only up to this many points */

/* Grid styling. */
#define SVG_GRID_MINOR "#f0f0f0"
#define SVG_GRID_MAJOR "#dcdcdc"
#define SVG_GRID_MINOR_W 0.5
#define SVG_GRID_MAJOR_W 1.0

/* A single data series is expected per plot; draw it in black. */
#define SVG_SERIES_COLOR "#000000"

static double nice_num(double range, int do_round);
static void nice_axis(double lo, double hi, int ntick, double *nlo, double *nhi, double *step);
static void nice_axis_log(double lo, double hi, double *nlo, double *nhi);
static void svg_escape(FILE *fp, char const *s);
static void svg_line(FILE *fp, double x0, double y0, double x1, double y1, char const *stroke, double width);
static void draw_titles(FILE *fp, struct osh_svg_plot const *p, char const *title, char const *subtitle);
static void draw_x_grid(FILE *fp, struct osh_svg_plot const *p);
static void draw_y_grid(FILE *fp, struct osh_svg_plot const *p);
static void draw_border_and_axis_titles(FILE *fp, struct osh_svg_plot const *p, char const *xlabel, char const *ylabel);

static inline double map_x(struct osh_svg_plot const *p, double x) {
    if (p->xlog) {
        double a = log10(p->xlo);
        double b = log10(p->xhi);
        return p->px0 + (log10(x) - a) / (b - a) * (p->px1 - p->px0);
    }
    return p->px0 + (x - p->xlo) / (p->xhi - p->xlo) * (p->px1 - p->px0);
}

static inline double map_y(struct osh_svg_plot const *p, double y) {
    /* Larger data values map to smaller pixel-y (top of the canvas). */
    return p->py1 - (y - p->ylo) / (p->yhi - p->ylo) * (p->py1 - p->py0);
}

void osh_svg_plot_init(struct osh_svg_plot *p, double xmin, double xmax, int xlog, double ymin, double ymax) {
    p->px0 = (double) SVG_MARGIN_L;
    p->px1 = (double) (OSH_SVG_W - SVG_MARGIN_R);
    p->py0 = (double) SVG_MARGIN_T;
    p->py1 = (double) (OSH_SVG_H - SVG_MARGIN_B);
    if (xlog) {
        p->xlog = 1;
    } else {
        p->xlog = 0;
    }
    if (p->xlog) {
        nice_axis_log(xmin, xmax, &p->xlo, &p->xhi);
        p->xstep = 0.0;
    } else {
        nice_axis(xmin, xmax, SVG_NTICK, &p->xlo, &p->xhi, &p->xstep);
    }
    nice_axis(ymin, ymax, SVG_NTICK, &p->ylo, &p->yhi, &p->ystep);
}

void osh_svg_begin(FILE *fp,
                   struct osh_svg_plot const *p,
                   char const *title,
                   char const *subtitle,
                   char const *xlabel,
                   char const *ylabel) {
    fprintf(fp, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    fprintf(fp,
            "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%d\" height=\"%d\" "
            "viewBox=\"0 0 %d %d\" font-family=\"sans-serif\">\n",
            OSH_SVG_W,
            OSH_SVG_H,
            OSH_SVG_W,
            OSH_SVG_H);
    fprintf(fp, "<rect width=\"%d\" height=\"%d\" fill=\"#ffffff\"/>\n", OSH_SVG_W, OSH_SVG_H);

    draw_titles(fp, p, title, subtitle);
    draw_x_grid(fp, p);
    draw_y_grid(fp, p);
    draw_border_and_axis_titles(fp, p, xlabel, ylabel);
}

void osh_svg_series(FILE *fp, struct osh_svg_plot const *p, double const *xs, double const *ys, size_t npts) {
    size_t i;

    fprintf(fp, "<polyline fill=\"none\" stroke=\"%s\" stroke-width=\"2\" points=\"", SVG_SERIES_COLOR);
    for (i = 0; i < npts; ++i) {
        fprintf(fp, "%.2f,%.2f ", map_x(p, xs[i]), map_y(p, ys[i]));
    }
    fputs("\"/>\n", fp);

    if (npts <= SVG_MARKER_MAX) {
        for (i = 0; i < npts; ++i) {
            fprintf(fp,
                    "<circle cx=\"%.2f\" cy=\"%.2f\" r=\"2.5\" fill=\"%s\"/>\n",
                    map_x(p, xs[i]),
                    map_y(p, ys[i]),
                    SVG_SERIES_COLOR);
        }
    }
}

void osh_svg_end(FILE *fp) {
    fputs("</svg>\n", fp);
}

static void draw_titles(FILE *fp, struct osh_svg_plot const *p, char const *title, char const *subtitle) {
    double cx = 0.5 * (p->px0 + p->px1);

    fprintf(fp, "<text x=\"%.1f\" y=\"22\" text-anchor=\"middle\" font-size=\"16\" font-weight=\"bold\">", cx);
    svg_escape(fp, title);
    fputs("</text>\n", fp);
    fprintf(fp, "<text x=\"%.1f\" y=\"39\" text-anchor=\"middle\" font-size=\"11\" fill=\"#555555\">", cx);
    svg_escape(fp, subtitle);
    fputs("</text>\n", fp);
}

/* Vertical grid: minor lines first (thin), then major lines + x tick labels. */
static void draw_x_grid(FILE *fp, struct osh_svg_plot const *p) {
    int k;
    int m;

    if (p->xlog) {
        int e0 = (int) floor(log10(p->xlo) + 0.5);
        int e1 = (int) floor(log10(p->xhi) + 0.5);
        for (k = e0; k < e1; ++k) {
            for (m = 2; m <= 9; ++m) {
                double px = map_x(p, (double) m * pow(10.0, (double) k));
                svg_line(fp, px, p->py0, px, p->py1, SVG_GRID_MINOR, SVG_GRID_MINOR_W);
            }
        }
        for (k = e0; k <= e1; ++k) {
            double v = pow(10.0, (double) k);
            double px = map_x(p, v);
            svg_line(fp, px, p->py0, px, p->py1, SVG_GRID_MAJOR, SVG_GRID_MAJOR_W);
            fprintf(fp,
                    "<text x=\"%.1f\" y=\"%.1f\" text-anchor=\"middle\" font-size=\"11\">%g</text>\n",
                    px,
                    p->py1 + 16.0,
                    v);
        }
    } else {
        int nx = (int) floor((p->xhi - p->xlo) / p->xstep + 0.5) + 1;
        for (k = 0; k + 1 < nx; ++k) {
            for (m = 1; m < SVG_MINOR_DIV; ++m) {
                double v = p->xlo + ((double) k + (double) m / SVG_MINOR_DIV) * p->xstep;
                double px = map_x(p, v);
                svg_line(fp, px, p->py0, px, p->py1, SVG_GRID_MINOR, SVG_GRID_MINOR_W);
            }
        }
        for (k = 0; k < nx; ++k) {
            double v = p->xlo + (double) k * p->xstep;
            double px = map_x(p, v);
            svg_line(fp, px, p->py0, px, p->py1, SVG_GRID_MAJOR, SVG_GRID_MAJOR_W);
            fprintf(fp,
                    "<text x=\"%.1f\" y=\"%.1f\" text-anchor=\"middle\" font-size=\"11\">%g</text>\n",
                    px,
                    p->py1 + 16.0,
                    v);
        }
    }
}

/* Horizontal grid: minor lines first (thin), then major lines + y tick labels. */
static void draw_y_grid(FILE *fp, struct osh_svg_plot const *p) {
    int ny = (int) floor((p->yhi - p->ylo) / p->ystep + 0.5) + 1;
    int k;
    int m;

    for (k = 0; k + 1 < ny; ++k) {
        for (m = 1; m < SVG_MINOR_DIV; ++m) {
            double v = p->ylo + ((double) k + (double) m / SVG_MINOR_DIV) * p->ystep;
            double py = map_y(p, v);
            svg_line(fp, p->px0, py, p->px1, py, SVG_GRID_MINOR, SVG_GRID_MINOR_W);
        }
    }
    for (k = 0; k < ny; ++k) {
        double v = p->ylo + (double) k * p->ystep;
        double py = map_y(p, v);
        svg_line(fp, p->px0, py, p->px1, py, SVG_GRID_MAJOR, SVG_GRID_MAJOR_W);
        fprintf(fp,
                "<text x=\"%.1f\" y=\"%.1f\" text-anchor=\"end\" font-size=\"11\">%g</text>\n",
                p->px0 - 8.0,
                py + 4.0,
                v);
    }
}

static void
draw_border_and_axis_titles(FILE *fp, struct osh_svg_plot const *p, char const *xlabel, char const *ylabel) {
    double cx = 0.5 * (p->px0 + p->px1);
    double cy = 0.5 * (p->py0 + p->py1);

    fprintf(fp,
            "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"none\" stroke=\"#333333\"/>\n",
            p->px0,
            p->py0,
            p->px1 - p->px0,
            p->py1 - p->py0);
    fprintf(fp, "<text x=\"%.1f\" y=\"%d\" text-anchor=\"middle\" font-size=\"13\">", cx, OSH_SVG_H - 16);
    svg_escape(fp, xlabel);
    fputs("</text>\n", fp);
    fprintf(fp,
            "<text x=\"18\" y=\"%.1f\" text-anchor=\"middle\" font-size=\"13\" transform=\"rotate(-90 18 %.1f)\">",
            cy,
            cy);
    svg_escape(fp, ylabel);
    fputs("</text>\n", fp);
}

static void svg_line(FILE *fp, double x0, double y0, double x1, double y1, char const *stroke, double width) {
    fprintf(fp,
            "<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"%s\" stroke-width=\"%.1f\"/>\n",
            x0,
            y0,
            x1,
            y1,
            stroke,
            width);
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
