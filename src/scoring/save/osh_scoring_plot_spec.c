/*
 * Plot spec builder (issue #238) — the scoring-model half of the SVG writer.
 *
 * Classifies one compiled Output and selects the pages that form a 1-D curve:
 * a spatial profile (MESH/CYL with one non-singleton axis) or a single-bin
 * spectrum (a Diff1 page over a 0-D voxel or a single Zone).  Mixed outputs keep
 * only the pages that fit the chosen shape.  Everything here operates on the
 * page/geometry model; the drawing lives in osh_scoring_svg.c.
 */

#include "scoring/save/osh_scoring_plot_spec.h"

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static enum osh_status build_profile_spec(struct osh_scoring_runtime const *rt,
                                          struct osh_scoring_output_runtime const *out,
                                          struct osh_scoring_geometry_runtime const *geo,
                                          struct osh_plot_spec *spec);
static enum osh_status build_spectrum_spec(struct osh_scoring_runtime const *rt,
                                           struct osh_scoring_output_runtime const *out,
                                           struct osh_plot_spec *spec);
static int spectrum_page_matches(struct osh_scoring_page_runtime const *ref,
                                 struct osh_scoring_page_runtime const *page);
static enum osh_status pick_profile_axis(struct osh_scoring_geometry_runtime const *geo, size_t *axis_idx_out);
static double diff_bin_center(struct osh_scoring_page_runtime const *page, size_t db);
static char const *plot_data_unit(struct osh_scoring_page_runtime const *page);
static char const *diff_kind_name(enum osh_scoring_diff_kind kind);
static char const *diff_kind_unit(enum osh_scoring_diff_kind kind);
static void page_value_unit_str(struct osh_scoring_page_runtime const *page, char *buf, size_t cap);
static void build_ylabel(char *buf, size_t cap, struct osh_scoring_runtime const *rt, size_t const *pages, size_t nsel);
static void upcase_copy(char *dst, size_t cap, char const *src);

enum osh_status osh_plot_spec_build(struct osh_scoring_runtime const *rt,
                                    struct osh_scoring_output_runtime const *out,
                                    struct osh_scoring_geometry_runtime const *geo,
                                    struct osh_plot_spec *spec) {
    enum osh_status rc;

    memset(spec, 0, sizeof(*spec));
    rc = build_profile_spec(rt, out, geo, spec);
    if (rc == OSH_ENOTSUP) {
        rc = build_spectrum_spec(rt, out, spec);
    }
    return rc;
}

void osh_plot_spec_free(struct osh_plot_spec *spec) {
    if (spec) {
        free(spec->xs);
        free(spec->pages);
        spec->xs = NULL;
        spec->pages = NULL;
    }
}

/* NORM/SUM quantities are divided by nstat; AVER quantities (DLET/TLET) already
 * hold a physical mean after postprocess and must not be. */
double osh_plot_page_scale(struct osh_scoring_page_runtime const *page, double inv_nstat) {
    return (page->postproc == OSH_SCORING_POSTPROC_AVER) ? 1.0 : inv_nstat;
}

void osh_plot_spec_series_label(
    struct osh_scoring_runtime const *rt, struct osh_plot_spec const *spec, size_t k, char *buf, size_t cap) {
    struct osh_scoring_page_runtime const *page = &rt->pages[spec->pages[k]];
    upcase_copy(buf, cap, page->quantity ? page->quantity : "?");
}

/* Spatial 1-D profile: MESH/CYL with exactly one non-singleton axis.  Selects
 * the plain (non-differential) pages and fills @p spec with the spatial bin
 * centres (cm); differential pages on the same Output are 2-D here and skipped.
 * Returns OSH_ENOTSUP when the geometry is not a 1-D profile or no plain page
 * exists (leaving the caller free to try the spectrum shape). */
static enum osh_status build_profile_spec(struct osh_scoring_runtime const *rt,
                                          struct osh_scoring_output_runtime const *out,
                                          struct osh_scoring_geometry_runtime const *geo,
                                          struct osh_plot_spec *spec) {
    struct osh_scoring_axis_runtime const *axis;
    size_t axis_idx;
    size_t i;
    size_t ip;
    double lo;
    double dx;
    enum osh_status rc;

    if (geo->geo_kind != OSH_SCORING_GEO_MESH && geo->geo_kind != OSH_SCORING_GEO_CYL) {
        return OSH_ENOTSUP; /* ZONE / voxel handled as a spectrum */
    }
    rc = pick_profile_axis(geo, &axis_idx);
    if (rc != OSH_OK) {
        return rc;
    }

    spec->pages = (size_t *) malloc(out->npages * sizeof(*spec->pages));
    if (!spec->pages) {
        return OSH_ENOMEM;
    }
    spec->nsel = 0u;
    for (ip = 0; ip < out->npages; ++ip) {
        size_t page_idx = out->page_indices[ip];
        if (rt->pages[page_idx].diff_nbins == 0u) {
            spec->pages[spec->nsel++] = page_idx;
        }
    }
    if (spec->nsel == 0u) {
        free(spec->pages);
        spec->pages = NULL;
        return OSH_ENOTSUP; /* every page is differential -> 2-D on this geometry */
    }

    axis = &geo->axes[axis_idx];
    spec->npts = (size_t) axis->nbins;
    spec->xs = (double *) malloc(spec->npts * sizeof(*spec->xs));
    if (!spec->xs) {
        free(spec->pages);
        spec->pages = NULL;
        return OSH_ENOMEM;
    }
    lo = axis->lo;
    dx = (axis->hi - axis->lo) / (double) spec->npts;
    for (i = 0; i < spec->npts; ++i) {
        spec->xs[i] = lo + ((double) i + 0.5) * dx;
    }
    spec->xlog = 0;
    snprintf(spec->xlabel, sizeof(spec->xlabel), "%s [cm]", axis->label);
    build_ylabel(spec->ylabel, sizeof(spec->ylabel), rt, spec->pages, spec->nsel);
    return OSH_OK;
}

/* Differential spectrum over a single spatial bin (0-D voxel or single Zone).
 * The first page that is a valid single-bin Diff1 (diff_stride == 1, no Diff2)
 * fixes the shared differential axis; every page matching that axis is plotted
 * and the rest are skipped.  Fills @p spec with the Diff1 bin centres,
 * log-scaled when the binning is logarithmic. */
static enum osh_status build_spectrum_spec(struct osh_scoring_runtime const *rt,
                                           struct osh_scoring_output_runtime const *out,
                                           struct osh_plot_spec *spec) {
    struct osh_scoring_page_runtime const *ref = NULL;
    char const *name;
    char const *unit;
    size_t i;
    size_t ip;

    for (ip = 0; ip < out->npages; ++ip) {
        struct osh_scoring_page_runtime const *page = &rt->pages[out->page_indices[ip]];
        if (page->diff_nbins > 0u && page->diff2_nbins == 0u && page->diff_stride == 1u) {
            ref = page;
            break;
        }
    }
    if (!ref) {
        return OSH_ENOTSUP; /* no 1-D single-bin spectrum page */
    }

    spec->pages = (size_t *) malloc(out->npages * sizeof(*spec->pages));
    if (!spec->pages) {
        return OSH_ENOMEM;
    }
    spec->nsel = 0u;
    for (ip = 0; ip < out->npages; ++ip) {
        size_t page_idx = out->page_indices[ip];
        if (spectrum_page_matches(ref, &rt->pages[page_idx])) {
            spec->pages[spec->nsel++] = page_idx;
        }
    }

    spec->npts = ref->diff_nbins;
    spec->xs = (double *) malloc(spec->npts * sizeof(*spec->xs));
    if (!spec->xs) {
        free(spec->pages);
        spec->pages = NULL;
        return OSH_ENOMEM;
    }
    for (i = 0; i < spec->npts; ++i) {
        spec->xs[i] = diff_bin_center(ref, i);
    }
    spec->xlog = ref->diff_log ? 1 : 0;
    name = diff_kind_name(ref->diff_kind);
    unit = diff_kind_unit(ref->diff_kind);
    if (unit[0]) {
        snprintf(spec->xlabel, sizeof(spec->xlabel), "%s [%s]", name, unit);
    } else {
        snprintf(spec->xlabel, sizeof(spec->xlabel), "%s", name);
    }
    build_ylabel(spec->ylabel, sizeof(spec->ylabel), rt, spec->pages, spec->nsel);
    return OSH_OK;
}

/* True when @p page can share @p ref's differential x-axis (same single-bin
 * Diff1 layout, no Diff2). */
static int spectrum_page_matches(struct osh_scoring_page_runtime const *ref,
                                 struct osh_scoring_page_runtime const *page) {
    return page->diff_nbins == ref->diff_nbins && page->diff2_nbins == 0u && page->diff_stride == 1u
           && page->diff_kind == ref->diff_kind && page->diff_lo == ref->diff_lo && page->diff_hi == ref->diff_hi
           && page->diff_log == ref->diff_log;
}

/* Find the single non-singleton spatial axis (the profile axis).  Returns
 * OSH_ENOTSUP unless exactly one axis has nbins > 1: zero means a single point
 * (handled as a spectrum), more than one means a 2-D/3-D map (deferred).
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

/* Build the y-axis title from the plotted pages.  A single page reads
 * "QUANTITY [unit]".  With several pages the unit is shown only when every page
 * agrees on it ("value [unit]"), and omitted otherwise ("value") so the label
 * never misrepresents a series. */
static void
build_ylabel(char *buf, size_t cap, struct osh_scoring_runtime const *rt, size_t const *pages, size_t nsel) {
    struct osh_scoring_page_runtime const *p0 = &rt->pages[pages[0]];
    char unit0[48];
    char unitk[48];
    int shared;
    size_t k;

    page_value_unit_str(p0, unit0, sizeof(unit0));
    shared = 1;
    for (k = 1; k < nsel; ++k) {
        page_value_unit_str(&rt->pages[pages[k]], unitk, sizeof(unitk));
        if (strcmp(unit0, unitk) != 0) {
            shared = 0;
            break;
        }
    }

    if (nsel == 1u) {
        char qty[64];
        upcase_copy(qty, sizeof(qty), p0->quantity ? p0->quantity : "value");
        snprintf(buf, cap, "%s [%s]", qty, unit0);
    } else if (shared) {
        snprintf(buf, cap, "value [%s]", unit0);
    } else {
        snprintf(buf, cap, "value");
    }
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
