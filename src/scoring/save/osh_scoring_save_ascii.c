/*
 * Plain-text (ASCII) scorer output writer.
 *
 * Normalisation contract
 * ----------------------
 * Values are normalised per primary particle at write time, making the output
 * immediately human-readable without post-processing:
 *
 *   NORM quantities (DOSE, FLUENCE, ENERGY, …):
 *       value_per_primary = data / nstat
 *
 *   AVER quantities (DLET, TLET, …):
 *       averaged_value = data     (already a physical mean after
 *                                  osh_scoring_postprocess(); nstat division
 *                                  would be incorrect here)
 *
 *   SUM quantities (COUNT, …):
 *       total = data / nstat      (normalised count)
 *
 * ASCII output is intended for quick single-run inspection.  It is NOT suited
 * for multi-run merging: once divided by nstat the absolute weight of each
 * run is lost.  Use BDO format for production work and multi-run accumulation.
 */

#include "scoring/save/osh_scoring_save_ascii.h"

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "common/osh_version.h"

static enum osh_status
mesh_axis_index(struct osh_scoring_geometry_runtime const *geo, char const *label, size_t *idx_out);
static enum osh_status validate_output(struct osh_scoring_workspace const *ws,
                                       struct osh_scoring_runtime const *rt,
                                       size_t output_idx,
                                       struct osh_scoring_output_runtime const **out_out,
                                       struct osh_scoring_geometry_runtime const **geo_out);
static void format_now_rfc2822(char *buf, size_t cap);
static void
fprint_quantity_names(FILE *fp, struct osh_scoring_runtime const *rt, struct osh_scoring_output_runtime const *out);
static double ascii_diff_center(struct osh_scoring_page_runtime const *page, size_t db);
static double ascii_diff2_center(struct osh_scoring_page_runtime const *page, size_t db);
static char const *diff_kind_label(enum osh_scoring_diff_kind kind);

enum osh_status osh_scoring_save_ascii_output(struct osh_scoring_workspace const *ws,
                                              struct osh_scoring_runtime const *rt,
                                              unsigned long long nstat,
                                              size_t output_idx) {
    FILE *fp;
    char datestr[128];
    enum osh_status rc;
    struct osh_scoring_output_runtime const *out;
    struct osh_scoring_geometry_runtime const *geo;
    double inv_nstat;
    size_t ip;

    rc = validate_output(ws, rt, output_idx, &out, &geo);
    if (rc != OSH_OK) {
        return rc;
    }
    if (nstat == 0ull) {
        return OSH_EINVAL;
    }

    fp = fopen(out->filename, "w");
    if (!fp) {
        return OSH_EIO;
    }

    inv_nstat = 1.0 / (double) nstat;
    format_now_rfc2822(datestr, sizeof(datestr));

    if (geo->geo_kind == OSH_SCORING_GEO_ZONE) {
        /* --- ZONE output --- */
        size_t izone;

        {
            struct osh_scoring_page_runtime const *p0;
            size_t diff_nbins;
            size_t diff2_nbins;

            if (out->npages == 0u) {
                fclose(fp);
                return OSH_OK;
            }
            p0 = &rt->pages[out->page_indices[0]];
            diff_nbins = p0->diff_nbins;
            diff2_nbins = p0->diff2_nbins;

            fprintf(fp, "# OpenShieldHIT version %s\n", OSH_VERSION);
            fprintf(fp, "# Calculated %s\n", datestr);
            fprintf(fp, "# DETECTOR OUTPUT ZONE\n");
            fprintf(fp, "# ZONE BIN: %5zu\n", geo->nzone_indices);
            fprintf(fp, "# DETECTOR TYPE:");
            fprint_quantity_names(fp, rt, out);
            fputc('\n', fp);
            fprintf(fp, "# PRIMARIES: %llu\n", nstat);
            fprintf(fp, "# COMPLETENESS: %s\n", osh_scoring_runtime_completeness_label(rt));
            fprintf(fp, "# Data written in explicit Zone order from detect.dat\n");
            fprintf(fp, "# Values: NORM/SUM quantities divided by nstat; AVER quantities written as physical mean\n");
            if (diff_nbins > 0u) {
                fprintf(fp,
                        "# Diff1Type: %s  lo=%g  hi=%g  nbins=%zu%s%s\n",
                        diff_kind_label(p0->diff_kind),
                        p0->diff_lo,
                        p0->diff_hi,
                        p0->diff_nbins,
                        p0->diff_log ? " LOG" : "",
                        p0->has_diff_sset ? "  (SP override active)" : "");
            }
            if (diff2_nbins > 0u) {
                fprintf(fp,
                        "# Diff2Type: %s  lo=%g  hi=%g  nbins=%zu%s%s\n",
                        diff_kind_label(p0->diff2_kind),
                        p0->diff2_lo,
                        p0->diff2_hi,
                        p0->diff2_nbins,
                        p0->diff2_log ? " LOG" : "",
                        p0->has_diff2_sset ? "  (SP override active)" : "");
            }
            fprintf(fp, "# ZONE");
            if (diff_nbins > 0u) {
                fprintf(fp, " %s", diff_kind_label(p0->diff_kind));
            }
            if (diff2_nbins > 0u) {
                fprintf(fp, " %s", diff_kind_label(p0->diff2_kind));
            }
            fprint_quantity_names(fp, rt, out);
            fputc('\n', fp);

            for (izone = 0u; izone < geo->nzone_indices; ++izone) {
                size_t db;
                size_t db2;
                size_t ndb = (diff_nbins > 0u) ? diff_nbins : 1u;
                size_t ndb2 = (diff2_nbins > 0u) ? diff2_nbins : 1u;
                for (db = 0u; db < ndb; ++db) {
                    for (db2 = 0u; db2 < ndb2; ++db2) {
                        fprintf(fp, " %zu", geo->zone_indices[izone]);
                        if (diff_nbins > 0u) {
                            fprintf(fp, " %.12e", ascii_diff_center(p0, db));
                        }
                        if (diff2_nbins > 0u) {
                            fprintf(fp, " %.12e", ascii_diff2_center(p0, db2));
                        }
                        for (ip = 0u; ip < out->npages; ++ip) {
                            size_t page_idx = out->page_indices[ip];
                            struct osh_scoring_page_runtime const *page = &rt->pages[page_idx];
                            double scale = (page->postproc == OSH_SCORING_POSTPROC_AVER) ? 1.0 : inv_nstat;
                            size_t data_idx = izone + (page->diff_nbins > 0u ? db * page->diff_stride : 0u)
                                              + (page->diff2_nbins > 0u ? db2 * page->diff2_stride : 0u);
                            fprintf(fp, " %.12e", page->acc.data[data_idx] * scale);
                        }
                        fprintf(fp, "\n");
                    }
                }
            }
        }
    } else if (geo->geo_kind == OSH_SCORING_GEO_CYL) {
        /* --- CYL output --- */
        size_t ir_axis;
        size_t iz_axis;
        size_t nr;
        size_t nz;
        double r0;
        double z0;
        double dr;
        double dz;
        size_t ir;
        size_t iz;

        rc = mesh_axis_index(geo, "R", &ir_axis);
        if (rc != OSH_OK) {
            fclose(fp);
            return rc;
        }
        rc = mesh_axis_index(geo, "Z", &iz_axis);
        if (rc != OSH_OK) {
            fclose(fp);
            return rc;
        }

        nr = (size_t) geo->axes[ir_axis].nbins;
        nz = (size_t) geo->axes[iz_axis].nbins;
        r0 = geo->axes[ir_axis].lo;
        z0 = geo->axes[iz_axis].lo;
        dr = (geo->axes[ir_axis].hi - r0) / (double) nr;
        dz = (geo->axes[iz_axis].hi - z0) / (double) nz;

        {
            struct osh_scoring_page_runtime const *p0;
            size_t diff_nbins;
            size_t diff2_nbins;

            if (out->npages == 0u) {
                fclose(fp);
                return OSH_OK;
            }
            p0 = &rt->pages[out->page_indices[0]];
            diff_nbins = p0->diff_nbins;
            diff2_nbins = p0->diff2_nbins;

            fprintf(fp, "# OpenShieldHIT version %s\n", OSH_VERSION);
            fprintf(fp, "# Calculated %s\n", datestr);
            fprintf(fp, "# DETECTOR OUTPUT CYL\n");
            fprintf(fp, "# R BIN: %5zu Z BIN: %5zu\n", nr, nz);
            fprintf(fp, "# DETECTOR TYPE:");
            fprint_quantity_names(fp, rt, out);
            fputc('\n', fp);
            fprintf(fp, "# R START: %12.6E Z START: %12.6E\n", r0, z0);
            fprintf(fp, "# R END  : %12.6E Z END  : %12.6E\n", geo->axes[ir_axis].hi, geo->axes[iz_axis].hi);
            fprintf(fp, "# PRIMARIES: %llu\n", nstat);
            fprintf(fp, "# COMPLETENESS: %s\n", osh_scoring_runtime_completeness_label(rt));
            fprintf(fp, "# Data written in canonical flat cyl order: idx = ir + nr * iz\n");
            fprintf(fp,
                    "# Values: NORM/SUM quantities divided by nstat; AVER quantities (DLET/TLET) written as physical "
                    "mean\n");
            if (diff_nbins > 0u) {
                fprintf(fp,
                        "# Diff1Type: %s  lo=%g  hi=%g  nbins=%zu%s%s\n",
                        diff_kind_label(p0->diff_kind),
                        p0->diff_lo,
                        p0->diff_hi,
                        p0->diff_nbins,
                        p0->diff_log ? " LOG" : "",
                        p0->has_diff_sset ? "  (SP override active)" : "");
            }
            if (diff2_nbins > 0u) {
                fprintf(fp,
                        "# Diff2Type: %s  lo=%g  hi=%g  nbins=%zu%s%s\n",
                        diff_kind_label(p0->diff2_kind),
                        p0->diff2_lo,
                        p0->diff2_hi,
                        p0->diff2_nbins,
                        p0->diff2_log ? " LOG" : "",
                        p0->has_diff2_sset ? "  (SP override active)" : "");
            }
            fprintf(fp, "# Z R");
            if (diff_nbins > 0u) {
                fprintf(fp, " %s", diff_kind_label(p0->diff_kind));
            }
            if (diff2_nbins > 0u) {
                fprintf(fp, " %s", diff_kind_label(p0->diff2_kind));
            }
            fprint_quantity_names(fp, rt, out);
            fputc('\n', fp);

            for (iz = 0; iz < nz; ++iz) {
                for (ir = 0; ir < nr; ++ir) {
                    size_t spatial_idx = ir + nr * iz;
                    size_t db;
                    size_t db2;
                    size_t ndb = (diff_nbins > 0u) ? diff_nbins : 1u;
                    size_t ndb2 = (diff2_nbins > 0u) ? diff2_nbins : 1u;
                    for (db = 0; db < ndb; ++db) {
                        for (db2 = 0; db2 < ndb2; ++db2) {
                            fprintf(fp, " %.12e %.12e", z0 + dz * ((double) iz + 0.5), r0 + dr * ((double) ir + 0.5));
                            if (diff_nbins > 0u) {
                                fprintf(fp, " %.12e", ascii_diff_center(p0, db));
                            }
                            if (diff2_nbins > 0u) {
                                fprintf(fp, " %.12e", ascii_diff2_center(p0, db2));
                            }
                            for (ip = 0; ip < out->npages; ++ip) {
                                size_t page_idx = out->page_indices[ip];
                                struct osh_scoring_page_runtime const *page = &rt->pages[page_idx];
                                double scale = (page->postproc == OSH_SCORING_POSTPROC_AVER) ? 1.0 : inv_nstat;
                                /* Use each page's own diff_nbins to gate the offset: non-differential
                                 * pages (diff_nbins==0) have len==geo_nbins and must not use db/db2
                                 * offsets, even when other pages in the same output are differential.
                                 * Their scalar value is repeated across all diff-bin rows. */
                                size_t data_idx = spatial_idx + (page->diff_nbins > 0u ? db * page->diff_stride : 0u)
                                                  + (page->diff2_nbins > 0u ? db2 * page->diff2_stride : 0u);
                                fprintf(fp, " %.12e", page->acc.data[data_idx] * scale);
                            }
                            fprintf(fp, "\n");
                        }
                    }
                }
            }
        }
    } else {
        /* --- MESH output --- */
        size_t ix_axis;
        size_t iy_axis;
        size_t iz_axis;
        size_t nx;
        size_t ny;
        size_t nz;
        double x0;
        double y0;
        double z0;
        double dx;
        double dy;
        double dz;
        size_t ix;
        size_t iy;
        size_t iz;

        rc = mesh_axis_index(geo, "X", &ix_axis);
        if (rc != OSH_OK) {
            fclose(fp);
            return rc;
        }
        rc = mesh_axis_index(geo, "Y", &iy_axis);
        if (rc != OSH_OK) {
            fclose(fp);
            return rc;
        }
        rc = mesh_axis_index(geo, "Z", &iz_axis);
        if (rc != OSH_OK) {
            fclose(fp);
            return rc;
        }

        nx = (size_t) geo->axes[ix_axis].nbins;
        ny = (size_t) geo->axes[iy_axis].nbins;
        nz = (size_t) geo->axes[iz_axis].nbins;
        x0 = geo->axes[ix_axis].lo;
        y0 = geo->axes[iy_axis].lo;
        z0 = geo->axes[iz_axis].lo;
        dx = (geo->axes[ix_axis].hi - x0) / (double) nx;
        dy = (geo->axes[iy_axis].hi - y0) / (double) ny;
        dz = (geo->axes[iz_axis].hi - z0) / (double) nz;

        {
            struct osh_scoring_page_runtime const *p0;
            size_t diff_nbins;
            size_t diff2_nbins;

            if (out->npages == 0u) {
                fclose(fp);
                return OSH_OK;
            }
            p0 = &rt->pages[out->page_indices[0]];
            diff_nbins = p0->diff_nbins;
            diff2_nbins = p0->diff2_nbins;

            fprintf(fp, "# OpenShieldHIT version %s\n", OSH_VERSION);
            fprintf(fp, "# Calculated %s\n", datestr);
            fprintf(fp, "# DETECTOR OUTPUT MSH\n");
            fprintf(fp, "# X BIN: %5zu Y BIN: %5zu Z BIN: %5zu\n", nx, ny, nz);
            fprintf(fp, "# DETECTOR TYPE:");
            fprint_quantity_names(fp, rt, out);
            fputc('\n', fp);
            fprintf(fp, "# X START: %12.6E Y START: %12.6E Z START: %12.6E\n", x0, y0, z0);
            fprintf(fp,
                    "# X END  : %12.6E Y END  : %12.6E Z END  : %12.6E\n",
                    geo->axes[ix_axis].hi,
                    geo->axes[iy_axis].hi,
                    geo->axes[iz_axis].hi);
            fprintf(fp, "# PRIMARIES: %llu\n", nstat);
            fprintf(fp, "# COMPLETENESS: %s\n", osh_scoring_runtime_completeness_label(rt));
            fprintf(fp, "# Data written in canonical flat mesh order: idx = ix + nx * (iy + ny * iz)\n");
            fprintf(fp,
                    "# Values: NORM/SUM quantities divided by nstat; AVER quantities (DLET/TLET) written as physical "
                    "mean\n");
            if (diff_nbins > 0u) {
                fprintf(fp,
                        "# Diff1Type: %s  lo=%g  hi=%g  nbins=%zu%s%s\n",
                        diff_kind_label(p0->diff_kind),
                        p0->diff_lo,
                        p0->diff_hi,
                        p0->diff_nbins,
                        p0->diff_log ? " LOG" : "",
                        p0->has_diff_sset ? "  (SP override active)" : "");
            }
            if (diff2_nbins > 0u) {
                fprintf(fp,
                        "# Diff2Type: %s  lo=%g  hi=%g  nbins=%zu%s%s\n",
                        diff_kind_label(p0->diff2_kind),
                        p0->diff2_lo,
                        p0->diff2_hi,
                        p0->diff2_nbins,
                        p0->diff2_log ? " LOG" : "",
                        p0->has_diff2_sset ? "  (SP override active)" : "");
            }
            fprintf(fp, "# X Y Z");
            if (diff_nbins > 0u) {
                fprintf(fp, " %s", diff_kind_label(p0->diff_kind));
            }
            if (diff2_nbins > 0u) {
                fprintf(fp, " %s", diff_kind_label(p0->diff2_kind));
            }
            fprint_quantity_names(fp, rt, out);
            fputc('\n', fp);

            for (iz = 0; iz < nz; ++iz) {
                for (iy = 0; iy < ny; ++iy) {
                    for (ix = 0; ix < nx; ++ix) {
                        size_t spatial_idx = ix + nx * (iy + ny * iz);
                        size_t db;
                        size_t db2;
                        size_t ndb = (diff_nbins > 0u) ? diff_nbins : 1u;
                        size_t ndb2 = (diff2_nbins > 0u) ? diff2_nbins : 1u;
                        for (db = 0; db < ndb; ++db) {
                            for (db2 = 0; db2 < ndb2; ++db2) {
                                fprintf(fp,
                                        " %.12e %.12e %.12e",
                                        x0 + dx * ((double) ix + 0.5),
                                        y0 + dy * ((double) iy + 0.5),
                                        z0 + dz * ((double) iz + 0.5));
                                if (diff_nbins > 0u) {
                                    fprintf(fp, " %.12e", ascii_diff_center(p0, db));
                                }
                                if (diff2_nbins > 0u) {
                                    fprintf(fp, " %.12e", ascii_diff2_center(p0, db2));
                                }
                                for (ip = 0; ip < out->npages; ++ip) {
                                    size_t page_idx = out->page_indices[ip];
                                    struct osh_scoring_page_runtime const *page = &rt->pages[page_idx];
                                    /* AVER pages (e.g. DLET/TLET) hold a physical mean after
                                     * postprocessing — do not normalise per primary. */
                                    double scale = (page->postproc == OSH_SCORING_POSTPROC_AVER) ? 1.0 : inv_nstat;
                                    /* Use each page's own diff_nbins to gate the offset: non-differential
                                     * pages have len==geo_nbins and must not have db/db2 offsets applied,
                                     * even when other pages in the same output are differential. */
                                    size_t data_idx = spatial_idx
                                                      + (page->diff_nbins > 0u ? db * page->diff_stride : 0u)
                                                      + (page->diff2_nbins > 0u ? db2 * page->diff2_stride : 0u);
                                    fprintf(fp, " %.12e", page->acc.data[data_idx] * scale);
                                }
                                fprintf(fp, "\n");
                            }
                        }
                    }
                }
            }
        }
    }

    if (fclose(fp) != 0) {
        return OSH_EIO;
    }
    return OSH_OK;
}

static void
fprint_quantity_names(FILE *fp, struct osh_scoring_runtime const *rt, struct osh_scoring_output_runtime const *out) {
    size_t ip;
    size_t page_idx;
    char const *qty;
    char const *c;

    for (ip = 0; ip < out->npages; ++ip) {
        page_idx = out->page_indices[ip];
        qty = rt->pages[page_idx].quantity ? rt->pages[page_idx].quantity : "?";
        fputc(' ', fp);
        for (c = qty; *c; ++c) {
            fputc(toupper((unsigned char) *c), fp);
        }
    }
}

static enum osh_status
mesh_axis_index(struct osh_scoring_geometry_runtime const *geo, char const *label, size_t *idx_out) {
    size_t i;

    if (!geo || !label || !idx_out) {
        return OSH_EINVAL;
    }
    for (i = 0; i < geo->naxes; ++i) {
        if (strcmp(geo->axes[i].label, label) == 0) {
            *idx_out = i;
            return OSH_OK;
        }
    }
    return OSH_ENOTSUP;
}

static enum osh_status validate_output(struct osh_scoring_workspace const *ws,
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
    if (geo->geo_kind != OSH_SCORING_GEO_MESH && geo->geo_kind != OSH_SCORING_GEO_CYL
        && geo->geo_kind != OSH_SCORING_GEO_ZONE) {
        return OSH_ENOTSUP;
    }
    /* Rotated MESH ASCII output would emit local-frame X/Y/Z which is ambiguous;
     * reject until a proper coordinate annotation is implemented.
     * CYL R/Z output is always in local frame regardless of rotation, so it is allowed. */
    if (geo->geo_kind == OSH_SCORING_GEO_MESH && geo->has_rotation) {
        return OSH_ENOTSUP;
    }
    for (ip = 0; ip < out->npages; ++ip) {
        struct osh_scoring_page_runtime const *page = &rt->pages[out->page_indices[ip]];
        if (!page->acc.data || page->variance || page->has_data2 || page->divide) {
            return OSH_ENOTSUP;
        }
    }

    *out_out = out;
    *geo_out = geo;
    return OSH_OK;
}

static void format_now_rfc2822(char *buf, size_t cap) {
    time_t now;

    time(&now);
    strftime(buf, cap, "%a, %d %b %Y %H:%M:%S %z", localtime(&now));
}

/* Compute the centre value of differential bin db (linear: arithmetic midpoint;
 * log: geometric midpoint). */
static double ascii_diff_center(struct osh_scoring_page_runtime const *page, size_t db) {
    double t0; /* normalised left edge */
    double t1; /* normalised right edge */

    t0 = (double) db / (double) page->diff_nbins;
    t1 = (double) (db + 1u) / (double) page->diff_nbins;
    if (page->diff_log) {
        double ratio = page->diff_hi / page->diff_lo;
        return page->diff_lo * sqrt(pow(ratio, t0) * pow(ratio, t1));
    }
    return page->diff_lo + 0.5 * (t0 + t1) * (page->diff_hi - page->diff_lo);
}

/* Compute the centre value of the second differential bin db2. */
static double ascii_diff2_center(struct osh_scoring_page_runtime const *page, size_t db) {
    double t0;
    double t1;

    t0 = (double) db / (double) page->diff2_nbins;
    t1 = (double) (db + 1u) / (double) page->diff2_nbins;
    if (page->diff2_log) {
        double ratio = page->diff2_hi / page->diff2_lo;
        return page->diff2_lo * sqrt(pow(ratio, t0) * pow(ratio, t1));
    }
    return page->diff2_lo + 0.5 * (t0 + t1) * (page->diff2_hi - page->diff2_lo);
}

/* Map diff_kind enum to a short ASCII label used in column headers. */
static char const *diff_kind_label(enum osh_scoring_diff_kind kind) {
    switch (kind) {
    case OSH_SCORING_DIFF_EKIN:
        return "EKIN";
    case OSH_SCORING_DIFF_ENUC:
        return "ENUC";
    case OSH_SCORING_DIFF_EAMU:
        return "EAMU";
    case OSH_SCORING_DIFF_LET:
        return "LET";
    case OSH_SCORING_DIFF_QEFF:
        return "QEFF";
    default:
        return "?";
    }
}
