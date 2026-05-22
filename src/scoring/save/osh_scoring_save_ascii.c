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

enum osh_status osh_scoring_save_ascii_output(struct osh_scoring_workspace const *ws,
                                              struct osh_scoring_runtime const *rt,
                                              unsigned long long nstat,
                                              size_t output_idx) {
    FILE *fp;
    char datestr[128];
    enum osh_status rc;
    struct osh_scoring_output_runtime const *out;
    struct osh_scoring_geometry_runtime const *geo;
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
    double inv_nstat;
    size_t ip;
    size_t ix;
    size_t iy;
    size_t iz;

    rc = validate_output(ws, rt, output_idx, &out, &geo);
    if (rc != OSH_OK) {
        return rc;
    }
    if (nstat == 0ull) {
        return OSH_EINVAL;
    }
    rc = mesh_axis_index(geo, "X", &ix_axis);
    if (rc != OSH_OK) {
        return rc;
    }
    rc = mesh_axis_index(geo, "Y", &iy_axis);
    if (rc != OSH_OK) {
        return rc;
    }
    rc = mesh_axis_index(geo, "Z", &iz_axis);
    if (rc != OSH_OK) {
        return rc;
    }

    fp = fopen(out->filename, "w");
    if (!fp) {
        return OSH_EIO;
    }

    inv_nstat = 1.0 / (double) nstat;

    format_now_rfc2822(datestr, sizeof(datestr));
    fprintf(fp, "# OpenShieldHIT version %s\n", OSH_VERSION);
    fprintf(fp, "# Calculated %s\n", datestr);
    fprintf(fp, "# Geometry %s\n", geo->name ? geo->name : "(unnamed)");
    fprintf(fp, "# nstat %llu\n", nstat);
    fprintf(fp, "# Data are written in canonical flat mesh order: idx = ix + nx * (iy + ny * iz)\n");
    fprintf(fp, "# Values: NORM/SUM pages divided by nstat; AVER pages (DLET/TLET/…) written as physical mean\n");
    fprintf(fp, "# X Y Z");
    for (ip = 0; ip < out->npages; ++ip) {
        size_t page_idx = out->page_indices[ip];
        fprintf(fp, " Page(%zu:%s)", ip, rt->pages[page_idx].quantity ? rt->pages[page_idx].quantity : "?");
    }
    fprintf(fp, "\n");

    nx = (size_t) geo->axes[ix_axis].nbins;
    ny = (size_t) geo->axes[iy_axis].nbins;
    nz = (size_t) geo->axes[iz_axis].nbins;
    x0 = geo->axes[ix_axis].lo;
    y0 = geo->axes[iy_axis].lo;
    z0 = geo->axes[iz_axis].lo;
    dx = (geo->axes[ix_axis].hi - geo->axes[ix_axis].lo) / (double) geo->axes[ix_axis].nbins;
    dy = (geo->axes[iy_axis].hi - geo->axes[iy_axis].lo) / (double) geo->axes[iy_axis].nbins;
    dz = (geo->axes[iz_axis].hi - geo->axes[iz_axis].lo) / (double) geo->axes[iz_axis].nbins;

    for (iz = 0; iz < nz; ++iz) {
        for (iy = 0; iy < ny; ++iy) {
            for (ix = 0; ix < nx; ++ix) {
                size_t idx = ix + nx * (iy + ny * iz);

                fprintf(fp,
                        " %.12e %.12e %.12e",
                        x0 + dx * ((double) ix + 0.5),
                        y0 + dy * ((double) iy + 0.5),
                        z0 + dz * ((double) iz + 0.5));
                for (ip = 0; ip < out->npages; ++ip) {
                    size_t page_idx = out->page_indices[ip];
                    /* AVER pages (e.g. DLET/TLET) hold a physical mean after
                     * postprocessing — do not normalise per primary. */
                    double scale = (rt->pages[page_idx].postproc == OSH_SCORING_POSTPROC_AVER) ? 1.0 : inv_nstat;
                    fprintf(fp, " %.12e", rt->pages[page_idx].data[idx] * scale);
                }
                fprintf(fp, "\n");
            }
        }
    }

    if (fclose(fp) != 0) {
        return OSH_EIO;
    }
    return OSH_OK;
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
    if (geo->geo_kind != OSH_SCORING_GEO_MESH || geo->has_rotation) {
        return OSH_ENOTSUP;
    }
    for (ip = 0; ip < out->npages; ++ip) {
        struct osh_scoring_page_runtime const *page = &rt->pages[out->page_indices[ip]];
        if (!page->data || page->variance || page->has_data2 || page->divide) {
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
