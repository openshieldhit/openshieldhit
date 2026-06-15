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

    if (geo->geo_kind == OSH_SCORING_GEO_CYL) {
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

        fprintf(fp, "# OpenShieldHIT version %s\n", OSH_VERSION);
        fprintf(fp, "# Calculated %s\n", datestr);
        fprintf(fp, "# DETECTOR OUTPUT CYL\n");
        fprintf(fp, "# R BIN: %5zu Z BIN: %5zu\n", nr, nz);
        fprintf(fp, "# DETECTOR TYPE:");
        for (ip = 0; ip < out->npages; ++ip) {
            size_t page_idx = out->page_indices[ip];
            char const *qty = rt->pages[page_idx].quantity ? rt->pages[page_idx].quantity : "?";
            char const *c;
            fputc(' ', fp);
            for (c = qty; *c; ++c) {
                fputc(toupper((unsigned char) *c), fp);
            }
        }
        fputc('\n', fp);
        fprintf(fp, "# R START: %12.6E Z START: %12.6E\n", r0, z0);
        fprintf(fp, "# R END  : %12.6E Z END  : %12.6E\n", geo->axes[ir_axis].hi, geo->axes[iz_axis].hi);
        fprintf(fp, "# PRIMARIES: %llu\n", nstat);
        fprintf(fp, "# Data written in canonical flat cyl order: idx = ir + nr * iz\n");
        fprintf(
            fp,
            "# Values: NORM/SUM quantities divided by nstat; AVER quantities (DLET/TLET) written as physical mean\n");
        fprintf(fp, "# R Z");
        for (ip = 0; ip < out->npages; ++ip) {
            size_t page_idx = out->page_indices[ip];
            char const *qty = rt->pages[page_idx].quantity ? rt->pages[page_idx].quantity : "?";
            char const *c;
            fputc(' ', fp);
            for (c = qty; *c; ++c) {
                fputc(toupper((unsigned char) *c), fp);
            }
        }
        fputc('\n', fp);

        for (iz = 0; iz < nz; ++iz) {
            for (ir = 0; ir < nr; ++ir) {
                size_t idx = ir + nr * iz;
                fprintf(fp, " %.12e %.12e", r0 + dr * ((double) ir + 0.5), z0 + dz * ((double) iz + 0.5));
                for (ip = 0; ip < out->npages; ++ip) {
                    size_t page_idx = out->page_indices[ip];
                    double scale = (rt->pages[page_idx].postproc == OSH_SCORING_POSTPROC_AVER) ? 1.0 : inv_nstat;
                    fprintf(fp, " %.12e", rt->pages[page_idx].data[idx] * scale);
                }
                fprintf(fp, "\n");
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

        fprintf(fp, "# OpenShieldHIT version %s\n", OSH_VERSION);
        fprintf(fp, "# Calculated %s\n", datestr);
        fprintf(fp, "# DETECTOR OUTPUT MSH\n");
        fprintf(fp, "# X BIN: %5zu Y BIN: %5zu Z BIN: %5zu\n", nx, ny, nz);
        fprintf(fp, "# DETECTOR TYPE:");
        for (ip = 0; ip < out->npages; ++ip) {
            size_t page_idx = out->page_indices[ip];
            char const *qty = rt->pages[page_idx].quantity ? rt->pages[page_idx].quantity : "?";
            char const *c;
            fputc(' ', fp);
            for (c = qty; *c; ++c) {
                fputc(toupper((unsigned char) *c), fp);
            }
        }
        fputc('\n', fp);
        fprintf(fp, "# X START: %12.6E Y START: %12.6E Z START: %12.6E\n", x0, y0, z0);
        fprintf(fp,
                "# X END  : %12.6E Y END  : %12.6E Z END  : %12.6E\n",
                geo->axes[ix_axis].hi,
                geo->axes[iy_axis].hi,
                geo->axes[iz_axis].hi);
        fprintf(fp, "# PRIMARIES: %llu\n", nstat);
        fprintf(fp, "# Data written in canonical flat mesh order: idx = ix + nx * (iy + ny * iz)\n");
        fprintf(
            fp,
            "# Values: NORM/SUM quantities divided by nstat; AVER quantities (DLET/TLET) written as physical mean\n");
        fprintf(fp, "# X Y Z");
        for (ip = 0; ip < out->npages; ++ip) {
            size_t page_idx = out->page_indices[ip];
            char const *qty = rt->pages[page_idx].quantity ? rt->pages[page_idx].quantity : "?";
            char const *c;
            fputc(' ', fp);
            for (c = qty; *c; ++c) {
                fputc(toupper((unsigned char) *c), fp);
            }
        }
        fputc('\n', fp);

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
    if ((geo->geo_kind != OSH_SCORING_GEO_MESH && geo->geo_kind != OSH_SCORING_GEO_CYL) || geo->has_rotation) {
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
