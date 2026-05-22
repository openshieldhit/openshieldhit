/*
 * BDO 2019 binary scorer output writer.
 *
 * Normalisation contract
 * ----------------------
 * Page data are written EXACTLY as accumulated in the runtime buffers (raw
 * un-normalised sums over all transport steps and all primary histories).
 * The total primary count is embedded as the OSHBDO_RT_NSTAT tag so that
 * readers and merge tools can apply the correct post-hoc normalisation:
 *
 *   NORM quantities (DOSE, FLUENCE, ENERGY, …):
 *       value_per_primary = data / nstat
 *
 *   AVER quantities (DLET, TLET, …):
 *       averaged_value = data / nstat     (data already holds nstat-weighted sum
 *                                          after osh_scoring_postprocess())
 *
 *   SUM quantities (COUNT, …):
 *       total = data          (no normalisation)
 *
 * Multi-run merging: when combining several BDO files from independent runs
 * (embarrassingly-parallel or multi-node setups) each file contributes its
 * own nstat, so the merge tool must weight accordingly.  The OSHBDO_PAG_NORMALIZE
 * tag records the postproc mode per page so a merge tool does not need to
 * re-derive it from the scorer kind.
 *
 * DOSE is stored in Gy (osh_scoring_postprocess() has already applied the
 * MeV/g → Gy conversion); DLET and TLET are stored in MeV/cm.
 */

#include "scoring/save/osh_scoring_save_bdo2019.h"

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "common/osh_version.h"
#include "scoring/save/osh_scoring_save_bdo2019_raw.h"

static char const *geometry_type_name(struct osh_scoring_geometry_runtime const *geo);
static char const *page_data_unit(struct osh_scoring_page_runtime const *page);
static enum osh_status validate_output(struct osh_scoring_workspace const *ws,
                                       struct osh_scoring_runtime const *rt,
                                       size_t output_idx,
                                       struct osh_scoring_output_runtime const **out_out,
                                       struct osh_scoring_geometry_runtime const **geo_out);
static void geometry_mesh_arrays(struct osh_scoring_geometry_runtime const *geo, double p[3], double q[3], int n[3]);
static void format_now_rfc2822(char *buf, size_t cap);
static int legacy_geo_kind(struct osh_scoring_geometry_runtime const *geo);
static int legacy_score_kind(struct osh_scoring_page_runtime const *page);

enum osh_status osh_scoring_save_bdo2019_output(struct osh_scoring_workspace const *ws,
                                                struct osh_scoring_runtime const *rt,
                                                unsigned long long nstat,
                                                size_t output_idx) {
    enum osh_status rc;
    FILE *fp;
    char datestr[128];
    struct osh_scoring_output_runtime const *out;
    struct osh_scoring_geometry_runtime const *geo;
    double p[3];
    double q[3];
    int n[3];
    int format_id;
    int est_count;
    int est_npages;
    double est_rescale_nstat;
    size_t ip;

    rc = validate_output(ws, rt, output_idx, &out, &geo);
    if (rc != OSH_OK) {
        return rc;
    }
    if (nstat == 0ull) {
        return OSH_EINVAL;
    }
    if (nstat > (unsigned long long) LLONG_MAX) {
        return OSH_EINVAL;
    }

    fp = fopen(out->filename, "wb");
    if (!fp) {
        return OSH_EIO;
    }

    rc = osh_scoring_bdo2019_write_preamble(fp, OSH_VERSION);
    if (rc != OSH_OK) {
        fclose(fp);
        return rc;
    }

    format_now_rfc2822(datestr, sizeof(datestr));
    format_id = OSH_SCORING_BDO2019_FORMAT_ID;
    est_count = (int) output_idx;
    est_npages = (int) out->npages;
    est_rescale_nstat = 1.0;
    geometry_mesh_arrays(geo, p, q, n);

    rc = osh_scoring_bdo2019_write_token_str(fp, OSHBDO_SHVERSION, OSH_VERSION);
    if (rc == OSH_OK) {
        rc = osh_scoring_bdo2019_write_token_int(fp, OSHBDO_FORMAT, &format_id, 1u);
    }
    if (rc == OSH_OK) {
        rc = osh_scoring_bdo2019_write_token_str(fp, OSHBDO_FILEDATE, datestr);
    }
    if (rc == OSH_OK) {
        long long int nstat_ll = (long long int) nstat;
        rc = osh_scoring_bdo2019_write_token_llint(fp, OSHBDO_RT_NSTAT, &nstat_ll, 1u);
    }
    if (rc == OSH_OK) {
        rc = osh_scoring_bdo2019_write_token_str(fp, OSHBDO_GEO_TYPE, geometry_type_name(geo));
    }
    if (rc == OSH_OK) {
        rc = osh_scoring_bdo2019_write_token_str(fp, OSHBDO_GEO_NAME, geo->name ? geo->name : "(unnamed)");
    }
    if (rc == OSH_OK) {
        rc = osh_scoring_bdo2019_write_token_double(fp, OSHBDO_GEO_P, p, 3u);
    }
    if (rc == OSH_OK) {
        rc = osh_scoring_bdo2019_write_token_double(fp, OSHBDO_GEO_Q, q, 3u);
    }
    if (rc == OSH_OK) {
        rc = osh_scoring_bdo2019_write_token_int(fp, OSHBDO_GEO_N, n, 3u);
    }
    if (rc == OSH_OK) {
        rc = osh_scoring_bdo2019_write_token_str(fp, OSHBDO_EST_FILENAME, out->filename);
    }
    if (rc == OSH_OK) {
        rc = osh_scoring_bdo2019_write_token_int(fp, OSHBDO_EST_COUNT, &est_count, 1u);
    }
    if (rc == OSH_OK) {
        rc = osh_scoring_bdo2019_write_token_int(fp, OSHBDO_EST_NPAGES, &est_npages, 1u);
    }
    if (rc == OSH_OK) {
        rc = osh_scoring_bdo2019_write_token_double(fp, OSHBDO_EST_RESCALE_NSTAT, &est_rescale_nstat, 1u);
    }

    for (ip = 0; rc == OSH_OK && ip < out->npages; ++ip) {
        int page_type;
        int page_count;
        int page_norm;
        double page_rescale;
        double page_offset;
        struct osh_scoring_page_runtime const *page = &rt->pages[out->page_indices[ip]];

        page_type = legacy_score_kind(page);
        page_count = (int) ip;
        page_norm = (int) page->postproc;
        page_rescale = 1.0;
        page_offset = 0.0;

        rc = osh_scoring_bdo2019_write_token_int(fp, OSHBDO_PAG_TYPE, &page_type, 1u);
        if (rc == OSH_OK) {
            rc = osh_scoring_bdo2019_write_token_int(fp, OSHBDO_PAG_COUNT, &page_count, 1u);
        }
        if (rc == OSH_OK) {
            rc = osh_scoring_bdo2019_write_token_int(fp, OSHBDO_PAG_NORMALIZE, &page_norm, 1u);
        }
        if (rc == OSH_OK) {
            rc = osh_scoring_bdo2019_write_token_double(fp, OSHBDO_PAG_RESCALE, &page_rescale, 1u);
        }
        if (rc == OSH_OK) {
            rc = osh_scoring_bdo2019_write_token_double(fp, OSHBDO_PAG_OFFSET, &page_offset, 1u);
        }
        if (rc == OSH_OK) {
            rc = osh_scoring_bdo2019_write_token_str(fp, OSHBDO_PAG_DATA_UNIT, page_data_unit(page));
        }
        if (rc == OSH_OK) {
            /* TODO: current runtime does not yet carry full normalization
             * metadata into save. Data are written exactly as stored in the
             * runtime page buffer, in canonical flat order. */
            rc = osh_scoring_bdo2019_write_token_double(fp, OSHBDO_PAG_DATA, page->data, page->len);
        }
    }

    if (fclose(fp) != 0 && rc == OSH_OK) {
        return OSH_EIO;
    }
    return rc;
}

static char const *geometry_type_name(struct osh_scoring_geometry_runtime const *geo) {
    switch (legacy_geo_kind(geo)) {
    case 1:
        return "MSH";
    case 2:
        return "CYL";
    case 3:
        return "ZONE";
    case 5:
        return "VOXSCORE";
    default:
        return "NONE";
    }
}

static char const *page_data_unit(struct osh_scoring_page_runtime const *page) {
    switch (page->score_kind) {
    case OSH_SCORING_SCORE_ENERGY:
        return "MeV";
    case OSH_SCORING_SCORE_FLUENCE:
        return "1/cm^2";
    case OSH_SCORING_SCORE_DOSE:
        return "Gy";
    case OSH_SCORING_SCORE_DLET:
    case OSH_SCORING_SCORE_TLET:
        return "MeV/cm";
    default:
        return "arb";
    }
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
    if (output_idx >= ws->noutputs || output_idx >= rt->noutputs) {
        return OSH_EINVAL;
    }

    out = &rt->outputs[output_idx];
    if (out->geometry_idx >= rt->ngeometries) {
        return OSH_ESTATE;
    }
    geo = &rt->geometries[out->geometry_idx];
    if (geo->geo_kind != OSH_SCORING_GEO_MESH) {
        return OSH_ENOTSUP;
    }
    for (ip = 0; ip < out->npages; ++ip) {
        struct osh_scoring_page_runtime const *page = &rt->pages[out->page_indices[ip]];
        if (!page->data || page->variance || page->has_data2 || page->divide) {
            return OSH_ENOTSUP;
        }
        if (legacy_score_kind(page) < 0) {
            return OSH_ENOTSUP;
        }
    }

    *out_out = out;
    *geo_out = geo;
    return OSH_OK;
}

static void geometry_mesh_arrays(struct osh_scoring_geometry_runtime const *geo, double p[3], double q[3], int n[3]) {
    size_t i;

    p[0] = p[1] = p[2] = 0.0;
    q[0] = q[1] = q[2] = 0.0;
    n[0] = n[1] = n[2] = 1;

    for (i = 0; i < geo->naxes; ++i) {
        if (strcmp(geo->axes[i].label, "X") == 0) {
            p[0] = geo->axes[i].lo;
            q[0] = geo->axes[i].hi;
            n[0] = geo->axes[i].nbins;
        } else if (strcmp(geo->axes[i].label, "Y") == 0) {
            p[1] = geo->axes[i].lo;
            q[1] = geo->axes[i].hi;
            n[1] = geo->axes[i].nbins;
        } else if (strcmp(geo->axes[i].label, "Z") == 0) {
            p[2] = geo->axes[i].lo;
            q[2] = geo->axes[i].hi;
            n[2] = geo->axes[i].nbins;
        }
    }
}

static void format_now_rfc2822(char *buf, size_t cap) {
    time_t now;

    time(&now);
    strftime(buf, cap, "%a, %d %b %Y %H:%M:%S %z", localtime(&now));
}

static int legacy_geo_kind(struct osh_scoring_geometry_runtime const *geo) {
    switch (geo->geo_kind) {
    case OSH_SCORING_GEO_MESH:
        return 1;
    case OSH_SCORING_GEO_CYL:
        return 2;
    case OSH_SCORING_GEO_ZONE:
        return 3;
    default:
        return 0;
    }
}

static int legacy_score_kind(struct osh_scoring_page_runtime const *page) {
    switch (page->score_kind) {
    case OSH_SCORING_SCORE_ENERGY:
        return 1;
    case OSH_SCORING_SCORE_FLUENCE:
        return 2;
    case OSH_SCORING_SCORE_DOSE:
        return 5;
    case OSH_SCORING_SCORE_LETFLU:
        return 4;
    case OSH_SCORING_SCORE_DLET:
        return 6;
    case OSH_SCORING_SCORE_TLET:
        return 7;
    case OSH_SCORING_SCORE_NORMCOUNT:
        return 14;
    case OSH_SCORING_SCORE_COUNT:
        return 55;
    case OSH_SCORING_SCORE_NKERMA:
        return 38;
    case OSH_SCORING_SCORE_ALANINE:
        return 13;
    case OSH_SCORING_SCORE_MCPL:
        return 63;
    default:
        return -1;
    }
}
