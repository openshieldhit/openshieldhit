#include "scoring/runtime/osh_scoring_compile.h"

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_diag.h"
#include "common/raytrace/osh_raytrace.h"
#include "openshieldhit/const.h"
#include "openshieldhit/scoring.h"
#include "scoring/runtime/osh_scoring_accumulator.h"
#include "scoring/runtime/osh_scoring_multiformat.h"
#include "scoring/runtime/osh_scoring_postprocess.h"

/* Scratch record built during the first compile pass; sorted by geometry+kind
 * so pages that share a geometry and scorer type end up in the same hot group. */
struct prepared_page_ref {
    size_t output_idx;
    size_t src_page_idx;
    size_t ordinal;
    size_t geometry_idx;
    enum osh_scoring_score_kind score_kind;
};

static enum osh_scoring_geo_kind geometry_kind_to_enum(char const *kind) {
    if (!kind) {
        return OSH_SCORING_GEO_UNKNOWN;
    }
    if (strcmp(kind, "mesh") == 0) {
        return OSH_SCORING_GEO_MESH;
    }
    if (strcmp(kind, "cyl") == 0) {
        return OSH_SCORING_GEO_CYL;
    }
    if (strcmp(kind, "zone") == 0) {
        return OSH_SCORING_GEO_ZONE;
    }
    /* For voxelised Cartesian structures (e.g. from DICOM), set kind to "mesh". */
    if (strcmp(kind, "all") == 0) {
        return OSH_SCORING_GEO_ALL;
    }
    return OSH_SCORING_GEO_UNKNOWN;
}

static size_t geometry_nbins(struct osh_scoring_geometry_def const *geo) {
    size_t i;
    size_t nbins = 1u;

    if (geo->naxes > 0u) {
        for (i = 0; i < geo->naxes; ++i) {
            if (geo->axes[i].nbins <= 0) {
                return 0u;
            }
            nbins *= (size_t) geo->axes[i].nbins;
        }
        return nbins;
    }

    if (geo->nzone_indices > 0u) {
        return geo->nzone_indices;
    }

    return 0u;
}

static enum osh_scoring_score_kind quantity_to_score_kind(char const *quantity) {
    if (!quantity) {
        return OSH_SCORING_SCORE_UNKNOWN;
    }
    if (strcmp(quantity, "energy") == 0) {
        return OSH_SCORING_SCORE_ENERGY;
    }
    if (strcmp(quantity, "fluence") == 0) {
        return OSH_SCORING_SCORE_FLUENCE;
    }
    if (strcmp(quantity, "dose") == 0) {
        return OSH_SCORING_SCORE_DOSE;
    }
    if (strcmp(quantity, "dosegy") == 0) {
        return OSH_SCORING_SCORE_DOSEGY;
    }
    if (strcmp(quantity, "dirtydose") == 0) {
        return OSH_SCORING_SCORE_DIRTYDOSE;
    }
    if (strcmp(quantity, "dirtydosegy") == 0) {
        return OSH_SCORING_SCORE_DIRTYDOSEGY;
    }
    if (strcmp(quantity, "dlet") == 0) {
        return OSH_SCORING_SCORE_DLET;
    }
    if (strcmp(quantity, "tlet") == 0) {
        return OSH_SCORING_SCORE_TLET;
    }
    if (strcmp(quantity, "dqeff") == 0) {
        return OSH_SCORING_SCORE_DQEFF;
    }
    if (strcmp(quantity, "tqeff") == 0) {
        return OSH_SCORING_SCORE_TQEFF;
    }
    if (strcmp(quantity, "davge") == 0) {
        return OSH_SCORING_SCORE_DAVGE;
    }
    if (strcmp(quantity, "tavge") == 0) {
        return OSH_SCORING_SCORE_TAVGE;
    }
    if (strcmp(quantity, "dbeta") == 0) {
        return OSH_SCORING_SCORE_DBETA;
    }
    if (strcmp(quantity, "tbeta") == 0) {
        return OSH_SCORING_SCORE_TBETA;
    }
    return OSH_SCORING_SCORE_UNKNOWN;
}

static enum osh_scoring_filter_field filter_field_to_enum(char const *field) {
    if (!field) {
        return OSH_SCORING_FILTER_FIELD_UNKNOWN;
    }
    if (strcmp(field, "ID") == 0) {
        return OSH_SCORING_FILTER_FIELD_ID;
    }
    if (strcmp(field, "Z") == 0) {
        return OSH_SCORING_FILTER_FIELD_Z;
    }
    if (strcmp(field, "A") == 0) {
        return OSH_SCORING_FILTER_FIELD_A;
    }
    if (strcmp(field, "AMASS") == 0) {
        return OSH_SCORING_FILTER_FIELD_AMASS;
    }
    if (strcmp(field, "AMU") == 0) {
        return OSH_SCORING_FILTER_FIELD_AMU;
    }
    if (strcmp(field, "E") == 0) {
        return OSH_SCORING_FILTER_FIELD_E;
    }
    if (strcmp(field, "ENUC") == 0) {
        return OSH_SCORING_FILTER_FIELD_ENUC;
    }
    if (strcmp(field, "EAMU") == 0) {
        return OSH_SCORING_FILTER_FIELD_EAMU;
    }
    if (strcmp(field, "GEN") == 0) {
        return OSH_SCORING_FILTER_FIELD_GEN;
    }
    if (strcmp(field, "NPRIM") == 0) {
        return OSH_SCORING_FILTER_FIELD_NPRIM;
    }
    return OSH_SCORING_FILTER_FIELD_UNKNOWN;
}

static enum osh_scoring_filter_op filter_op_to_enum(char const *op) {
    if (!op || op[0] == '\0') {
        return OSH_SCORING_FILTER_OP_INVALID;
    }

    switch (op[0]) {
    case '<':
        if (op[1] == '\0') {
            return OSH_SCORING_FILTER_OP_LT;
        }
        if (op[1] == '=' && op[2] == '\0') {
            return OSH_SCORING_FILTER_OP_LE;
        }
        break;
    case '>':
        if (op[1] == '\0') {
            return OSH_SCORING_FILTER_OP_GT;
        }
        if (op[1] == '=' && op[2] == '\0') {
            return OSH_SCORING_FILTER_OP_GE;
        }
        break;
    case '=':
        if (op[1] == '\0') {
            return OSH_SCORING_FILTER_OP_EQ;
        }
        if (op[1] == '=' && op[2] == '\0') {
            return OSH_SCORING_FILTER_OP_EQ;
        }
        break;
    case '!':
        if (op[1] == '=' && op[2] == '\0') {
            return OSH_SCORING_FILTER_OP_NE;
        }
        break;
    default:
        break;
    }

    return OSH_SCORING_FILTER_OP_INVALID;
}

/* LET scorers use a two-pass accumulator: data = weighted sum, data2 = weight sum.
 * osh_scoring_postprocess() finalises the ratio data/data2 in-place. */
static char score_kind_uses_data2(enum osh_scoring_score_kind score_kind) {
    switch (score_kind) {
    case OSH_SCORING_SCORE_DLET:
    case OSH_SCORING_SCORE_TLET:
    case OSH_SCORING_SCORE_DQEFF:
    case OSH_SCORING_SCORE_TQEFF:
    case OSH_SCORING_SCORE_DAVGE:
    case OSH_SCORING_SCORE_TAVGE:
    case OSH_SCORING_SCORE_DBETA:
    case OSH_SCORING_SCORE_TBETA:
        return 1;
    default:
        return 0;
    }
}

/* Does any Settings block referenced on this page enable variance tracking?
 * Mirrors the per-page resolution in osh_scoring_compile() so the memory
 * estimate cannot drift from the real allocation.  Settings references share the
 * page's filter_names[] list (they are classified filter-vs-settings at compile). */
static int page_settings_enable_variance(struct osh_scoring_workspace const *ws,
                                         struct osh_scoring_page_def const *page) {
    size_t s;
    size_t t;

    for (s = 0u; s < page->nfilter_names; ++s) {
        for (t = 0u; t < ws->nsettings; ++t) {
            if (ws->settings[t].name && strcmp(ws->settings[t].name, page->filter_names[s]) == 0
                && ws->settings[t].has_variance && ws->settings[t].variance) {
                return 1;
            }
        }
    }
    return 0;
}

/* divide=1 means postproc divides data by data2 (LET average) -> AVER mode.
 * All other scorers default to NORM (÷nstat) except raw-count kinds. */
static enum osh_scoring_postproc score_kind_postproc(enum osh_scoring_score_kind score_kind, char divide) {
    if (divide) {
        return OSH_SCORING_POSTPROC_AVER;
    }

    switch (score_kind) {
    case OSH_SCORING_SCORE_COUNT:
        return OSH_SCORING_POSTPROC_SUM;
    case OSH_SCORING_SCORE_MCPL:
        return OSH_SCORING_POSTPROC_APPEND;
    default:
        return OSH_SCORING_POSTPROC_NORM;
    }
}

enum osh_status osh_scoring_estimate_memory(struct osh_scoring_workspace const *ws,
                                            struct osh_scoring_mem_estimate *out) {
    size_t i;
    size_t j;

    if (!ws || !out) {
        return OSH_EINVAL;
    }

    out->accum_bytes = 0u;
    out->shadow_bytes = 0u;
    out->npages = 0u;
    out->largest_page_bytes = 0u;
    out->largest_geometry[0] = '\0';

    /* Mirror exactly what osh_scoring_compile() allocates: one page per
     * (output, quantity); each page owns `bins` doubles for its primary
     * accumulator, plus a second `bins`-double weight accumulator for the
     * "average" quantities (LET/Qeff).  A per-estimator "Variance On" Settings
     * (issue #209) doubles this: each sum array gains a companion Welford M2 array.
     * bins is the geometry's bin product, computed by the same geometry_nbins() the
     * compiler uses, so the estimate cannot drift from the real allocation. */
    for (i = 0u; i < ws->noutputs; ++i) {
        struct osh_scoring_output_def const *output = &ws->outputs[i];
        struct osh_scoring_geometry_def const *geo = osh_scoring_geometry_by_name(ws, output->geometry_name);
        size_t const bins = geo ? geometry_nbins(geo) : 0u;

        for (j = 0u; j < output->npages; ++j) {
            struct osh_scoring_page_def const *page = &output->pages[j];
            enum osh_scoring_score_kind const kind = quantity_to_score_kind(page->quantity);
            char const has_data2 = (char) score_kind_uses_data2(kind);
            struct osh_scoring_page_runtime page_rt;
            unsigned sum_arrays = 1u; /* # of sum arrays: data, plus data2 for two-pass averages */
            unsigned arrays;          /* sum arrays, plus their M2 companions when variance is on */
            uint64_t len = (uint64_t) bins;
            uint64_t page_bytes;

            if (has_data2) {
                sum_arrays = 2u;
            }
            /* Variance doubles the arrays: data_var mirrors data, data2_var mirrors data2. */
            arrays = sum_arrays;
            if (page_settings_enable_variance(ws, page)) {
                arrays = sum_arrays * 2u;
            }

            if (len > 0u && page->diff_nbins > 0u) {
                uint64_t const d1 = (uint64_t) page->diff_nbins;
                uint64_t const d2 = (uint64_t) ((page->diff2_nbins > 0u) ? page->diff2_nbins : 1u);
                if (d1 > 0u) {
                    len = (len <= UINT64_MAX / d1) ? (len * d1) : UINT64_MAX;
                }
                if (d2 > 0u) {
                    len = (len <= UINT64_MAX / d2) ? (len * d2) : UINT64_MAX;
                }
            }

            {
                /* Size this page exactly as osh_scoring_accumulator_alloc() would,
                 * so the estimate cannot drift from the real allocation (an
                 * invariant the header documents). */
                uint64_t n = len;
                /* alloc never makes a zero-length array: a degenerate (bins==0)
                 * page still gets one element per array, so round len up to one. */
                if (n == 0u) {
                    n = 1u;
                }
                uint64_t const bytes_per_bin = (uint64_t) sizeof(double) * (uint64_t) arrays;
                page_bytes = (n <= UINT64_MAX / bytes_per_bin) ? (n * bytes_per_bin) : UINT64_MAX;
            }

            out->accum_bytes =
                (out->accum_bytes <= UINT64_MAX - page_bytes) ? (out->accum_bytes + page_bytes) : UINT64_MAX;
            out->npages += 1u;

            /* A mid-run snapshot copies only the `data` array (one, never data2)
             * of pages whose postprocess writes data — including page-level
             * differential additive pages such as Energy vs Ekin — at the same
             * per-page bin count used above. */
            memset(&page_rt, 0, sizeof(page_rt));
            page_rt.score_kind = kind;
            page_rt.postproc = score_kind_postproc(kind, has_data2);
            page_rt.diff_nbins = page->diff_nbins;
            page_rt.diff2_nbins = page->diff2_nbins;
            if (osh_scoring_postprocess_page_writes_data(&page_rt)) {
                uint64_t shadow_n = len;
                uint64_t shadow_page;
                /* osh_scoring_shadow_refresh() grabs one double even for a len==0
                 * page, so round up to keep the estimate matched to it. */
                if (shadow_n == 0u) {
                    shadow_n = 1u;
                }
                /* Bytes for this page's single shadow data array.  Guard the
                 * multiply: saturate at UINT64_MAX rather than wrap on overflow. */
                if (shadow_n <= UINT64_MAX / (uint64_t) sizeof(double)) {
                    shadow_page = shadow_n * (uint64_t) sizeof(double);
                } else {
                    shadow_page = UINT64_MAX;
                }
                /* Fold into the running total, again saturating instead of
                 * wrapping if the addition would overflow. */
                if (out->shadow_bytes <= UINT64_MAX - shadow_page) {
                    out->shadow_bytes += shadow_page;
                } else {
                    out->shadow_bytes = UINT64_MAX;
                }
            }

            if (page_bytes > out->largest_page_bytes) {
                char const *name = (geo && geo->name)      ? geo->name
                                   : output->geometry_name ? output->geometry_name
                                                           : "(unnamed)";
                out->largest_page_bytes = page_bytes;
                (void) snprintf(out->largest_geometry, sizeof(out->largest_geometry), "%s", name);
            }
        }
    }

    return OSH_OK;
}

static enum osh_scoring_diff_kind diff_kind_from_str(char const *s) {
    if (!s || strcmp(s, "ekin") == 0 || strcmp(s, "e") == 0) {
        return OSH_SCORING_DIFF_EKIN;
    }
    if (strcmp(s, "enuc") == 0) {
        return OSH_SCORING_DIFF_ENUC;
    }
    if (strcmp(s, "eamu") == 0) {
        return OSH_SCORING_DIFF_EAMU;
    }
    if (strcmp(s, "let") == 0 || strcmp(s, "dedx") == 0) {
        /* LET/DEDX is a differential-axis value (Diff1Type/Diff2Type), not a
         * standalone score quantity.  The averaged score quantities are DLET
         * and TLET. */
        return OSH_SCORING_DIFF_LET;
    }
    if (strcmp(s, "qeff") == 0 || strcmp(s, "zeff2beta2") == 0) {
        /* Same distinction as LET: QEFF is an axis value; DQEFF/TQEFF are the
         * averaged scored quantities. */
        return OSH_SCORING_DIFF_QEFF;
    }
    return OSH_SCORING_DIFF_NONE;
}

static int runtime_supports_geometry(struct osh_scoring_geometry_runtime const *geo) {
    if (!geo) {
        return 0;
    }
    return geo->geo_kind == OSH_SCORING_GEO_MESH || geo->geo_kind == OSH_SCORING_GEO_CYL
           || geo->geo_kind == OSH_SCORING_GEO_ZONE;
}

static int runtime_supports_score_kind(enum osh_scoring_score_kind score_kind) {
    switch (score_kind) {
    case OSH_SCORING_SCORE_ENERGY:
    case OSH_SCORING_SCORE_FLUENCE:
    case OSH_SCORING_SCORE_DOSE:
    case OSH_SCORING_SCORE_DOSEGY:
    case OSH_SCORING_SCORE_DIRTYDOSE:
    case OSH_SCORING_SCORE_DIRTYDOSEGY:
    case OSH_SCORING_SCORE_DLET:
    case OSH_SCORING_SCORE_TLET:
    case OSH_SCORING_SCORE_DQEFF:
    case OSH_SCORING_SCORE_TQEFF:
    case OSH_SCORING_SCORE_DAVGE:
    case OSH_SCORING_SCORE_TAVGE:
    case OSH_SCORING_SCORE_DBETA:
    case OSH_SCORING_SCORE_TBETA:
        return 1;
    default:
        return 0;
    }
}

/* Sort key: geometry first (locality), then score_kind (groups same-kind pages
 * together so the hot path can iterate a contiguous run), then original ordinal
 * to preserve user-specified page order within a group. */
static int compare_prepared_pages(void const *a, void const *b) {
    struct prepared_page_ref const *pa;
    struct prepared_page_ref const *pb;

    pa = (struct prepared_page_ref const *) a;
    pb = (struct prepared_page_ref const *) b;

    if (pa->geometry_idx < pb->geometry_idx) {
        return -1;
    }
    if (pa->geometry_idx > pb->geometry_idx) {
        return 1;
    }
    if ((int) pa->score_kind < (int) pb->score_kind) {
        return -1;
    }
    if ((int) pa->score_kind > (int) pb->score_kind) {
        return 1;
    }
    if (pa->ordinal < pb->ordinal) {
        return -1;
    }
    if (pa->ordinal > pb->ordinal) {
        return 1;
    }
    return 0;
}

static enum osh_status copy_filter_runtime(struct osh_scoring_filter_runtime *dst,
                                           struct osh_diag_sink const *diag,
                                           struct osh_scoring_filter_def const *src) {
    size_t i;
    enum osh_scoring_filter_field field;
    enum osh_scoring_filter_op op;

    memset(dst, 0, sizeof(*dst));
    dst->name = strdup(src->name);
    if (!dst->name) {
        return OSH_ENOMEM;
    }
    if (src->nrules == 0u) {
        return OSH_OK;
    }

    dst->rules = (struct osh_scoring_filter_runtime_rule *) calloc(src->nrules, sizeof(*dst->rules));
    if (!dst->rules) {
        return OSH_ENOMEM;
    }
    dst->nrules = src->nrules;
    for (i = 0; i < src->nrules; ++i) {
        field = filter_field_to_enum(src->rules[i].field);
        op = filter_op_to_enum(src->rules[i].op);
        if (field == OSH_SCORING_FILTER_FIELD_UNKNOWN) {
            OSH_DIAG_ERRORF(diag,
                            "Scoring filter '%s' uses unknown field '%s'",
                            src->name ? src->name : "(unnamed)",
                            src->rules[i].field);
            return OSH_EINVAL;
        }
        if (op == OSH_SCORING_FILTER_OP_INVALID) {
            OSH_DIAG_ERRORF(diag,
                            "Scoring filter '%s' uses unknown operator '%s'",
                            src->name ? src->name : "(unnamed)",
                            src->rules[i].op);
            return OSH_EINVAL;
        }
        dst->rules[i].field = field;
        dst->rules[i].op = op;
        dst->rules[i].value = src->rules[i].value;
    }
    return OSH_OK;
}

static enum osh_status copy_settings_runtime(struct osh_scoring_settings_runtime *dst,
                                             struct osh_scoring_settings_def const *src) {
    memset(dst, 0, sizeof(*dst));
    dst->name = strdup(src->name);
    if (!dst->name) {
        return OSH_ENOMEM;
    }
    dst->rescale = src->rescale;
    dst->offset = src->offset;
    dst->site_diameter_um = src->site_diameter_um;
    dst->density_g_cm3 = src->density_g_cm3;
    dst->npart = src->npart;
    dst->medium = src->medium;
    dst->nkmedium = src->nkmedium;
    dst->variance = src->variance;
    dst->has_rescale = src->has_rescale;
    dst->has_offset = src->has_offset;
    dst->has_site_diameter_um = src->has_site_diameter_um;
    dst->has_density_g_cm3 = src->has_density_g_cm3;
    dst->has_npart = src->has_npart;
    dst->has_medium = src->has_medium;
    dst->has_nkmedium = src->has_nkmedium;
    dst->has_variance = src->has_variance;
    return OSH_OK;
}

static enum osh_status copy_geometry_runtime(struct osh_scoring_geometry_runtime *dst,
                                             struct osh_scoring_geometry_def const *src) {
    size_t i;

    memset(dst, 0, sizeof(*dst));
    dst->kind = strdup(src->kind);
    if (!dst->kind) {
        return OSH_ENOMEM;
    }
    dst->name = strdup(src->name);
    if (!dst->name) {
        return OSH_ENOMEM;
    }
    if (src->naxes > 0u) {
        dst->axes = (struct osh_scoring_axis_runtime *) calloc(src->naxes, sizeof(*dst->axes));
        if (!dst->axes) {
            return OSH_ENOMEM;
        }
        dst->naxes = src->naxes;
        for (i = 0; i < src->naxes; ++i) {
            memcpy(dst->axes[i].label, src->axes[i].label, sizeof(dst->axes[i].label));
            dst->axes[i].lo = src->axes[i].lo;
            dst->axes[i].hi = src->axes[i].hi;
            dst->axes[i].nbins = src->axes[i].nbins;
        }
    }
    memcpy(dst->t, src->t, sizeof(dst->t));
    if (src->nzone_indices > 0u) {
        if (!src->zone_indices) {
            /* Zone names were parsed but never resolved to transport indices.
             * osh_scoring_resolve_zone_names() must run (at app level, with the
             * geometry table) before compile; fail loudly rather than deref NULL. */
            return OSH_ESTATE;
        }
        dst->zone_indices = (size_t *) calloc(src->nzone_indices, sizeof(*dst->zone_indices));
        if (!dst->zone_indices) {
            return OSH_ENOMEM;
        }
        dst->nzone_indices = src->nzone_indices;
        for (i = 0; i < src->nzone_indices; ++i) {
            dst->zone_indices[i] = src->zone_indices[i];
        }
    }
    dst->has_rotation = src->has_rotation;
    dst->geo_kind = geometry_kind_to_enum(src->kind);
    dst->nbins = geometry_nbins(src);
    if (dst->nbins == 0u) {
        return OSH_EINVAL;
    }
    if (dst->geo_kind == OSH_SCORING_GEO_ZONE) {
        /* Zone bin volume comes from the user's per-zone Volume, so build the sole
         * volume source (bin_vol_inv) here where zone_volumes is in scope; Mesh/Cyl
         * bin_vol_inv is built from axes in osh_scoring_compile.  A missing/zero
         * Volume defaults to 1.0 cm3 (osh_scoring_resolve_zone_names already warns). */
        dst->bin_vol_inv = (double *) malloc(dst->nbins * sizeof(*dst->bin_vol_inv));
        if (!dst->bin_vol_inv) {
            return OSH_ENOMEM;
        }
        for (i = 0; i < dst->nbins; ++i) {
            if (i < src->nzone_indices && src->zone_volumes && src->zone_volumes[i] > 0.0) {
                dst->bin_vol_inv[i] = 1.0 / src->zone_volumes[i];
            } else {
                dst->bin_vol_inv[i] = 1.0;
            }
        }
    }
    if (src->vox_rtdose_path) {
        dst->rtdose_template_path = strdup(src->vox_rtdose_path);
        if (!dst->rtdose_template_path) {
            return OSH_ENOMEM;
        }
    }
    return OSH_OK;
}

static long find_geometry_index(struct osh_scoring_workspace const *ws, char const *name) {
    size_t i;

    for (i = 0; i < ws->ngeometries; ++i) {
        if (ws->geometries[i].name && strcmp(ws->geometries[i].name, name) == 0) {
            return (long) i;
        }
    }
    return -1;
}

static long find_filter_index(struct osh_scoring_runtime const *rt, char const *name) {
    size_t i;

    for (i = 0; i < rt->nfilters; ++i) {
        if (rt->filters[i].name && strcmp(rt->filters[i].name, name) == 0) {
            return (long) i;
        }
    }
    return -1;
}

static long find_settings_index(struct osh_scoring_runtime const *rt, char const *name) {
    size_t i;

    for (i = 0; i < rt->nsettings; ++i) {
        if (rt->settings[i].name && strcmp(rt->settings[i].name, name) == 0) {
            return (long) i;
        }
    }
    return -1;
}

/* Populate a page_override struct from one settings entry. */
static void apply_settings_to_override(struct osh_scoring_page_override *ovr,
                                       struct osh_scoring_settings_runtime const *s) {
    if (s->has_medium) {
        ovr->medium = s->medium;
        ovr->has_medium = 1;
    }
    if (s->has_density_g_cm3) {
        ovr->density_g_cm3 = s->density_g_cm3;
        ovr->has_density_g_cm3 = 1;
    }
}

void osh_scoring_runtime_finalize_ssets(struct osh_scoring_runtime *rt) {
    size_t ip;
    size_t k;
    struct osh_scoring_settings_runtime const *s;

    if (!rt) {
        return;
    }
    for (ip = 0; ip < rt->npages; ++ip) {
        struct osh_scoring_page_runtime *page = &rt->pages[ip];

        /* Quantity-level Settings override (e.g. "Quantity Dose in_Water"). */
        if (page->nsettings > 0u) {
            memset(&page->sset, 0, sizeof(page->sset));
            page->has_sset = 1;
            for (k = 0; k < page->nsettings; ++k) {
                s = &rt->settings[page->settings[k].settings_idx];
                apply_settings_to_override(&page->sset, s);
            }
        }

        /* Diff1 axis Settings override (e.g. "Diff1Type DEDX in_Si"). */
        if (page->has_diff_sset) {
            memset(&page->diff_sset, 0, sizeof(page->diff_sset));
            s = &rt->settings[page->diff_sset_idx];
            apply_settings_to_override(&page->diff_sset, s);
        }

        /* Diff2 axis Settings override (e.g. "Diff2Type LET in_Water"). */
        if (page->has_diff2_sset) {
            memset(&page->diff2_sset, 0, sizeof(page->diff2_sset));
            s = &rt->settings[page->diff2_sset_idx];
            apply_settings_to_override(&page->diff2_sset, s);
        }
    }
}

void osh_scoring_runtime_free(struct osh_scoring_runtime *rt) {
    size_t i;

    if (!rt) {
        return;
    }

    if (rt->filters) {
        for (i = 0; i < rt->nfilters; ++i) {
            free(rt->filters[i].name);
            free(rt->filters[i].rules);
        }
    }
    free(rt->filters);

    if (rt->settings) {
        for (i = 0; i < rt->nsettings; ++i) {
            free(rt->settings[i].name);
        }
    }
    free(rt->settings);

    if (rt->geometries) {
        for (i = 0; i < rt->ngeometries; ++i) {
            free(rt->geometries[i].kind);
            free(rt->geometries[i].name);
            free(rt->geometries[i].axes);
            free(rt->geometries[i].zone_indices);
            free(rt->geometries[i].groups);
            free(rt->geometries[i].rtdose_template_path);
            free(rt->geometries[i].bin_vol_inv);
        }
    }
    free(rt->geometries);

    if (rt->pages) {
        for (i = 0; i < rt->npages; ++i) {
            free(rt->pages[i].quantity);
            free(rt->pages[i].flat_rules);
            free(rt->pages[i].settings);
            osh_scoring_accumulator_free(&rt->pages[i].acc);
        }
    }
    free(rt->pages);

    if (rt->outputs) {
        for (i = 0; i < rt->noutputs; ++i) {
            free(rt->outputs[i].filename);
            free(rt->outputs[i].fileformat);
            free(rt->outputs[i].page_indices);
        }
    }
    free(rt->outputs);
    free(rt->master_scratch.crossing_buf);
    /* master_acc only shallow-aliases the per-page accumulators (freed above);
     * release the view array itself, not the arrays it points into. */
    free(rt->master_acc);

    memset(rt, 0, sizeof(*rt));
}

/* ---- Private deposit-target cloning (issue #230) ------------------------- */

void osh_scoring_runtime_free_accumulator_set(struct osh_scoring_accumulator *set, size_t npages) {
    size_t i;

    if (!set) {
        return;
    }
    /* Free only the per-page arrays; the set array itself is caller-owned. */
    for (i = 0; i < npages; ++i) {
        osh_scoring_accumulator_free(&set[i]);
    }
}

enum osh_status osh_scoring_runtime_alloc_accumulator_set(struct osh_scoring_runtime const *rt,
                                                          struct osh_scoring_accumulator *set) {
    size_t i;
    enum osh_status rc;

    if (!rt || (rt->npages > 0u && !set)) {
        return OSH_EINVAL;
    }
    for (i = 0; i < rt->npages; ++i) {
        /* Match the master page exactly: same len, same data2 presence, same
         * variance (Welford M2) presence — so a later merge always agrees on
         * optional-array presence and the private batch accumulates M2 alongside
         * the master. */
        rc = osh_scoring_accumulator_alloc_variance(
            &set[i], rt->pages[i].acc.len, rt->pages[i].acc.data2 != NULL, rt->pages[i].acc.data_var != NULL);
        if (rc != OSH_OK) {
            osh_scoring_runtime_free_accumulator_set(set, i); /* release the pages built so far */
            return rc;
        }
    }
    return OSH_OK;
}

void osh_scoring_runtime_free_scratch(struct osh_scoring_scratch *scratch) {
    if (!scratch) {
        return;
    }
    free(scratch->crossing_buf);
    scratch->crossing_buf = NULL;
    scratch->crossing_cap = 0u;
}

enum osh_status osh_scoring_runtime_clone_scratch(struct osh_scoring_runtime const *rt,
                                                  struct osh_scoring_scratch *scratch_out) {
    if (!rt || !scratch_out) {
        return OSH_EINVAL;
    }
    /* @p scratch_out is a pure output: reject a struct that still owns a buffer so
     * a mistaken reuse surfaces as an error instead of silently leaking it.  The
     * caller must hand us a zero-initialised (or already-freed) scratch. */
    if (scratch_out->crossing_buf != NULL) {
        return OSH_EINVAL;
    }
    scratch_out->crossing_buf = NULL;
    scratch_out->crossing_cap = 0u;
    /* Mirror the master scratch: a runtime with no crossing geometry keeps a NULL
     * buffer (cap 0), so the clone does too. */
    if (rt->master_scratch.crossing_cap == 0u) {
        return OSH_OK;
    }
    scratch_out->crossing_buf =
        (struct osh_voxel_crossing *) malloc(rt->master_scratch.crossing_cap * sizeof(*scratch_out->crossing_buf));
    if (!scratch_out->crossing_buf) {
        return OSH_ENOMEM;
    }
    scratch_out->crossing_cap = rt->master_scratch.crossing_cap;
    return OSH_OK;
}

enum osh_status osh_scoring_compile(struct osh_scoring_workspace const *ws,
                                    struct osh_diag_sink const *diag,
                                    struct osh_scoring_runtime *rt) {
    enum osh_status rc;
    size_t i;
    size_t j;
    size_t k;
    size_t total_pages;
    size_t ordinal; /* insertion-order counter; preserved inside each sort group */
    size_t page_idx;
    size_t current_group;
    long gidx; /* geometry index; -1 = not found */
    long fidx; /* filter index;   -1 = not found */
    long sidx; /* settings index; -1 = not found */
    size_t nf; /* number of filter references on current page */
    size_t ns; /* number of settings references on current page */
    enum osh_scoring_score_kind score_kind;
    long *output_geom_idx = NULL;                    /* per-output resolved geometry index */
    size_t *geom_page_counts = NULL;                 /* how many pages each geometry owns */
    struct prepared_page_ref *prepared_pages = NULL; /* scratch sort buffer */
    struct osh_scoring_output_runtime *out;
    struct osh_scoring_page_def const *src_page;
    struct osh_scoring_page_runtime *dst_page;

    if (!ws || !rt) {
        return OSH_EINVAL;
    }

    memset(rt, 0, sizeof(*rt));
    rt->nfilters = ws->nfilters;
    rt->nsettings = ws->nsettings;
    rt->ngeometries = ws->ngeometries;
    rt->noutputs = ws->noutputs;

    /* --- Phase 1: copy filters, settings, and geometries into runtime arrays. --- */

    if (rt->nfilters > 0u) {
        rt->filters = (struct osh_scoring_filter_runtime *) calloc(rt->nfilters, sizeof(*rt->filters));
        if (!rt->filters) {
            rc = OSH_ENOMEM;
            goto fail;
        }
        for (i = 0; i < rt->nfilters; ++i) {
            rc = copy_filter_runtime(&rt->filters[i], diag, &ws->filters[i]);
            if (rc != OSH_OK) {
                goto fail;
            }
        }
    }

    if (rt->nsettings > 0u) {
        rt->settings = (struct osh_scoring_settings_runtime *) calloc(rt->nsettings, sizeof(*rt->settings));
        if (!rt->settings) {
            rc = OSH_ENOMEM;
            goto fail;
        }
        for (i = 0; i < rt->nsettings; ++i) {
            rc = copy_settings_runtime(&rt->settings[i], &ws->settings[i]);
            if (rc != OSH_OK) {
                goto fail;
            }
        }
    }

    if (rt->ngeometries > 0u) {
        rt->geometries = (struct osh_scoring_geometry_runtime *) calloc(rt->ngeometries, sizeof(*rt->geometries));
        if (!rt->geometries) {
            rc = OSH_ENOMEM;
            goto fail;
        }
        for (i = 0; i < rt->ngeometries; ++i) {
            rc = copy_geometry_runtime(&rt->geometries[i], &ws->geometries[i]);
            if (rc != OSH_OK) {
                if (rc == OSH_EINVAL) {
                    OSH_DIAG_ERRORF(diag, "Scoring geometry '%s' has no valid runtime binning", ws->geometries[i].name);
                }
                goto fail;
            }
        }
    }

    /* --- Phase 2: validate outputs — resolve geometry names, reject unsupported
     *              score kinds, and count how many pages each geometry owns. --- */

    output_geom_idx = (long *) calloc(ws->noutputs ? ws->noutputs : 1u, sizeof(*output_geom_idx));
    geom_page_counts = (size_t *) calloc(rt->ngeometries ? rt->ngeometries : 1u, sizeof(*geom_page_counts));
    if (!output_geom_idx || !geom_page_counts) {
        rc = OSH_ENOMEM;
        goto fail;
    }

    total_pages = 0u;
    for (i = 0; i < ws->noutputs; ++i) {
        gidx = find_geometry_index(ws, ws->outputs[i].geometry_name);
        if (gidx < 0) {
            OSH_DIAG_ERRORF(diag,
                            "Scoring output '%s' references unknown geometry '%s'",
                            ws->outputs[i].filename ? ws->outputs[i].filename : "(unnamed)",
                            ws->outputs[i].geometry_name ? ws->outputs[i].geometry_name : "(null)");
            rc = OSH_EINVAL;
            goto fail;
        }
        output_geom_idx[i] = gidx;
        if (!runtime_supports_geometry(&rt->geometries[gidx])) {
            OSH_DIAG_ERRORF(diag,
                            "Scoring output '%s' uses unsupported runtime geometry '%s' (kind=%s)",
                            ws->outputs[i].filename ? ws->outputs[i].filename : "(unnamed)",
                            ws->outputs[i].geometry_name ? ws->outputs[i].geometry_name : "(null)",
                            ws->geometries[gidx].kind ? ws->geometries[gidx].kind : "(unknown)");
            rc = OSH_ENOTSUP;
            goto fail;
        }
        for (j = 0; j < ws->outputs[i].npages; ++j) {
            score_kind = quantity_to_score_kind(ws->outputs[i].pages[j].quantity);
            if (!runtime_supports_score_kind(score_kind)) {
                OSH_DIAG_ERRORF(diag,
                                "Scoring output '%s' uses unsupported quantity '%s' for runtime scoring",
                                ws->outputs[i].filename ? ws->outputs[i].filename : "(unnamed)",
                                ws->outputs[i].pages[j].quantity ? ws->outputs[i].pages[j].quantity : "(null)");
                rc = OSH_ENOTSUP;
                goto fail;
            }
            if (rt->geometries[gidx].geo_kind == OSH_SCORING_GEO_ZONE && score_kind != OSH_SCORING_SCORE_ENERGY
                && score_kind != OSH_SCORING_SCORE_FLUENCE && score_kind != OSH_SCORING_SCORE_DOSE
                && score_kind != OSH_SCORING_SCORE_DOSEGY && score_kind != OSH_SCORING_SCORE_DIRTYDOSE
                && score_kind != OSH_SCORING_SCORE_DIRTYDOSEGY) {
                OSH_DIAG_ERRORF(diag,
                                "Scoring output '%s' uses quantity '%s' on Zone geometry '%s'; only Energy, Fluence, "
                                "Dose, DoseGy, DirtyDose, and DirtyDoseGy are supported for Zone scoring",
                                ws->outputs[i].filename ? ws->outputs[i].filename : "(unnamed)",
                                ws->outputs[i].pages[j].quantity ? ws->outputs[i].pages[j].quantity : "(null)",
                                ws->outputs[i].geometry_name ? ws->outputs[i].geometry_name : "(null)");
                rc = OSH_ENOTSUP;
                goto fail;
            }
        }
        geom_page_counts[gidx] += ws->outputs[i].npages;
        total_pages += ws->outputs[i].npages;
    }

    /* --- Phase 3: allocate flat page and output arrays; assign each geometry a
     *              contiguous slice of the page array (first_page + npages). --- */

    rt->npages = total_pages;
    if (rt->npages > 0u) {
        rt->pages = (struct osh_scoring_page_runtime *) calloc(rt->npages, sizeof(*rt->pages));
        prepared_pages = (struct prepared_page_ref *) calloc(rt->npages, sizeof(*prepared_pages));
        if (!rt->pages || !prepared_pages) {
            rc = OSH_ENOMEM;
            goto fail;
        }
    }
    if (rt->noutputs > 0u) {
        rt->outputs = (struct osh_scoring_output_runtime *) calloc(rt->noutputs, sizeof(*rt->outputs));
        if (!rt->outputs) {
            rc = OSH_ENOMEM;
            goto fail;
        }
    }

    total_pages = 0u;
    for (i = 0; i < rt->ngeometries; ++i) {
        rt->geometries[i].first_page = total_pages;
        rt->geometries[i].npages = geom_page_counts[i];
        total_pages += geom_page_counts[i];
    }

    /* --- Phase 4: build prepared_pages scratch array used for sorting. --- */

    ordinal = 0u;
    for (i = 0; i < ws->noutputs; ++i) {
        char const *primary_format;

        out = &rt->outputs[i];

        out->filename = strdup(ws->outputs[i].filename);
        if (!out->filename) {
            rc = OSH_ENOMEM;
            goto fail;
        }
        /* The block's first requested format (lowercase, as stored by the
         * parser); an empty list defaults to BDO.  Additional formats become
         * extra runtime outputs in Phase 7 (issue #308). */
        primary_format = "bdo";
        if (ws->outputs[i].nfileformats > 0u) {
            primary_format = ws->outputs[i].fileformats[0];
        }
        out->fileformat = strdup(primary_format);
        if (!out->fileformat) {
            rc = OSH_ENOMEM;
            goto fail;
        }
        out->geometry_idx = (size_t) output_geom_idx[i];
        out->npages = ws->outputs[i].npages;
        if (out->npages > 0u) {
            out->page_indices = (size_t *) calloc(out->npages, sizeof(*out->page_indices));
            if (!out->page_indices) {
                rc = OSH_ENOMEM;
                goto fail;
            }
        }

        for (j = 0; j < ws->outputs[i].npages; ++j) {
            src_page = &ws->outputs[i].pages[j];
            prepared_pages[ordinal].output_idx = i;
            prepared_pages[ordinal].src_page_idx = j;
            prepared_pages[ordinal].ordinal = ordinal;
            prepared_pages[ordinal].geometry_idx = out->geometry_idx;
            prepared_pages[ordinal].score_kind = quantity_to_score_kind(src_page->quantity);
            ordinal++;
        }
    }

    /* Sort so pages with the same geometry+score_kind are contiguous, enabling
     * the hot path to walk a single group in one tight inner loop. */
    if (rt->npages > 1u) {
        qsort(prepared_pages, rt->npages, sizeof(*prepared_pages), compare_prepared_pages);
    }

    /* --- Phase 5: populate each runtime page from the sorted scratch records. --- */

    for (page_idx = 0; page_idx < rt->npages; ++page_idx) {
        out = &rt->outputs[prepared_pages[page_idx].output_idx];
        src_page = &ws->outputs[prepared_pages[page_idx].output_idx].pages[prepared_pages[page_idx].src_page_idx];
        dst_page = &rt->pages[page_idx];

        if (!out->page_indices || prepared_pages[page_idx].src_page_idx >= out->npages) {
            rc = OSH_ESTATE;
            goto fail;
        }
        out->page_indices[prepared_pages[page_idx].src_page_idx] = page_idx;
        dst_page->quantity = strdup(src_page->quantity);
        if (!dst_page->quantity) {
            rc = OSH_ENOMEM;
            goto fail;
        }
        dst_page->output_idx = prepared_pages[page_idx].output_idx;
        dst_page->geometry_idx = prepared_pages[page_idx].geometry_idx;
        dst_page->score_kind = prepared_pages[page_idx].score_kind;
        dst_page->has_data2 = score_kind_uses_data2(dst_page->score_kind);
        dst_page->divide = dst_page->has_data2;
        dst_page->postproc = score_kind_postproc(dst_page->score_kind, dst_page->divide);
        /* Monte-Carlo standard-error tracking (issue #209) is enabled per
         * estimator by attaching a Settings block carrying "Variance On" to the
         * Quantity line (e.g. "Quantity Dose withErr").  Resolve it here — before
         * the accumulator is allocated below — because it decides the accumulator's
         * memory layout (the batch-means M2 companion arrays).  The general Settings
         * override (medium/density) is applied later in finalize_ssets(), but that
         * runs after this allocation, so variance is resolved directly against the
         * page's Settings references now.  0 keeps the accumulator exactly as before
         * (no var arrays); the per-batch fold and finalize happen in transport and
         * postprocess. */
        dst_page->variance = 0;
        for (k = 0; k < src_page->nfilter_names; ++k) {
            sidx = find_settings_index(rt, src_page->filter_names[k]);
            if (sidx >= 0 && rt->settings[(size_t) sidx].has_variance && rt->settings[(size_t) sidx].variance) {
                dst_page->variance = 1;
            }
        }

        /* Differential axis — expand data[] to geo_nbins × diff_nbins. */
        dst_page->diff_stride = rt->geometries[dst_page->geometry_idx].nbins;
        if (src_page->diff_nbins > 0u) {
            if (dst_page->has_data2) {
                /* Two-pass averaged quantities (DLET, TLET, DQEFF, TQEFF, DAVGE,
                 * TAVGE, DBETA, TBETA) cannot be used as the primary quantity of a
                 * differential scorer because the two-pass postprocessor operates
                 * on the full flat array. */
                OSH_DIAG_ERRORF(diag,
                                "Scoring page '%s': differential axis not supported for averaged quantities",
                                src_page->quantity ? src_page->quantity : "(unnamed)");
                rc = OSH_ENOTSUP;
                goto fail;
            }
            dst_page->diff_nbins = src_page->diff_nbins;
            dst_page->diff_lo = src_page->diff_lo;
            dst_page->diff_hi = src_page->diff_hi;
            dst_page->diff_log = src_page->diff_log;
            dst_page->diff_kind = diff_kind_from_str(src_page->diff_kind_str);
            if (dst_page->diff_kind == OSH_SCORING_DIFF_NONE) {
                OSH_DIAG_ERRORF(diag,
                                "Scoring page '%s': unknown Diff1Type '%s'",
                                src_page->quantity ? src_page->quantity : "(unnamed)",
                                src_page->diff_kind_str ? src_page->diff_kind_str : "(null)");
                rc = OSH_ENOTSUP;
                goto fail;
            }
            dst_page->len = dst_page->diff_stride * dst_page->diff_nbins;

            /* Optional per-axis Settings override — e.g. "Diff1Type DEDX in_Si".
             * Resolve the name to a settings index now; the medium index is
             * filled by osh_scoring_runtime_finalize_ssets() after material names
             * are resolved by the simulation layer. */
            if (src_page->diff_kind_sset_name) {
                long sidx = find_settings_index(rt, src_page->diff_kind_sset_name);
                if (sidx < 0L) {
                    OSH_DIAG_ERRORF(diag,
                                    "Scoring page '%s': Diff1Type references unknown Settings '%s'",
                                    src_page->quantity ? src_page->quantity : "(unnamed)",
                                    src_page->diff_kind_sset_name);
                    rc = OSH_ENOTSUP;
                    goto fail;
                }
                dst_page->diff_sset_idx = (size_t) sidx;
                dst_page->has_diff_sset = 1;
            }

            /* Double-differential axis (requires diff1 to be set). */
            if (src_page->diff2_nbins > 0u) {
                dst_page->diff2_nbins = src_page->diff2_nbins;
                dst_page->diff2_lo = src_page->diff2_lo;
                dst_page->diff2_hi = src_page->diff2_hi;
                dst_page->diff2_log = src_page->diff2_log;
                dst_page->diff2_kind = diff_kind_from_str(src_page->diff2_kind_str);
                if (dst_page->diff2_kind == OSH_SCORING_DIFF_NONE) {
                    OSH_DIAG_ERRORF(diag,
                                    "Scoring page '%s': unknown Diff2Type '%s'",
                                    src_page->quantity ? src_page->quantity : "(unnamed)",
                                    src_page->diff2_kind_str ? src_page->diff2_kind_str : "(null)");
                    rc = OSH_ENOTSUP;
                    goto fail;
                }
                /* diff2_stride = diff_nbins * diff_stride = geo_nbins * diff_nbins. */
                dst_page->diff2_stride = dst_page->diff_nbins * dst_page->diff_stride;
                dst_page->len = dst_page->diff2_stride * dst_page->diff2_nbins;

                /* Optional per-axis Settings override for the diff2 axis. */
                if (src_page->diff2_kind_sset_name) {
                    long sidx = find_settings_index(rt, src_page->diff2_kind_sset_name);
                    if (sidx < 0L) {
                        OSH_DIAG_ERRORF(diag,
                                        "Scoring page '%s': Diff2Type references unknown Settings '%s'",
                                        src_page->quantity ? src_page->quantity : "(unnamed)",
                                        src_page->diff2_kind_sset_name);
                        rc = OSH_ENOTSUP;
                        goto fail;
                    }
                    dst_page->diff2_sset_idx = (size_t) sidx;
                    dst_page->has_diff2_sset = 1;
                }
            }
        } else {
            dst_page->diff_nbins = 0u;
            dst_page->len = dst_page->diff_stride;
        }

        rc = osh_scoring_accumulator_alloc_variance(
            &dst_page->acc, dst_page->len, dst_page->has_data2, dst_page->variance);
        if (rc != OSH_OK) {
            goto fail;
        }

        if (src_page->nfilter_names > 0u) {
            size_t total_rules;
            size_t r;
            /* Pass 1: classify names; count flat rules by summing each filter's rule count. */
            nf = 0u;
            ns = 0u;
            total_rules = 0u;
            for (k = 0; k < src_page->nfilter_names; ++k) {
                fidx = find_filter_index(rt, src_page->filter_names[k]);
                if (fidx >= 0) {
                    nf++;
                    total_rules += rt->filters[(size_t) fidx].nrules;
                } else {
                    sidx = find_settings_index(rt, src_page->filter_names[k]);
                    if (sidx >= 0) {
                        ns++;
                    } else {
                        OSH_DIAG_ERRORF(diag,
                                        "Scoring page '%s' references unknown filter or settings '%s'",
                                        src_page->quantity ? src_page->quantity : "(unnamed)",
                                        src_page->filter_names[k]);
                        rc = OSH_EINVAL;
                        goto fail;
                    }
                }
            }
            /* Pass 2: allocate flat_rules and settings, then fill them. */
            if (total_rules > 0u) {
                dst_page->flat_rules =
                    (struct osh_scoring_filter_runtime_rule *) calloc(total_rules, sizeof(*dst_page->flat_rules));
                if (!dst_page->flat_rules) {
                    rc = OSH_ENOMEM;
                    goto fail;
                }
            }
            if (ns > 0u) {
                dst_page->settings = (struct osh_scoring_page_settings_ref *) calloc(ns, sizeof(*dst_page->settings));
                if (!dst_page->settings) {
                    rc = OSH_ENOMEM;
                    goto fail;
                }
            }
            dst_page->nflat_rules = 0u;
            dst_page->nsettings = ns;
            ns = 0u;
            for (k = 0; k < src_page->nfilter_names; ++k) {
                fidx = find_filter_index(rt, src_page->filter_names[k]);
                if (fidx >= 0) {
                    struct osh_scoring_filter_runtime const *f = &rt->filters[(size_t) fidx];
                    for (r = 0; r < f->nrules; ++r) {
                        dst_page->flat_rules[dst_page->nflat_rules++] = f->rules[r];
                    }
                } else {
                    sidx = find_settings_index(rt, src_page->filter_names[k]);
                    dst_page->settings[ns++].settings_idx = (size_t) sidx;
                }
            }

            /* sset is built after all pages are processed; see the call below. */
        }
    }

    /* --- Phase 6: build score groups — one group per contiguous run of pages
     *              that share the same score_kind within a geometry.  The hot
     *              path iterates groups so it can hoist per-kind setup once. --- */

    for (i = 0; i < rt->ngeometries; ++i) {
        if (rt->geometries[i].npages == 0u) {
            continue;
        }

        rt->geometries[i].ngroups = 1u;
        for (j = 1; j < rt->geometries[i].npages; ++j) {
            page_idx = rt->geometries[i].first_page + j;
            if (rt->pages[page_idx - 1u].score_kind != rt->pages[page_idx].score_kind) {
                rt->geometries[i].ngroups++;
            }
        }

        rt->geometries[i].groups = (struct osh_scoring_geometry_score_group *) calloc(
            rt->geometries[i].ngroups, sizeof(*rt->geometries[i].groups));
        if (!rt->geometries[i].groups) {
            rc = OSH_ENOMEM;
            goto fail;
        }

        current_group = 0u;
        rt->geometries[i].groups[0].first_page = rt->geometries[i].first_page;
        rt->geometries[i].groups[0].npages = 1u;
        rt->geometries[i].groups[0].score_kind = rt->pages[rt->geometries[i].first_page].score_kind;

        for (j = 1; j < rt->geometries[i].npages; ++j) {
            page_idx = rt->geometries[i].first_page + j;
            if (rt->pages[page_idx].score_kind == rt->geometries[i].groups[current_group].score_kind) {
                rt->geometries[i].groups[current_group].npages++;
            } else {
                current_group++;
                rt->geometries[i].groups[current_group].first_page = page_idx;
                rt->geometries[i].groups[current_group].npages = 1u;
                rt->geometries[i].groups[current_group].score_kind = rt->pages[page_idx].score_kind;
            }
        }
    }

    /* --- Phase 7: expand multi-format outputs (issue #308). ---
     * A block that requests several formats writes the *same* scored pages once
     * per format; osh_scoring_expand_multiformat_outputs() adds one extra runtime
     * output per additional format, all sharing the block's page indices (a cheap
     * size_t copy — nothing in rt->pages[] is duplicated), so scoring memory is
     * independent of the format count.  See osh_scoring_multiformat.c. */
    rc = osh_scoring_expand_multiformat_outputs(ws, diag, rt);
    if (rc != OSH_OK) {
        goto fail;
    }

    free(output_geom_idx);
    free(geom_page_counts);
    free(prepared_pages);

    /* Pre-allocate the serial driver's per-step voxel-crossing scratch buffer.
     * Size it to the largest per-geometry cap, computed with the same formula
     * each traversal path uses in osh_scoring_step.c: mesh needs n[0]+n[1]+n[2]
     * (sum of axis nbins), while CYL needs 2*nr+nz (== 2*n[0]+n[2]).  The two
     * cases are handled separately below; keep both in sync with score_step so
     * the buffer never overflows.  One buffer is reused for every geometry on every
     * physics step, eliminating the calloc/free pair that was spending ~60% of
     * total CPU time on memset.  It lives in rt->master_scratch and is handed to
     * osh_scoring_score_step() by the serial driver via
     * osh_scoring_runtime_master_scratch(); a parallel worker owns its own scratch
     * instead, so the buffer is never shared mutable state on the deposit path. */
    {
        size_t max_cap = 0;
        for (i = 0; i < rt->ngeometries; ++i) {
            size_t cap_i = 0;
            size_t a;
            if (rt->geometries[i].geo_kind == OSH_SCORING_GEO_CYL) {
                /* CYL: max crossings = 2*nr + nz. Look up R/Z by label so that
                 * axis declaration order in detect.dat does not matter. */
                size_t nr = 0, nz = 0;
                for (a = 0; a < rt->geometries[i].naxes; ++a) {
                    if (strcmp(rt->geometries[i].axes[a].label, "R") == 0)
                        nr = (size_t) rt->geometries[i].axes[a].nbins;
                    else if (strcmp(rt->geometries[i].axes[a].label, "Z") == 0)
                        nz = (size_t) rt->geometries[i].axes[a].nbins;
                }
                cap_i = 2u * nr + nz;
            } else {
                for (a = 0; a < rt->geometries[i].naxes; ++a) {
                    cap_i += (size_t) rt->geometries[i].axes[a].nbins;
                }
            }
            if (cap_i > max_cap) {
                max_cap = cap_i;
            }
        }
        if (max_cap > 0) {
            rt->master_scratch.crossing_buf =
                (struct osh_voxel_crossing *) malloc(max_cap * sizeof(*rt->master_scratch.crossing_buf));
            if (!rt->master_scratch.crossing_buf) {
                osh_scoring_runtime_free(rt);
                return OSH_ENOMEM;
            }
            rt->master_scratch.crossing_cap = max_cap;
        }
    }

    /* Precompute the geometry-agnostic per-spatial-bin 1/volume that the
     * volume-normalised estimators (DOSE, FLUENCE, ...) apply in postprocess: the
     * scorer deposits the extensive quantity and postprocess divides by volume once
     * per bin.  bin_vol_inv is the single volume source.  Zone geometries already
     * built it in copy_geometry_runtime (from the user's per-zone Volume); here we
     * build Mesh (uniform) and Cyl (per-R shell). */
    for (i = 0; i < rt->ngeometries; ++i) {
        struct osh_scoring_geometry_runtime *g = &rt->geometries[i];
        size_t nb = g->nbins;
        size_t b;

        if (nb == 0u || g->geo_kind == OSH_SCORING_GEO_ZONE) {
            continue;
        }
        g->bin_vol_inv = (double *) malloc(nb * sizeof(*g->bin_vol_inv));
        if (!g->bin_vol_inv) {
            osh_scoring_runtime_free(rt);
            return OSH_ENOMEM;
        }
        if (g->geo_kind == OSH_SCORING_GEO_MESH) {
            /* Uniform voxel volume: V = (product of axis extents) / nbins, so
             * 1/V = nbins / extent, identical for every bin. */
            double extent = 1.0;
            double vinv;
            size_t a;
            for (a = 0u; a < g->naxes; ++a) {
                extent *= (g->axes[a].hi - g->axes[a].lo);
            }
            if (extent > 0.0) {
                vinv = (double) nb / extent;
            } else {
                vinv = 1.0;
            }
            for (b = 0u; b < nb; ++b) {
                g->bin_vol_inv[b] = vinv;
            }
        } else if (g->geo_kind == OSH_SCORING_GEO_CYL) {
            /* Cylindrical shell volume depends only on r_bin (flat bin =
             * z_bin*nr + r_bin): V = pi * (r1^2 - r0^2) * dz.  Resolve R/Z by label. */
            size_t a;
            size_t nr = 0u;
            size_t nz_bins = 0u;
            size_t r_bin;
            double r_lo = 0.0;
            double r_hi = 0.0;
            double z_lo = 0.0;
            double z_hi = 0.0;
            double dr;
            double dz;
            double r0;
            double r1;
            for (a = 0u; a < g->naxes; ++a) {
                if (strcmp(g->axes[a].label, "R") == 0) {
                    nr = (size_t) g->axes[a].nbins;
                    r_lo = g->axes[a].lo;
                    r_hi = g->axes[a].hi;
                } else if (strcmp(g->axes[a].label, "Z") == 0) {
                    nz_bins = (size_t) g->axes[a].nbins;
                    z_lo = g->axes[a].lo;
                    z_hi = g->axes[a].hi;
                }
            }
            if (nr == 0u || nz_bins == 0u) {
                for (b = 0u; b < nb; ++b) {
                    g->bin_vol_inv[b] = 1.0; /* malformed; compile already rejected */
                }
            } else {
                dr = (r_hi - r_lo) / (double) nr;
                dz = (z_hi - z_lo) / (double) nz_bins;
                for (b = 0u; b < nb; ++b) {
                    r_bin = b % nr;
                    r0 = r_lo + (double) r_bin * dr;
                    r1 = r0 + dr;
                    g->bin_vol_inv[b] = 1.0 / (OSH_M_PI * (r1 * r1 - r0 * r0) * dz);
                }
            }
        } else {
            for (b = 0u; b < nb; ++b) {
                g->bin_vol_inv[b] = 1.0;
            }
        }
    }

    /* Build the master accumulator view: one entry per page, shallow-aliasing
     * each page's accumulator storage.  The serial driver hands this to
     * osh_scoring_score_step() so the deposit path has the same shape as a
     * parallel worker's private set, with no per-step allocation. */
    if (rt->npages > 0u) {
        rt->master_acc = (struct osh_scoring_accumulator *) calloc(rt->npages, sizeof(*rt->master_acc));
        if (!rt->master_acc) {
            osh_scoring_runtime_free(rt);
            return OSH_ENOMEM;
        }
        for (i = 0; i < rt->npages; ++i) {
            rt->master_acc[i] = rt->pages[i].acc;
        }
    }

    osh_scoring_runtime_finalize_ssets(rt);
    return OSH_OK;

fail:
    free(output_geom_idx);
    free(geom_page_counts);
    free(prepared_pages);
    osh_scoring_runtime_free(rt);
    return rc;
}
