#include "scoring/runtime/osh_scoring_compile.h"

#include <stdlib.h>
#include <string.h>

#include "common/osh_diag.h"
#include "common/raytrace/osh_raytrace.h"

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

    if (geo->zone_stop >= geo->zone_start && geo->zone_start > 0) {
        return (size_t) (geo->zone_stop - geo->zone_start + 1);
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
        return 1;
    default:
        return 0;
    }
}

/* divide=1 means postproc divides data by data2 (LET average) → AVER mode.
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

/* Only Cartesian mesh geometry is wired to the hot-path scorer for now. */
static int runtime_supports_geometry(struct osh_scoring_geometry_runtime const *geo) {
    if (!geo) {
        return 0;
    }
    return geo->geo_kind == OSH_SCORING_GEO_MESH;
}

static int runtime_supports_score_kind(enum osh_scoring_score_kind score_kind) {
    switch (score_kind) {
    case OSH_SCORING_SCORE_ENERGY:
    case OSH_SCORING_SCORE_FLUENCE:
    case OSH_SCORING_SCORE_DOSE:
    case OSH_SCORING_SCORE_DLET:
    case OSH_SCORING_SCORE_TLET:
    case OSH_SCORING_SCORE_DQEFF:
    case OSH_SCORING_SCORE_TQEFF:
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
    dst->has_rescale = src->has_rescale;
    dst->has_offset = src->has_offset;
    dst->has_site_diameter_um = src->has_site_diameter_um;
    dst->has_density_g_cm3 = src->has_density_g_cm3;
    dst->has_npart = src->has_npart;
    dst->has_medium = src->has_medium;
    dst->has_nkmedium = src->has_nkmedium;
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
    dst->zone_start = src->zone_start;
    dst->zone_stop = src->zone_stop;
    dst->has_rotation = src->has_rotation;
    dst->geo_kind = geometry_kind_to_enum(src->kind);
    dst->nbins = geometry_nbins(src);
    if (dst->nbins == 0u) {
        return OSH_EINVAL;
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

void osh_scoring_runtime_finalize_ssets(struct osh_scoring_runtime *rt) {
    size_t ip;
    size_t k;

    if (!rt) {
        return;
    }
    for (ip = 0; ip < rt->npages; ++ip) {
        struct osh_scoring_page_runtime *page = &rt->pages[ip];
        if (page->nsettings == 0u) {
            continue;
        }
        memset(&page->sset, 0, sizeof(page->sset));
        page->has_sset = 1;
        for (k = 0; k < page->nsettings; ++k) {
            struct osh_scoring_settings_runtime const *s = &rt->settings[page->settings[k].settings_idx];
            if (s->has_medium) {
                page->sset.medium = s->medium;
                page->sset.has_medium = 1;
            }
            if (s->has_density_g_cm3) {
                page->sset.density_g_cm3 = s->density_g_cm3;
                page->sset.has_density_g_cm3 = 1;
            }
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
            free(rt->geometries[i].groups);
            free(rt->geometries[i].rtdose_template_path);
        }
    }
    free(rt->geometries);

    if (rt->pages) {
        for (i = 0; i < rt->npages; ++i) {
            free(rt->pages[i].quantity);
            free(rt->pages[i].flat_rules);
            free(rt->pages[i].settings);
            free(rt->pages[i].data);
            free(rt->pages[i].data_var);
            free(rt->pages[i].data2);
            free(rt->pages[i].data2_var);
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
    free(rt->crossing_buf);

    memset(rt, 0, sizeof(*rt));
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
        out = &rt->outputs[i];

        out->filename = strdup(ws->outputs[i].filename);
        if (!out->filename) {
            rc = OSH_ENOMEM;
            goto fail;
        }
        out->fileformat = strdup(ws->outputs[i].fileformat ? ws->outputs[i].fileformat : "BDO");
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
        dst_page->len = rt->geometries[dst_page->geometry_idx].nbins;
        dst_page->score_kind = prepared_pages[page_idx].score_kind;
        dst_page->has_data2 = score_kind_uses_data2(dst_page->score_kind);
        dst_page->divide = dst_page->has_data2;
        dst_page->postproc = score_kind_postproc(dst_page->score_kind, dst_page->divide);
        dst_page->data = (double *) calloc(dst_page->len ? dst_page->len : 1u, sizeof(*dst_page->data));
        if (!dst_page->data) {
            rc = OSH_ENOMEM;
            goto fail;
        }
        if (dst_page->has_data2) {
            dst_page->data2 = (double *) calloc(dst_page->len ? dst_page->len : 1u, sizeof(*dst_page->data2));
            if (!dst_page->data2) {
                rc = OSH_ENOMEM;
                goto fail;
            }
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

    free(output_geom_idx);
    free(geom_page_counts);
    free(prepared_pages);

    /* Pre-allocate the per-step voxel-crossing scratch buffer.  The maximum
     * capacity needed is the sum of nbins across all axes of the largest
     * geometry (matching the cap = n[0]+n[1]+n[2] formula in osh_scoring_step.c).
     * One buffer is reused for every geometry on every physics step, eliminating
     * the calloc/free pair that was spending ~60% of total CPU time on memset. */
    {
        size_t max_cap = 0;
        for (i = 0; i < rt->ngeometries; ++i) {
            size_t cap_i = 0;
            size_t a;
            for (a = 0; a < rt->geometries[i].naxes; ++a) {
                cap_i += (size_t) rt->geometries[i].axes[a].nbins;
            }
            if (cap_i > max_cap) {
                max_cap = cap_i;
            }
        }
        if (max_cap > 0) {
            rt->crossing_buf = (struct osh_voxel_crossing *) malloc(max_cap * sizeof(*rt->crossing_buf));
            if (!rt->crossing_buf) {
                osh_scoring_runtime_free(rt);
                return OSH_ENOMEM;
            }
            rt->crossing_cap = max_cap;
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
