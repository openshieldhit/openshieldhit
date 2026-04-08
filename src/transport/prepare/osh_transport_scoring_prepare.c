#include "transport/prepare/osh_transport_scoring_prepare.h"

#include <stdlib.h>
#include <string.h>

#include "common/osh_logger.h"

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

static enum osh_status copy_filter_runtime(struct osh_transport_scoring_filter_runtime *dst,
                                           struct osh_scoring_filter_def const *src) {
    size_t i;

    memset(dst, 0, sizeof(*dst));
    dst->name = strdup(src->name);
    if (!dst->name) {
        return OSH_ENOMEM;
    }
    if (src->nrules == 0u) {
        return OSH_OK;
    }

    dst->rules = (struct osh_transport_scoring_filter_rule *) calloc(src->nrules, sizeof(*dst->rules));
    if (!dst->rules) {
        return OSH_ENOMEM;
    }
    dst->nrules = src->nrules;
    for (i = 0; i < src->nrules; ++i) {
        memcpy(dst->rules[i].field, src->rules[i].field, sizeof(dst->rules[i].field));
        memcpy(dst->rules[i].op, src->rules[i].op, sizeof(dst->rules[i].op));
        dst->rules[i].value = src->rules[i].value;
    }
    return OSH_OK;
}

static enum osh_status copy_settings_runtime(struct osh_transport_scoring_settings_runtime *dst,
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

static enum osh_status copy_geometry_runtime(struct osh_transport_scoring_geometry_runtime *dst,
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
        dst->axes = (struct osh_transport_scoring_axis_runtime *) calloc(src->naxes, sizeof(*dst->axes));
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
    dst->rot_theta_deg = src->rot_theta_deg;
    dst->rot_phi_deg = src->rot_phi_deg;
    dst->zone_start = src->zone_start;
    dst->zone_stop = src->zone_stop;
    dst->has_rotation = src->has_rotation;
    dst->nbins = geometry_nbins(src);
    if (dst->nbins == 0u) {
        return OSH_EINVAL;
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

static long find_filter_index(struct osh_transport_scoring_runtime const *rt, char const *name) {
    size_t i;

    for (i = 0; i < rt->nfilters; ++i) {
        if (rt->filters[i].name && strcmp(rt->filters[i].name, name) == 0) {
            return (long) i;
        }
    }
    return -1;
}

void osh_transport_scoring_runtime_free(struct osh_transport_scoring_runtime *rt) {
    size_t i;

    if (!rt) {
        return;
    }

    for (i = 0; i < rt->nfilters; ++i) {
        free(rt->filters[i].name);
        free(rt->filters[i].rules);
    }
    free(rt->filters);

    for (i = 0; i < rt->nsettings; ++i) {
        free(rt->settings[i].name);
    }
    free(rt->settings);

    for (i = 0; i < rt->ngeometries; ++i) {
        free(rt->geometries[i].kind);
        free(rt->geometries[i].name);
        free(rt->geometries[i].axes);
    }
    free(rt->geometries);

    for (i = 0; i < rt->npages; ++i) {
        free(rt->pages[i].quantity);
        free(rt->pages[i].filters);
        free(rt->pages[i].settings);
        free(rt->pages[i].data);
        free(rt->pages[i].data_var);
        free(rt->pages[i].data2);
        free(rt->pages[i].data2_var);
    }
    free(rt->pages);

    for (i = 0; i < rt->noutputs; ++i) {
        free(rt->outputs[i].filename);
        free(rt->outputs[i].fileformat);
        free(rt->outputs[i].page_indices);
    }
    free(rt->outputs);

    memset(rt, 0, sizeof(*rt));
}

enum osh_status osh_transport_scoring_prepare(struct osh_scoring_workspace const *ws,
                                              struct osh_transport_scoring_runtime *rt) {
    enum osh_status rc;
    size_t i;
    size_t j;
    size_t total_pages;
    size_t *geom_page_counts = NULL;
    size_t *geom_next_page = NULL;

    if (!ws || !rt) {
        return OSH_EINVAL;
    }

    memset(rt, 0, sizeof(*rt));
    rt->nfilters = ws->nfilters;
    rt->nsettings = ws->nsettings;
    rt->ngeometries = ws->ngeometries;
    rt->noutputs = ws->noutputs;

    if (rt->nfilters > 0u) {
        rt->filters = (struct osh_transport_scoring_filter_runtime *) calloc(rt->nfilters, sizeof(*rt->filters));
        if (!rt->filters) {
            rc = OSH_ENOMEM;
            goto fail;
        }
        for (i = 0; i < rt->nfilters; ++i) {
            rc = copy_filter_runtime(&rt->filters[i], &ws->filters[i]);
            if (rc != OSH_OK) {
                goto fail;
            }
        }
    }

    if (rt->nsettings > 0u) {
        rt->settings =
            (struct osh_transport_scoring_settings_runtime *) calloc(rt->nsettings, sizeof(*rt->settings));
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
        rt->geometries =
            (struct osh_transport_scoring_geometry_runtime *) calloc(rt->ngeometries, sizeof(*rt->geometries));
        if (!rt->geometries) {
            rc = OSH_ENOMEM;
            goto fail;
        }
        for (i = 0; i < rt->ngeometries; ++i) {
            rc = copy_geometry_runtime(&rt->geometries[i], &ws->geometries[i]);
            if (rc != OSH_OK) {
                if (rc == OSH_EINVAL) {
                    osh_error("Scoring geometry '%s' has no valid runtime binning", ws->geometries[i].name);
                }
                goto fail;
            }
        }
    }

    geom_page_counts = (size_t *) calloc(rt->ngeometries ? rt->ngeometries : 1u, sizeof(*geom_page_counts));
    geom_next_page = (size_t *) calloc(rt->ngeometries ? rt->ngeometries : 1u, sizeof(*geom_next_page));
    if (!geom_page_counts || !geom_next_page) {
        rc = OSH_ENOMEM;
        goto fail;
    }

    total_pages = 0u;
    for (i = 0; i < ws->noutputs; ++i) {
        long gidx = find_geometry_index(ws, ws->outputs[i].geometry_name);

        if (gidx < 0) {
            osh_error("Scoring output '%s' references unknown geometry '%s'",
                      ws->outputs[i].filename ? ws->outputs[i].filename : "(unnamed)",
                      ws->outputs[i].geometry_name ? ws->outputs[i].geometry_name : "(null)");
            rc = OSH_EINVAL;
            goto fail;
        }
        geom_page_counts[gidx] += ws->outputs[i].npages;
        total_pages += ws->outputs[i].npages;
    }

    rt->npages = total_pages;
    if (rt->npages > 0u) {
        rt->pages = (struct osh_transport_scoring_page_runtime *) calloc(rt->npages, sizeof(*rt->pages));
        if (!rt->pages) {
            rc = OSH_ENOMEM;
            goto fail;
        }
    }
    if (rt->noutputs > 0u) {
        rt->outputs = (struct osh_transport_scoring_output_runtime *) calloc(rt->noutputs, sizeof(*rt->outputs));
        if (!rt->outputs) {
            rc = OSH_ENOMEM;
            goto fail;
        }
    }

    total_pages = 0u;
    for (i = 0; i < rt->ngeometries; ++i) {
        rt->geometries[i].first_page = total_pages;
        rt->geometries[i].npages = geom_page_counts[i];
        geom_next_page[i] = total_pages;
        total_pages += geom_page_counts[i];
    }

    for (i = 0; i < ws->noutputs; ++i) {
        long gidx = find_geometry_index(ws, ws->outputs[i].geometry_name);
        struct osh_transport_scoring_output_runtime *out = &rt->outputs[i];

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
        out->geometry_idx = (size_t) gidx;
        out->npages = ws->outputs[i].npages;
        if (out->npages > 0u) {
            out->page_indices = (size_t *) calloc(out->npages, sizeof(*out->page_indices));
            if (!out->page_indices) {
                rc = OSH_ENOMEM;
                goto fail;
            }
        }

        for (j = 0; j < ws->outputs[i].npages; ++j) {
            struct osh_scoring_page_def const *src_page = &ws->outputs[i].pages[j];
            struct osh_transport_scoring_page_runtime *dst_page;
            size_t page_idx = geom_next_page[gidx]++;

            out->page_indices[j] = page_idx;
            dst_page = &rt->pages[page_idx];
            dst_page->quantity = strdup(src_page->quantity);
            if (!dst_page->quantity) {
                rc = OSH_ENOMEM;
                goto fail;
            }
            dst_page->output_idx = i;
            dst_page->geometry_idx = (size_t) gidx;
            dst_page->len = rt->geometries[gidx].nbins;
            dst_page->data = (double *) calloc(dst_page->len ? dst_page->len : 1u, sizeof(*dst_page->data));
            if (!dst_page->data) {
                rc = OSH_ENOMEM;
                goto fail;
            }

            if (src_page->nfilter_names > 0u) {
                dst_page->filters = (struct osh_transport_scoring_page_filter_ref *) calloc(
                    src_page->nfilter_names, sizeof(*dst_page->filters));
                if (!dst_page->filters) {
                    rc = OSH_ENOMEM;
                    goto fail;
                }
                dst_page->nfilters = src_page->nfilter_names;
                for (size_t k = 0; k < src_page->nfilter_names; ++k) {
                    long fidx = find_filter_index(rt, src_page->filter_names[k]);

                    if (fidx < 0) {
                        osh_error("Scoring page '%s' references unknown filter '%s'",
                                  src_page->quantity ? src_page->quantity : "(unnamed)",
                                  src_page->filter_names[k]);
                        rc = OSH_EINVAL;
                        goto fail;
                    }
                    dst_page->filters[k].filter_idx = (size_t) fidx;
                }
            }
        }
    }

    free(geom_page_counts);
    free(geom_next_page);
    return OSH_OK;

fail:
    free(geom_page_counts);
    free(geom_next_page);
    osh_transport_scoring_runtime_free(rt);
    return rc;
}
