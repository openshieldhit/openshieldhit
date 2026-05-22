#include "scoring/osh_scoring.h"

#include <stdlib.h>
#include <string.h>

enum osh_status osh_scoring_workspace_create(struct osh_scoring_workspace **ws_out) {
    struct osh_scoring_workspace *ws;

    if (!ws_out) {
        return OSH_EINVAL;
    }
    ws = (struct osh_scoring_workspace *) calloc(1, sizeof(*ws));
    if (!ws) {
        return OSH_ENOMEM;
    }
    *ws_out = ws;
    return OSH_OK;
}

enum osh_status osh_scoring_workspace_free(struct osh_scoring_workspace *ws) {
    size_t i;
    size_t j;

    if (!ws) {
        return OSH_OK;
    }

    free(ws->fname);

    for (i = 0; i < ws->nfilters; ++i) {
        free(ws->filters[i].name);
        free(ws->filters[i].rules);
    }
    free(ws->filters);

    for (i = 0; i < ws->nsettings; ++i) {
        free(ws->settings[i].name);
        free(ws->settings[i].material_name);
    }
    free(ws->settings);

    for (i = 0; i < ws->ngeometries; ++i) {
        free(ws->geometries[i].kind);
        free(ws->geometries[i].name);
        free(ws->geometries[i].axes);
        free(ws->geometries[i].vox_rtdose_path);
        free(ws->geometries[i].vox_body_name);
    }
    free(ws->geometries);

    for (i = 0; i < ws->noutputs; ++i) {
        free(ws->outputs[i].filename);
        free(ws->outputs[i].geometry_name);
        free(ws->outputs[i].fileformat);
        for (j = 0; j < ws->outputs[i].npages; ++j) {
            size_t k;

            free(ws->outputs[i].pages[j].quantity);
            for (k = 0; k < ws->outputs[i].pages[j].nfilter_names; ++k) {
                free(ws->outputs[i].pages[j].filter_names[k]);
            }
            free(ws->outputs[i].pages[j].filter_names);
        }
        free(ws->outputs[i].pages);
    }
    free(ws->outputs);

    free(ws);
    return OSH_OK;
}

struct osh_scoring_filter_def const *osh_scoring_filter_by_name(struct osh_scoring_workspace const *ws,
                                                                char const *name) {
    size_t i;

    if (!ws || !name) {
        return NULL;
    }
    for (i = 0; i < ws->nfilters; ++i) {
        if (ws->filters[i].name && strcmp(ws->filters[i].name, name) == 0) {
            return &ws->filters[i];
        }
    }
    return NULL;
}

struct osh_scoring_settings_def const *osh_scoring_settings_by_name(struct osh_scoring_workspace const *ws,
                                                                    char const *name) {
    size_t i;

    if (!ws || !name) {
        return NULL;
    }
    for (i = 0; i < ws->nsettings; ++i) {
        if (ws->settings[i].name && strcmp(ws->settings[i].name, name) == 0) {
            return &ws->settings[i];
        }
    }
    return NULL;
}

struct osh_scoring_geometry_def const *osh_scoring_geometry_by_name(struct osh_scoring_workspace const *ws,
                                                                    char const *name) {
    size_t i;

    if (!ws || !name) {
        return NULL;
    }
    for (i = 0; i < ws->ngeometries; ++i) {
        if (ws->geometries[i].name && strcmp(ws->geometries[i].name, name) == 0) {
            return &ws->geometries[i];
        }
    }
    return NULL;
}

struct osh_scoring_output_def const *osh_scoring_output_by_filename(struct osh_scoring_workspace const *ws,
                                                                    char const *filename) {
    size_t i;

    if (!ws || !filename) {
        return NULL;
    }
    for (i = 0; i < ws->noutputs; ++i) {
        if (ws->outputs[i].filename && strcmp(ws->outputs[i].filename, filename) == 0) {
            return &ws->outputs[i];
        }
    }
    return NULL;
}
