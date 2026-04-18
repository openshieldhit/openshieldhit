#include "apps/osh/osh_app_osh.h"

#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_beam_parse.h"
#include "apps/osh/osh_geometry_parse.h"
#include "apps/osh/osh_material_parse.h"
#include "apps/osh/osh_scoring_parse.h"
#include "beam/osh_beam_spots.h"
#include "openshieldhit/file.h"
#include "openshieldhit/logger.h"
#include "openshieldhit/material.h"
#include "openshieldhit/status.h"

enum osh_status osh_beam_setup_from_path(char const *path, struct osh_logger *lg, struct osh_beam_workspace **wb_out) {
    enum osh_status rc = OSH_OK;
    struct oshfile *sf = NULL;
    struct osh_beam_workspace *wb = NULL;
    char *spotlist_path = NULL;

    (void) lg;

    if (!path || !wb_out) {
        return OSH_EINVAL;
    }
    *wb_out = NULL;

    sf = osh_fopen(path);
    if (!sf) {
        return OSH_EIO;
    }

    rc = osh_beam_workspace_create(&wb);
    if (rc != OSH_OK) {
        osh_fclose(sf);
        return rc;
    }

    rc = osh_beam_parse(sf, wb, &spotlist_path);
    osh_fclose(sf);
    sf = NULL;
    if (rc != OSH_OK) {
        free(spotlist_path);
        osh_beam_workspace_free(wb);
        return rc;
    }

    if (spotlist_path) {
        osh_info("Loading spotlist before beam post-parse");
        rc = osh_beam_spotlist_load(wb, spotlist_path);
        if (rc != OSH_OK) {
            free(spotlist_path);
            osh_beam_workspace_free(wb);
            return rc;
        }
    }
    free(spotlist_path);

    rc = osh_beam_workspace_prepare(wb);
    if (rc != OSH_OK) {
        osh_beam_workspace_free(wb);
        return rc;
    }

    *wb_out = wb;
    return OSH_OK;
}

enum osh_status
osh_geometry_setup_from_path(char const *path, struct osh_logger *lg, struct osh_geometry_workspace **ws_out) {
    enum osh_status rc = OSH_OK;
    struct oshfile *sf = NULL;
    struct osh_geometry_workspace *ws = NULL;

    (void) lg;

    if (!path || !ws_out) {
        return OSH_EINVAL;
    }
    *ws_out = NULL;

    sf = osh_fopen(path);
    if (!sf) {
        return OSH_EIO;
    }

    rc = osh_geometry_workspace_create(&ws);
    if (rc != OSH_OK) {
        osh_fclose(sf);
        return rc;
    }

    rc = osh_geometry_parse(sf, ws);
    osh_fclose(sf);
    sf = NULL;
    if (rc != OSH_OK) {
        osh_geometry_workspace_free(ws);
        return rc;
    }

    rc = osh_geometry_workspace_prepare(ws);
    if (rc != OSH_OK) {
        osh_geometry_workspace_free(ws);
        return rc;
    }

    *ws_out = ws;
    return OSH_OK;
}

enum osh_status
osh_material_setup_from_path(char const *path, struct osh_logger *lg, struct osh_material_workspace **wm_out) {
    enum osh_status rc = OSH_OK;
    struct oshfile *sf = NULL;
    struct osh_material_workspace *wm = NULL;
    char *wdir = NULL;

    (void) lg;

    if (!path || !wm_out) {
        return OSH_EINVAL;
    }
    *wm_out = NULL;

    sf = osh_fopen(path);
    if (!sf) {
        return OSH_EIO;
    }

    rc = osh_material_workspace_create(&wm);
    if (rc != OSH_OK) {
        osh_fclose(sf);
        return rc;
    }

    wdir = osh_path_dirname(path);
    wm->wdir = wdir;
    wdir = NULL;
    wm->fname = strdup(path);
    if (!wm->fname) {
        osh_material_workspace_free(wm);
        osh_fclose(sf);
        return OSH_ENOMEM;
    }

    rc = osh_material_parse(sf, wm);
    osh_fclose(sf);
    sf = NULL;
    if (rc != OSH_OK) {
        osh_material_workspace_free(wm);
        return rc;
    }

    rc = osh_material_workspace_prepare(wm);
    if (rc != OSH_OK) {
        osh_material_workspace_free(wm);
        return rc;
    }

    *wm_out = wm;
    return OSH_OK;
}

enum osh_status
osh_scoring_setup_from_path(char const *path, struct osh_logger *lg, struct osh_scoring_workspace **ws_out) {
    enum osh_status rc = OSH_OK;
    struct oshfile *sf = NULL;
    struct osh_scoring_workspace *ws = NULL;

    (void) lg;

    if (!path || !ws_out) {
        return OSH_EINVAL;
    }
    *ws_out = NULL;

    sf = osh_fopen(path);
    if (!sf) {
        return OSH_EIO;
    }

    rc = osh_scoring_workspace_create(&ws);
    if (rc != OSH_OK) {
        osh_fclose(sf);
        return rc;
    }

    ws->fname = strdup(path);
    if (!ws->fname) {
        osh_fclose(sf);
        osh_scoring_workspace_free(ws);
        return OSH_ENOMEM;
    }

    rc = osh_scoring_parse(sf, ws);
    osh_fclose(sf);
    sf = NULL;
    if (rc != OSH_OK) {
        osh_scoring_workspace_free(ws);
        return rc;
    }

    *ws_out = ws;
    return OSH_OK;
}
