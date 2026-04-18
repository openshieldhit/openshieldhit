#include "apps/osh/osh_app_osh.h"

#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_beam_parse.h"
#include "apps/osh/osh_beam_spotlist.h"
#include "apps/osh/osh_geometry_parse.h"
#include "apps/osh/osh_material_parse.h"
#include "apps/osh/osh_scoring_parse.h"
#include "common/osh_file.h"
#include "common/osh_logger.h"
#include "openshieldhit/material.h"
#include "openshieldhit/status.h"

/**
 * @brief Load and finalize a beam workspace from a legacy OpenShieldHIT beam file.
 *
 * @details
 * The app layer owns all file I/O and format-specific policy. This function
 * parses `beam.dat`, accumulates one template spot from inline beam cards,
 * optionally imports an external spotlist referenced by `USECBEAM`, then
 * populates the public beam workspace exclusively through `osh_beam_spots_set`
 * before calling `osh_beam_workspace_prepare()`.
 */
enum osh_status
osh_beam_setup_from_path(char const *path, struct osh_diag_sink const *diag, struct osh_beam_workspace **wb_out) {
    enum osh_status rc = OSH_OK;
    struct oshfile *sf = NULL;
    struct osh_beam_workspace *wb = NULL;
    struct osh_beam_spot template_spot;
    struct osh_beam_spot *spots = NULL;
    size_t nspots = 0u;
    char *spotlist_path = NULL;

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

    memset(&template_spot, 0, sizeof(template_spot));
    template_spot.shape = OSH_BEAM_SHAPE_PENCIL;

    rc = osh_beam_parse(sf, diag, wb, &template_spot, &spotlist_path);
    osh_fclose(sf);
    sf = NULL;
    if (rc != OSH_OK) {
        free(spotlist_path);
        osh_beam_workspace_free(wb);
        return rc;
    }

    if (spotlist_path) {
        OSH_DIAG_INFOF(diag, "Loading spotlist before beam post-parse");
        rc = osh_beam_spotlist_import(spotlist_path, diag, &template_spot, &spots, &nspots);
        if (rc != OSH_OK) {
            free(spotlist_path);
            osh_beam_workspace_free(wb);
            return rc;
        }
        rc = osh_beam_spots_set(wb, spots, nspots);
        free(spots);
        if (rc != OSH_OK) {
            free(spotlist_path);
            osh_beam_workspace_free(wb);
            return rc;
        }
    } else {
        rc = osh_beam_spots_set(wb, &template_spot, 1u);
        if (rc != OSH_OK) {
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

/**
 * @brief Load and finalize a geometry workspace from a legacy OpenShieldHIT geometry file.
 *
 * @details
 * This is the file-oriented wrapper around `osh_geometry_parse()` plus the
 * required geometry workspace preparation step.
 */
enum osh_status osh_geometry_setup_from_path(char const *path,
                                             struct osh_diag_sink const *diag,
                                             struct osh_geometry_workspace **ws_out) {
    enum osh_status rc = OSH_OK;
    struct oshfile *sf = NULL;
    struct osh_geometry_workspace *ws = NULL;

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

    rc = osh_geometry_parse(sf, diag, ws);
    osh_fclose(sf);
    sf = NULL;
    if (rc != OSH_OK) {
        osh_geometry_workspace_free(ws);
        return rc;
    }

    rc = osh_geometry_workspace_prepare(ws, diag);
    if (rc != OSH_OK) {
        osh_geometry_workspace_free(ws);
        return rc;
    }

    *ws_out = ws;
    return OSH_OK;
}

/**
 * @brief Load and finalize a material workspace from a legacy OpenShieldHIT material file.
 *
 * @details
 * The app wrapper still records the source path and base directory on the
 * material workspace because `LOADDEDX` path resolution remains file-relative
 * during parsing. The parser itself owns file I/O for imported tables and
 * converts them into material-owned in-memory overrides.
 */
enum osh_status osh_material_setup_from_path(char const *path,
                                             struct osh_diag_sink const *diag,
                                             struct osh_material_workspace **wm_out) {
    enum osh_status rc = OSH_OK;
    struct oshfile *sf = NULL;
    struct osh_material_workspace *wm = NULL;
    char *wdir = NULL;

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

    rc = osh_material_parse(sf, diag, wm);
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

/**
 * @brief Load a scoring workspace from a legacy OpenShieldHIT detect file.
 *
 * @details
 * Scoring has no separate workspace prepare step at present, so this wrapper
 * only allocates the workspace, records the source filename, and parses the
 * cold scoring definitions from `detect.dat`.
 */
enum osh_status
osh_scoring_setup_from_path(char const *path, struct osh_diag_sink const *diag, struct osh_scoring_workspace **ws_out) {
    enum osh_status rc = OSH_OK;
    struct oshfile *sf = NULL;
    struct osh_scoring_workspace *ws = NULL;

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

    rc = osh_scoring_parse(sf, diag, ws);
    osh_fclose(sf);
    sf = NULL;
    if (rc != OSH_OK) {
        osh_scoring_workspace_free(ws);
        return rc;
    }

    *ws_out = ws;
    return OSH_OK;
}
