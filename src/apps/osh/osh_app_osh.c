#include "apps/osh/osh_app_osh.h"

#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_beam_parse.h"
#include "apps/osh/osh_beam_spotlist.h"
#include "apps/osh/osh_geometry_parse.h"
#include "apps/osh/osh_material_parse.h"
#include "apps/osh/osh_scoring_parse.h"
#include "common/osh_diag.h"
#include "common/osh_file.h"
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

        if (wb->shared.use_sad) {
            /* USECBEAM files follow the SH12A / DICOM RT Plan convention:
             * spot x/y are isocenter coordinates (lateral offset at z = 0 in
             * the beam-local PZALIGN frame).  osh_beam_spot.p[] must hold
             * physical beam-start coordinates, so back-project each spot to
             * the beam-entrance plane (spot->p[2] = BEAMPOS z, negative =
             * upstream).
             *
             * The projection from isocenter to beam start follows the line
             * from the virtual point source (at z = -SAD upstream of
             * isocenter) through the isocenter position:
             *
             *   x_bs = x_iso * (sad + z_start) / sad
             *
             * SAD is always positive (it is a distance, not a signed z
             * coordinate).  z_start is negative for a beam that enters
             * upstream of the isocenter, so (sad + z_start) < sad and
             * the beam-start offset is smaller than the isocenter offset.
             *
             * Example: x_iso = 5 cm, z_start = -50 cm, SAD = 200 cm
             *   factor = (200 - 50) / 200 = 0.75  →  x_bs = 3.75 cm
             * Without this step the SAD formula would compute an angle that
             * is too steep, making the field ~33% wider at isocenter. */
            size_t i;
            for (i = 0u; i < nspots; i++) {
                double z = spots[i].p[2]; /* beam-start z (negative = upstream) */
                spots[i].p[0] *= (wb->shared.sad[0] + z) / wb->shared.sad[0];
                spots[i].p[1] *= (wb->shared.sad[1] + z) / wb->shared.sad[1];
            }
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

    rc = osh_beam_workspace_prepare(wb, diag);
    if (rc != OSH_OK) {
        osh_beam_workspace_free(wb);
        return rc;
    }

    *wb_out = wb;
    return OSH_OK;
}

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
    if (rc != OSH_OK) {
        osh_material_workspace_free(wm);
        return rc;
    }

    rc = osh_material_workspace_prepare(wm, diag);
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

/* ---- Zone scoring name resolution ---------------------------------------- */

/**
 * @brief Is @p quantity finalised by dividing each bin by its volume?
 *
 * @details
 * The set whose estimator registers `postprocess_volume`/`postprocess_dosegy`
 * (see `src/scoring/runtime/osh_scoring_estimator.c`), i.e. the quantities for
 * which a missing per-zone `Volume` card silently changes the saved number.
 * `Energy` is extensive and `MCPL` is a phase-space dump — neither reads a bin
 * volume at all, so a Zone geometry carrying only those must not be warned
 * about (issue #328).  Quantity keywords are lowercased once at parse time.
 * `NKERMA` also postprocesses by volume but has no `detect.dat` spelling yet
 * (it is a registry placeholder, see docs/dev/scoring.md), so it is not listed.
 */
static int quantity_is_volume_normalized(char const *quantity) {
    if (!quantity) {
        return 0;
    }
    return strcmp(quantity, "fluence") == 0 || strcmp(quantity, "dose") == 0 || strcmp(quantity, "dosegy") == 0
           || strcmp(quantity, "dirtydose") == 0 || strcmp(quantity, "dirtydosegy") == 0;
}

/**
 * @brief Does any Output page attached to geometry @p gname need a bin volume?
 *
 * @details
 * Scans every Output that references the geometry by name — the same match
 * `osh_scoring_compile()` makes — and reports whether at least one of its pages
 * is volume-normalized.  A geometry no Output references yet also yields 0:
 * `osh_scoring_compile()` is what rejects that, and warning about a volume the
 * run will never read would be noise either way.
 */
static int geometry_needs_zone_volume(struct osh_scoring_workspace const *scoring, char const *gname) {
    size_t o;
    size_t p;

    if (!scoring || !gname) {
        return 0;
    }
    for (o = 0u; o < scoring->noutputs; ++o) {
        if (!scoring->outputs[o].geometry_name || strcmp(scoring->outputs[o].geometry_name, gname) != 0) {
            continue;
        }
        for (p = 0u; p < scoring->outputs[o].npages; ++p) {
            if (quantity_is_volume_normalized(scoring->outputs[o].pages[p].quantity)) {
                return 1;
            }
        }
    }
    return 0;
}

enum osh_status osh_scoring_resolve_zone_names(struct osh_scoring_workspace *scoring,
                                               struct osh_geometry_workspace const *geom,
                                               struct osh_diag_sink const *diag) {
    size_t g;

    if (!scoring || !geom) {
        return OSH_EINVAL;
    }

    for (g = 0u; g < scoring->ngeometries; ++g) {
        struct osh_scoring_geometry_def *geo;
        char const *gname;
        size_t i;
        int needs_volume;

        geo = &scoring->geometries[g];
        if (geo->nzone_indices == 0u || geo->zone_names == NULL) {
            continue; /* not a Zone geometry */
        }
        gname = geo->name;
        if (!gname) {
            gname = "(unnamed)";
        }
        needs_volume = geometry_needs_zone_volume(scoring, geo->name);

        /* Names -> dense 0-based transport indices; drop any earlier resolution. */
        free(geo->zone_indices);
        geo->zone_indices = (size_t *) calloc(geo->nzone_indices, sizeof(*geo->zone_indices));
        if (!geo->zone_indices) {
            return OSH_ENOMEM;
        }

        for (i = 0u; i < geo->nzone_indices; ++i) {
            size_t z;
            int found;

            found = 0;
            for (z = 0u; z < geom->nzones; ++z) {
                if (geom->zones[z].name && strcmp(geom->zones[z].name, geo->zone_names[i]) == 0) {
                    geo->zone_indices[i] = z;
                    found = 1;
                    break;
                }
            }
            if (!found) {
                OSH_DIAG_ERRORF(diag,
                                "Scoring Zone geometry '%s': unknown zone name '%s' (not defined in geo.dat)",
                                gname,
                                geo->zone_names[i]);
                return OSH_EINVAL;
            }
            if (needs_volume && !(geo->zone_volumes && geo->zone_volumes[i] > 0.0)) {
                OSH_DIAG_WARNF(diag,
                               "Scoring Zone geometry '%s': zone '%s' has no Volume card; "
                               "volume-normalized quantities (Fluence/Dose/DoseGy/DirtyDose/DirtyDoseGy) "
                               "use 1.0 cm3",
                               gname,
                               geo->zone_names[i]);
            }
        }
    }

    return OSH_OK;
}
