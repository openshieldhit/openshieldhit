#include <stdio.h>
#include <stdlib.h>

#include "osh_beam.h"
#include "osh_beam_parse.h"
#include "osh_beam_spots.h"
#include "osh_file.h"

#include "osh_logger.h"
#include "osh_rc.h"

static void _wb_defaults(struct beam_workspace *wb);
static int _wb_validate(const struct beam_workspace *wb);

int osh_beam_setup(const char *filename, const char *wdir, struct beam_workspace **wb_out) {
    int rc = OSH_OK;
    struct oshfile *sf = NULL;
    struct beam_workspace *wb = NULL;

    if (!filename || !wb_out)
        return OSH_EINVAL;
    *wb_out = NULL;

    sf = osh_fopen(filename);
    if (!sf)
        return OSH_EIO;

    wb = (struct beam_workspace *)calloc(1, sizeof *wb);
    if (!wb) {
        osh_fclose(sf);
        return OSH_ENOMEM;
    }

    rc = osh_beam_spots_init(&wb->spots, 1);
    if (rc != OSH_OK) {
        osh_fclose(sf);
        free(wb);
        return rc;
    }

    rc = osh_beam_shared_init(wb->shared);
    if (rc != OSH_OK) {
        osh_fclose(sf);
        osh_beam_spots_free(wb->spots);
        free(wb);
        return rc;
    }

    _wb_defaults(wb);

    if (wdir && *wdir) {
        size_t L = strlen(wdir);
        wb->wdir = (char *)malloc(L + 1);
        if (!wb->wdir) {
            osh_fclose(sf);
            osh_beam_workspace_free(wb);
            return OSH_ENOMEM;
        }
        memcpy(wb->wdir, wdir, L + 1);
    }

    /* Parse: fills wb->shared, wb->spotlist (and spot->_tm), or wb->phsp, sets
     * wb->beam_mode */
    rc = osh_beam_parse(sf, wb);
    osh_fclose(sf);
    sf = NULL;
    if (rc != OSH_OK) {
        /* Optional: a non-fatal log could show diag.message/line/col here */
        osh_beam_workspace_free(wb);
        return rc;
    }

    rc = _wb_validate(wb);
    if (rc != OSH_OK) {
        osh_beam_workspace_free(wb);
        return rc;
    }

    *wb_out = wb;
    return OSH_OK;
}

int osh_beam_workspace_free(struct beam_workspace *wb) {
    if (!wb)
        return OSH_OK;
    if (wb->wdir)
        free(wb->wdir);
    if (wb->spots) {
        osh_beam_spots_free(wb->spots);
    }
    if (wb->shared) {
        osh_beam_shared_free(wb->shared);
    }
    if (wb->phsp) {
        // TODO osh_beam_phsp_free(wb->phsp);
        free(wb->phsp);
    }
    if (wb->rifi) {
        // TODO osh_beam_rifi_free(wb->rifi);
        free(wb->rifi);
    }
    if (wb->parlev) {
        // TODO osh_beam_parlev_free(wb->parlev);
        free(wb->parlev);
    }
    if (wb->fname)
        free(wb->fname);
    free(wb);
    return OSH_OK;
}

static void _wb_defaults(struct beam_workspace *wb) {
    if (!wb)
        return;

    wb->wdir = NULL;
    wb->fname = NULL;
    wb->fname_spotlist = NULL;
    wb->spots = NULL;
    wb->shared = NULL;
    wb->phsp = NULL;
    wb->rifi = NULL;
    wb->parlev = NULL;

    wb->nspots = 0;
    wb->nstat = 0;
    wb->nsave = 0;
    wb->rndseed = 0;
    wb->rndoffset = 0;
    wb->deltae = 0.0f;
    wb->oln = 0.0f;
    wb->demin = 0.0f;
    wb->straggl = 0;
    wb->scatter = 0;
    wb->nuclear = 0;
    wb->emtrans = 0;
    wb->apcorr = 0;
    wb->beam_mode = 0;
    wb->makeln = 0;
    wb->neutrfast = 0;
}

static int _wb_validate(const struct beam_workspace *wb) {
    if (!wb)
        return OSH_EINVAL;

    // Check required fields
    if (!wb->spots && !wb->phsp)
        return OSH_EINVAL;
    if (wb->spots && wb->nspots == 0)
        return OSH_EINVAL;

    return OSH_OK;
}