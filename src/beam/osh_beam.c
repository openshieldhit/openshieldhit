#include "beam/osh_beam.h"

#include <stdio.h>
#include <stdlib.h>

#include "beam/osh_beam_parse.h"
#include "beam/osh_beam_spots.h"
#include "common/osh_file.h"
#include "common/osh_logger.h"
#include "common/osh_rc.h"

static void _wb_defaults(struct beam_workspace *wb);
static int _wb_validate(const struct beam_workspace *wb);

int osh_beam_setup_from_path(char const *path, struct osh_logger *lg, struct beam_workspace **wb_out) {
    int rc = OSH_OK;
    struct oshfile *sf = NULL;
    struct beam_workspace *wb = NULL;
    char *wdir = NULL;

    /* lg may be NULL — callers pass NULL to use the global default logger.
     * Internal code should call osh_log_default() when lg is NULL. */
    (void) lg; /* TODO: thread through to parser diagnostics */

    if (!path || !wb_out) {
        return OSH_EINVAL;
    }
    *wb_out = NULL;

    /* Derive working directory. Path must be pre-normalized to '/' by the
     * caller (use osh_path_normalize() in main.c before any library calls). */
    wdir = osh_path_dirname(path);

    sf = osh_fopen(path);
    if (!sf) {
        free(wdir);
        return OSH_EIO;
    }

    wb = (struct beam_workspace *) calloc(1, sizeof *wb);
    if (!wb) {
        osh_fclose(sf);
        free(wdir);
        return OSH_ENOMEM;
    }

    _wb_defaults(wb);
    wb->wdir = wdir; /* ownership transferred */
    wdir = NULL;

    rc = osh_beam_spots_init(&wb->spots, 1);
    if (rc != OSH_OK) {
        osh_fclose(sf);
        osh_beam_workspace_free(wb);
        return rc;
    }

    rc = osh_beam_shared_init(&wb->shared);
    if (rc != OSH_OK) {
        osh_fclose(sf);
        osh_beam_workspace_free(wb);
        return rc;
    }

    rc = osh_beam_parse(sf, wb);
    osh_fclose(sf);
    sf = NULL;
    if (rc != OSH_OK) {
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
    if (!wb) {
        return OSH_OK;
    }
    if (wb->wdir) {
        free(wb->wdir);
    }
    if (wb->spots) {
        osh_beam_spots_free(wb->spots);
    }
    /* wb->shared is embedded by value — no free needed */
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
    if (wb->fname) {
        free(wb->fname);
    }
    free(wb);
    return OSH_OK;
}

static void _wb_defaults(struct beam_workspace *wb) {
    if (!wb) {
        return;
    }

    wb->wdir = NULL;
    wb->fname = NULL;
    wb->fname_spotlist = NULL;
    wb->spots = NULL;
    /* wb->shared is embedded by value (see osh_beam.h); zero-initialised by
     * calloc above, then given sensible defaults by osh_beam_shared_init(). */
    wb->phsp = NULL;
    wb->rifi = NULL;
    wb->parlev = NULL;

    wb->nspots = 0;
    wb->nstat = 0;
    wb->nsave = 0;
    wb->rndseed = 0;
    wb->rndoffset = 0;
    wb->deltae = 0.0f;
    wb->ncut = 0.0f;
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
    if (!wb) {
        return OSH_EINVAL;
    }

    // Check required fields
    if (!wb->spots && !wb->phsp) {
        return OSH_EINVAL;
    }
    if (wb->spots && wb->nspots == 0) {
        return OSH_EINVAL;
    }

    return OSH_OK;
}

/* TODO: wire osh_beam_print / osh_beam_print_spot to the logger (osh_logger)
 * rather than stdout once the logger is threaded through the workspace. */
void osh_beam_print(struct beam_workspace const *wb) {
    if (!wb) {
        return;
    }
    printf("beam_workspace: nspots=%zu nstat=%zu mode=%d\n", wb->nspots, wb->nstat, (int) wb->beam_mode);
}

void osh_beam_print_spot(struct beam_spot const *spot) {
    if (!spot) {
        return;
    }
    printf("beam_spot: t0=%.4f MeV p=[%.2f %.2f %.2f] cm id=%u layer=%u\n",
           spot->t0,
           spot->p[0],
           spot->p[1],
           spot->p[2],
           spot->spot_id,
           spot->layer_id);
}
