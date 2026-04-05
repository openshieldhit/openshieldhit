#include "beam/osh_beam.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "beam/osh_beam_parse.h"
#include "beam/osh_beam_spots.h"
#include "common/osh_const.h"
#include "common/osh_coord.h"
#include "common/osh_file.h"
#include "common/osh_logger.h"
#include "common/osh_physics.h"
#include "common/osh_rc.h"
#include "common/osh_vect.h"

static void _wb_defaults(struct beam_workspace *wb);
static int _wb_validate(const struct beam_workspace *wb);
static int _wb_postparse(struct beam_workspace *wb);

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
    wb->nspots = 1;

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

    if (wb->fname_spotlist) {
        osh_info("Loading spotlist before beam post-parse");
        rc = osh_beam_spotlist_load(wb);
        if (rc != OSH_OK) {
            osh_beam_workspace_free(wb);
            return rc;
        }
    }

    rc = _wb_validate(wb);
    if (rc != OSH_OK) {
        osh_beam_workspace_free(wb);
        return rc;
    }

    rc = _wb_postparse(wb);
    if (rc != OSH_OK) {
        osh_beam_workspace_free(wb);
        return rc;
    }

    osh_beam_print(wb);

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
    if (wb->cum_wt) {
        free(wb->cum_wt);
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
    if (wb->fname_spotlist) {
        free(wb->fname_spotlist);
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
    wb->cum_wt = NULL;
    /* wb->shared is embedded by value (see osh_beam.h); zero-initialised by
     * calloc above, then given sensible defaults by osh_beam_shared_init(). */
    wb->phsp = NULL;
    wb->rifi = NULL;
    wb->parlev = NULL;
    wb->has_primary = 0;

    wb->nspots = 0;
    wb->wt_sum = 0.0;
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

    if (!wb->spots && !wb->phsp) {
        return OSH_EINVAL;
    }
    if (wb->spots && wb->nspots == 0) {
        return OSH_EINVAL;
    }

    return OSH_OK;
}

/* Derive missing t0/p0 and tsigma/psigma for one spot.
 *
 * The parser stores either t0 (positive input) or p0 (negative input, fabs).
 * Post-parse fills in whichever was not given, using the particle rest mass.
 * psigma is derived from tsigma via the first-order Jacobian dp/dT = E/p,
 * i.e. psigma = tsigma * (t0 + mass) / p0, valid when sigma << T0.
 *
 * Skipped silently when no primary particle was resolved yet. */
static void _postparse_spot_energy(struct beam_spot *spot) {
    double mass;

    if (!spot->part) {
        return;
    }

    mass = spot->part->mass;

    if (spot->t0 > 0.0 && spot->p0 == 0.0) {
        spot->p0 = osh_physics_momentum(spot->t0, mass);
    } else if (spot->p0 > 0.0 && spot->t0 == 0.0) {
        spot->t0 = osh_physics_tkin(spot->p0, mass);
    }

    /* first-order sigma conversion: dp/dT = E/p,  dT/dp = p/E */
    if (spot->tsigma > 0.0 && spot->psigma == 0.0 && spot->p0 > 0.0) {
        spot->psigma = spot->tsigma * (spot->t0 + mass) / spot->p0;
    } else if (spot->psigma > 0.0 && spot->tsigma == 0.0 && spot->p0 > 0.0) {
        spot->tsigma = spot->psigma * spot->p0 / (spot->t0 + mass);
    }
}

/* Build the beam-sampling affine matrix for one spot using the standard
 * local->universe helper in osh_vect.
 *
 * Matrix contract (row-major):
 *   p_universe = R * p_local + R * spot->p
 *
 * spot->p is interpreted as a local beam-frame offset before gantry/table
 * rotation. */
static void _build_spot_tm(struct beam_spot *spot, const struct beam_shared *sh) {
    double cs[3];
    double r[3];

    if (fabs(sh->theta) > 0.0) {
        cs[0] = cos(sh->theta);
        cs[1] = sin(sh->phi);
        cs[2] = cos(sh->phi);
    } else {
        cs[0] = 1.0; /* cos(0) */
        cs[1] = 0.0; /* sin(0) */
        cs[2] = 1.0; /* cos(0) */
    }

    osh_coord_c2v(cs, r);
    osh_vect_setup_tmatrix_bzalign_affine(spot->p, r, spot->_tm);
}

static int _wb_postparse(struct beam_workspace *wb) {
    size_t i;
    double wt_acc;

    if (!wb) {
        return OSH_EINVAL;
    }

    wb->shared.emax = 0.0;
    wb->shared.pmax = 0.0;
    wb->wt_sum = 0.0;

    free(wb->cum_wt);
    wb->cum_wt = NULL;

    if (wb->nspots > 0) {
        wb->cum_wt = (double *) calloc(wb->nspots, sizeof(double));
        if (!wb->cum_wt) {
            return OSH_ENOMEM;
        }
    }

    wt_acc = 0.0;
    for (i = 0; i < wb->nspots; i++) {
        if (wb->has_primary) {
            wb->spots[i].part = &wb->primary;
        }
        _postparse_spot_energy(&wb->spots[i]);
        _build_spot_tm(&wb->spots[i], &wb->shared);

        if (wb->spots[i].wt < 0.0) {
            return OSH_EINVAL;
        }
        wt_acc += wb->spots[i].wt;
        wb->cum_wt[i] = wt_acc;

        if (wb->spots[i].t0 > wb->shared.emax) {
            wb->shared.emax = wb->spots[i].t0;
        }
        if (wb->spots[i].p0 > wb->shared.pmax) {
            wb->shared.pmax = wb->spots[i].p0;
        }
    }

    wb->wt_sum = wt_acc;

    return OSH_OK;
}

/* osh_beamdef.h name arrays (osh_beam_mscat_names etc.) included via osh_beam.h */


void osh_beam_print(struct beam_workspace const *wb) {
    if (!wb) {
        return;
    }

    osh_info("%s", "");
    osh_info("%s", "");
    osh_info("Beam configuration:");
    osh_info(OSH_LOG_HLINE);
    osh_info("%-40s : %li", "Requested primaries NSTAT", (long int) wb->nstat);
    if (wb->nsave > 0) {
        osh_info("%-40s : %li", "Save interval NSAVE", (long int) wb->nsave);
    } else {
        osh_info("%-40s : %s", "Save interval NSAVE", "OFF");
    }
    osh_info("%-40s : %i", "Random seed RNDSEED", wb->rndseed);
    osh_info("%-40s : %i", "Random seed offset", wb->rndoffset);
    osh_info("%s", "");
    osh_info("%-18s : %f", "DeltaE/E", wb->deltae);
    osh_info("%-18s : %f MeV", "Neutron cut", wb->ncut);
    osh_info("%-18s : %f MeV/n", "DeltaE min", wb->demin);
    osh_info("%s", "");
    osh_info("%-18s : %s", "Scatter mode", osh_beam_mscat_names[(int) wb->scatter]);
    osh_info("%-18s : %s", "Straggling mode", osh_beam_stragg_names[(int) wb->straggl]);
    osh_info("%-18s : %s", "Nuclear react.", osh_log_offon[(int) wb->nuclear]);
    osh_info("%-18s : %s", "Apcorr mode", osh_log_offon[(int) wb->apcorr]);
    osh_info("%-18s : %s", "Beam mode", osh_beam_mode_names[(int) wb->beam_mode]);
    osh_info("%-18s : %s", "Make LN", osh_log_offon[(int) wb->makeln]);
    osh_info("%-18s : %s", "Fast neutrons", osh_log_offon[(int) wb->neutrfast]);
    osh_info("%s", "");
    osh_info("%-18s : %.3f  %.3f  cm", "SAD (x,y)", wb->shared.sad[0], wb->shared.sad[1]);
    osh_info("%-18s : %.3f  cm", "Focus", wb->shared.focus);
    osh_info("%-18s : %.3f  deg", "Theta", wb->shared.theta * 180.0 * OSH_M_1_PI);
    osh_info("%-18s : %.3f  deg", "Phi", wb->shared.phi * 180.0 * OSH_M_1_PI);
    osh_info("%-18s : %.3f  MeV", "Emax", wb->shared.emax);
    osh_info("%-18s : %.3f  MeV/c", "Pmax", wb->shared.pmax);

    if (wb->spots) {
        osh_beam_print_spot(&wb->spots[0]);
    }
}

void osh_beam_print_spot(struct beam_spot const *spot) {
    if (!spot) {
        return;
    }

    if (spot->part) {
        osh_info("%s", "");
        osh_print_particle(spot->part);
    }

    osh_info("%s", "");
    osh_info("%-18s : %.3f  %.3f  %.3f  cm", "Position", spot->p[0], spot->p[1], spot->p[2]);
    osh_info("%-18s : %.3f  %.3f  cm", "Size/sigma", spot->size[0], spot->size[1]);
    osh_info("%-18s : %.3f  %.3f  mrad", "Divergence", spot->div[0] * 1000.0, spot->div[1] * 1000.0);
    osh_info("%-18s : %.3f  %.3f", "Correlation", spot->cor[0], spot->cor[1]);
    osh_info("%s", "");
    osh_info("%-18s : %.3f  MeV", "T0", spot->t0);
    osh_info("%-18s : %.3f  MeV", "TSigma", spot->tsigma);
    osh_info("%-18s : %.3f  MeV/c", "P0", spot->p0);
    osh_info("%-18s : %.3f  MeV/c", "PSigma", spot->psigma);
    osh_info("%-18s : %i", "TSigma type", (int) spot->tsigma_type);
    osh_info("%s", "");
    osh_info("%-18s : %.3f", "Stat.weight", spot->wt);
    osh_info("%-18s : %u", "Spot ID", spot->spot_id);
    osh_info("%-18s : %u", "Layer ID", spot->layer_id);
    if (spot->shape >= 0 && (int) spot->shape < 4) {
        osh_info("%-18s : %s", "Shape", osh_beam_shape_names[(int) spot->shape]);
    } else {
        osh_info("%-18s : invalid (%i)", "Shape", (int) spot->shape);
    }
}
