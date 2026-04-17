#include "beam/osh_beam.h"

#include <math.h>
#include <stdlib.h>

#include "beam/osh_beam_prepared.h"
#include "beam/osh_beam_spots.h"
#include "common/osh_const.h"
#include "common/osh_coord.h"
#include "common/osh_logger.h"
#include "common/osh_physics.h"
#include "common/osh_rc.h"
#include "common/osh_vect.h"
#include "particle/osh_particle.h"

static void _wb_defaults(struct beam_workspace *wb);
static int _wb_validate(const struct beam_workspace *wb);
static int _wb_postparse(struct beam_workspace *wb);
static int _resolve_primary_particle(struct beam_workspace const *wb, struct particle *part_out);
static void _postparse_spot_energy(struct beam_spot *spot, struct particle const *part);
static void _build_spot_tm(struct beam_spot const *spot, struct beam_shared const *sh, double tm_out[16]);

enum osh_status osh_beam_workspace_create(struct beam_workspace **wb_out) {
    enum osh_status rc = OSH_OK;
    struct beam_workspace *wb = NULL;

    if (!wb_out) {
        return OSH_EINVAL;
    }
    *wb_out = NULL;

    wb = (struct beam_workspace *) calloc(1, sizeof *wb);
    if (!wb) {
        return OSH_ENOMEM;
    }

    _wb_defaults(wb);

    rc = osh_beam_spots_init(&wb->spots, 1);
    if (rc != OSH_OK) {
        osh_beam_workspace_free(wb);
        return rc;
    }
    wb->nspots = 1;

    rc = osh_beam_shared_init(&wb->shared);
    if (rc != OSH_OK) {
        osh_beam_workspace_free(wb);
        return rc;
    }

    *wb_out = wb;
    return OSH_OK;
}

enum osh_status osh_beam_workspace_prepare(struct beam_workspace *wb) {
    enum osh_status rc;

    rc = _wb_validate(wb);
    if (rc != OSH_OK) {
        return rc;
    }

    rc = _wb_postparse(wb);
    if (rc != OSH_OK) {
        return rc;
    }

    if (osh_log_get_level() <= OSH_LOG_INFO) {
        osh_beam_print(wb);
    }

    return OSH_OK;
}

enum osh_status osh_beam_workspace_free(struct beam_workspace *wb) {
    if (!wb) {
        return OSH_OK;
    }
    if (wb->spots) {
        osh_beam_spots_free(wb->spots);
    }
    if (wb->prepared) {
        free(wb->prepared->cum_wt);
        free(wb->prepared->tm);
        free(wb->prepared);
    }
    /* wb->shared is embedded by value — no free needed */
    if (wb->phsp) {
        osh_beam_phsp_free(wb->phsp);
    }
    if (wb->rifi) {
        /* TODO: replace with osh_beam_rifi_free() once struct ripple_filter is defined */
        free(wb->rifi);
    }
    if (wb->parlev) {
        /* TODO: replace with osh_beam_parlev_free() once struct parlev is defined */
        free(wb->parlev);
    }
    free(wb);
    return OSH_OK;
}

enum osh_status osh_beam_phsp_free(struct beam_phsp *phsp) {
    int axis;

    if (!phsp) {
        return OSH_OK;
    }
    for (axis = 0; axis < 3; axis++) {
        free(phsp->p[axis]);
        free(phsp->d[axis]);
    }
    free(phsp->e);
    free(phsp->wt);
    free(phsp->pdg);
    free(phsp->fname);
    free(phsp);
    return OSH_OK;
}

/**
 * @brief Set all pointer and scalar fields in a freshly-allocated workspace to
 *        safe zero/NULL values.
 *
 * @details
 * Called immediately after calloc(), which already zero-initialises memory.
 * The explicit assignments here serve as documentation of every field the
 * workspace owns and make the intended initial state clear to the reader,
 * independent of the allocator used.
 *
 * @param[in,out] wb  Workspace to initialise; silently ignored when NULL.
 */
static void _wb_defaults(struct beam_workspace *wb) {
    if (!wb) {
        return;
    }

    wb->spots = NULL;
    /* wb->shared is embedded by value (see osh_beam.h); zero-initialised by
     * calloc above, then given sensible defaults by osh_beam_shared_init(). */
    wb->phsp = NULL;
    wb->prepared = NULL;
    wb->rifi = NULL;
    wb->parlev = NULL;
    wb->has_primary = 0;

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

/**
 * @brief Check that the workspace is in a consistent state after parsing.
 *
 * @details
 * Verifies that exactly one source type (spots or phsp) is present, that the
 * spot array is non-empty when spots mode is active, and that a PRIMARY
 * particle was resolved for spot-based sources.
 *
 * @param[in] wb  Workspace to validate.
 *
 * @returns OSH_OK if valid, OSH_EINVAL otherwise.
 */
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
    if (wb->spots && !wb->has_primary) {
        osh_warn("Beam setup requires PRIMARY for spot-based sources");
        return OSH_EINVAL;
    }

    return OSH_OK;
}

/**
 * @brief Derive whichever of t0/p0 and tsigma/psigma was not given by the parser.
 *
 * @details
 * The parser stores either t0 (positive TMAX0 input or spot-list energy) or p0
 * (negative TMAX0 input, stored as fabs). Positive user-facing ion energies
 * are specific energies; this function first converts marked t0 and tsigma values
 * to total MeV by multiplying by A for ions with A > 1. It then fills in the
 * missing quantity using the particle rest mass via the relativistic
 * energy-momentum relation:
 *
 *   p = sqrt(T² + 2*T*m)    (momentum from kinetic energy)
 *   T = sqrt(p² + m²) - m   (kinetic energy from momentum)
 *
 * The momentum spread is derived from the energy spread using the first-order
 * Jacobian dp/dT = E/p = (T + m)/p, which gives:
 *
 *   psigma = tsigma * (t0 + mass) / p0
 *
 * This linear approximation is valid when sigma << T0.
 * Silently skipped when @p part is NULL (PRIMARY not yet resolved).
 *
 * @param[in,out] spot  Beam spot whose energy/momentum fields are completed.
 * @param[in]     part  Resolved primary particle species.
 */
static void _postparse_spot_energy(struct beam_spot *spot, struct particle const *part) {
    double mass;
    double scale;

    if (!spot || !part) {
        return;
    }

    mass = part->mass;

    if (spot->t0_per_nucleon || spot->tsigma_per_nucleon) {
        scale = 1.0;
        if (part->is_nucleus && part->a > 1u) {
            scale = (double) part->a;
        }
        if (spot->t0_per_nucleon) {
            spot->t0 *= scale;
            spot->t0_per_nucleon = 0;
        }
        if (spot->tsigma_per_nucleon) {
            spot->tsigma *= scale;
            spot->tsigma_per_nucleon = 0;
        }
    }

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

static int _resolve_primary_particle(struct beam_workspace const *wb, struct particle *part_out) {
    if (!wb || !part_out || !wb->has_primary) {
        return OSH_EINVAL;
    }
    if (!osh_particle_from_pdg(part_out, wb->primary.pdg)) {
        return OSH_EINVAL;
    }
    return OSH_OK;
}

/**
 * @brief Build the affine sampling matrix spot->_tm for one beam spot.
 *
 * @details
 * Converts the beam direction (theta, phi) to a unit direction vector r,
 * then calls osh_vect_setup_tmatrix_bzalign_affine() to build a 4x4
 * row-major matrix that maps beam-local PZALIGN coordinates to UNIVERSE:
 *
 *   p_universe = R * p_local + R * spot->p
 *
 * where R is the 3x3 rotation derived from r and spot->p is the BEAMPOS
 * offset in the local beam frame. The theta == 0 branch skips the
 * trigonometric evaluation for the common no-rotation case.
 *
 * @param[in,out] spot  Spot whose _tm[16] matrix is written.
 * @param[in]     sh    Shared beam parameters providing theta and phi [rad].
 */
static void _build_spot_tm(struct beam_spot const *spot, const struct beam_shared *sh, double tm_out[16]) {
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
    osh_vect_setup_tmatrix_bzalign_affine(spot->p, r, tm_out);
}

/**
 * @brief Run all post-parse derivations on a fully-parsed beam workspace.
 *
 * @details
 * Iterates over all spots and for each:
 *   - derives missing energy/momentum via _postparse_spot_energy()
 *   - builds the PZALIGN->UNIVERSE affine matrix via _build_spot_tm()
 *   - accumulates the cumulative weight array for SOBP spot selection
 *   - tracks the maximum energy and momentum across all spots
 *
 * The cumulative weight array is reallocated on every call so that
 * _wb_postparse() can be called again after a spot list is replaced.
 *
 * @param[in,out] wb  Workspace to finalise; must not be NULL.
 *
 * @returns OSH_OK on success, OSH_ENOMEM on allocation failure,
 *          OSH_EINVAL if any spot weight is negative.
 */
static int _wb_postparse(struct beam_workspace *wb) {
    size_t i;
    double wt_acc;
    struct particle primary_part;
    struct particle const *part;
    int rc;
    struct osh_beam_prepared *prepared;

    if (!wb) {
        return OSH_EINVAL;
    }

    part = NULL;
    if (wb->has_primary) {
        rc = _resolve_primary_particle(wb, &primary_part);
        if (rc != OSH_OK) {
            return rc;
        }
        part = &primary_part;
    }

    prepared = wb->prepared;
    if (!prepared) {
        prepared = (struct osh_beam_prepared *) calloc(1, sizeof(*prepared));
        if (!prepared) {
            return OSH_ENOMEM;
        }
        wb->prepared = prepared;
    }

    free(prepared->cum_wt);
    free(prepared->tm);
    prepared->cum_wt = NULL;
    prepared->tm = NULL;
    prepared->wt_sum = 0.0;
    prepared->emax = 0.0;
    prepared->pmax = 0.0;
    prepared->nspots = wb->nspots;

    if (wb->nspots > 0) {
        prepared->cum_wt = (double *) calloc(wb->nspots, sizeof(double));
        prepared->tm = (double *) calloc(wb->nspots * 16u, sizeof(double));
        if (!prepared->cum_wt || !prepared->tm) {
            free(prepared->cum_wt);
            free(prepared->tm);
            prepared->cum_wt = NULL;
            prepared->tm = NULL;
            prepared->nspots = 0u;
            return OSH_ENOMEM;
        }
    }

    wt_acc = 0.0;
    for (i = 0; i < wb->nspots; i++) {
        _postparse_spot_energy(&wb->spots[i], part);
        _build_spot_tm(&wb->spots[i], &wb->shared, &prepared->tm[i * 16u]);

        if (wb->spots[i].wt < 0.0) {
            return OSH_EINVAL;
        }
        wt_acc += wb->spots[i].wt;
        prepared->cum_wt[i] = wt_acc;

        if (wb->spots[i].t0 > prepared->emax) {
            prepared->emax = wb->spots[i].t0;
        }
        if (wb->spots[i].p0 > prepared->pmax) {
            prepared->pmax = wb->spots[i].p0;
        }
    }

    prepared->wt_sum = wt_acc;

    return OSH_OK;
}

/* osh_beamdef.h name arrays (osh_beam_mscat_names etc.) included via osh_beam.h */

void osh_beam_print(struct beam_workspace const *wb) {
    struct particle primary_part;
    struct osh_beam_prepared const *prepared;
    int rc;

    if (!wb) {
        return;
    }
    prepared = wb->prepared;

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
    osh_info("%-18s : %f MeV/nucleon", "DeltaE min", wb->demin);
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
    osh_info("%-18s : %.3f  MeV", "Emax", prepared ? prepared->emax : 0.0);
    osh_info("%-18s : %.3f  MeV/c", "Pmax", prepared ? prepared->pmax : 0.0);

    if (wb->has_primary) {
        rc = _resolve_primary_particle(wb, &primary_part);
        if (rc == OSH_OK) {
            osh_info("%s", "");
            osh_print_particle(&primary_part);
        }
    }

    if (wb->spots) {
        osh_beam_print_spot(&wb->spots[0]);
    }
}

void osh_beam_print_spot(struct beam_spot const *spot) {
    if (!spot) {
        return;
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
