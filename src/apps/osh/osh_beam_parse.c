#include "apps/osh/osh_beam_parse.h"

#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_beam_parse_keys.h"
#include "common/osh_diag.h"
#include "common/osh_duration.h"
#include "common/osh_file.h"
#include "common/osh_readline.h"
#include "openshieldhit/beam.h"
#include "openshieldhit/const.h"
#include "openshieldhit/status.h"
#include "particle/osh_isotope_db.h"
#include "particle/osh_particle.h"
#include "particle/osh_particle_pdg.h"

/*
 * Shared parser context passed to all beam.dat key handlers.
 *
 * Today this bundles both parser-control context (diag) and parser-owned
 * output scratch/result pointers (single-spot template + deferred USECBEAM
 * path). That keeps the handler signatures small, but the roles are mixed.
 * TODO: split this into clearer "context" vs "parse outputs" structs if the
 * parser grows further.
 */
struct beam_parse_state {
    struct osh_diag_sink const *diag;
    struct osh_beam_spot *spot_out;
    char **spotlist_path_out;
};

/*
 * Shared parse-handler signatures intentionally keep some parameters present
 * even when a specific handler does not use them. For those callback-style
 * signatures, an `unused` attribute is less noisy than sprinkling `(void)x`
 * casts through every trivial handler body.
 */
#if defined(__GNUC__) || defined(__clang__)
#define OSH_PARSE_UNUSED __attribute__((unused))
#else
#define OSH_PARSE_UNUSED
#endif

#define PARSE_HANDLER_ARGS                                                                                             \
    struct osh_beam_workspace *beam, struct oshfile *oshf, char const *args,                                           \
        struct beam_parse_state *state OSH_PARSE_UNUSED

/* ---- Parse-phase handlers ------------------------------------------------
 *
 * Each handler fills raw values from the file into the beam workspace or the
 * parser-owned template spot.
 * No derived quantities (p0 from t0, _tm matrix, emax) are computed here —
 * that is the job of the post-parse step in osh_beam_setup_from_path().
 *
 * Design: keeping parse and post-parse as separate phases makes each
 * independently testable and clearly separates "what the file says" from
 * "what the simulation needs". */

static int _resolve_input_relative_path(struct oshfile *oshf, char const *input_path, char **resolved_path_out);
static struct osh_beam_spot *spot_from_state(struct beam_parse_state *state);
static int _parse_apcorr(PARSE_HANDLER_ARGS);
static int _parse_beamdir(PARSE_HANDLER_ARGS);
static int _parse_beamdiv(PARSE_HANDLER_ARGS);
static int _parse_beampos(PARSE_HANDLER_ARGS);
static int _parse_beamsad(PARSE_HANDLER_ARGS);
static int _parse_beamsigma(PARSE_HANDLER_ARGS);
static int _parse_bmodmc(PARSE_HANDLER_ARGS);
static int _parse_bmodtrans(PARSE_HANDLER_ARGS);
static int _parse_deltae(PARSE_HANDLER_ARGS);
static int _parse_demin(PARSE_HANDLER_ARGS);
static int _parse_emtrans(PARSE_HANDLER_ARGS);
static int _parse_extspec(PARSE_HANDLER_ARGS);
static int _parse_makeln(PARSE_HANDLER_ARGS);
static int _parse_dumpevery(PARSE_HANDLER_ARGS);
static int _parse_maxtime(PARSE_HANDLER_ARGS);
static int _parse_mscat(PARSE_HANDLER_ARGS);
static int _parse_neutrfast(PARSE_HANDLER_ARGS);
static int _parse_neutrlcut(PARSE_HANDLER_ARGS);
static int _parse_nstat(PARSE_HANDLER_ARGS);
static int _parse_nucre(PARSE_HANDLER_ARGS);
static int _parse_primary(PARSE_HANDLER_ARGS);
static int _parse_rndseed(PARSE_HANDLER_ARGS);
static int _parse_stragg(PARSE_HANDLER_ARGS);
static int _parse_tmax0(PARSE_HANDLER_ARGS);
static int _parse_tcut0(PARSE_HANDLER_ARGS);
static int _parse_usebmod(PARSE_HANDLER_ARGS);
static int _parse_usecbeam(PARSE_HANDLER_ARGS);
static int _parse_useparlev(PARSE_HANDLER_ARGS);

/* ---- Dispatch table ------------------------------------------------------ */

struct _beam_dispatch_entry {
    char const *key;
    int (*handler)(PARSE_HANDLER_ARGS);
};

/* One dispatch entry per line for readable diffs; without this clang-format
 * bin-packs them into a two-column grid.  Same idiom as the key table in
 * osh_beam_parse_keys.h. */
// clang-format off
static struct _beam_dispatch_entry _dispatch_table[] = {
    {OSH_BEAM_KEY_APCORR, _parse_apcorr},
    {OSH_BEAM_KEY_BEAMDIR, _parse_beamdir},
    {OSH_BEAM_KEY_BEAMDIV, _parse_beamdiv},
    {OSH_BEAM_KEY_BEAMPOS, _parse_beampos},
    {OSH_BEAM_KEY_BEAMSAD, _parse_beamsad},
    {OSH_BEAM_KEY_BEAMSIGMA, _parse_beamsigma},
    {OSH_BEAM_KEY_BMODMC, _parse_bmodmc},
    {OSH_BEAM_KEY_BMODTRANS, _parse_bmodtrans},
    {OSH_BEAM_KEY_DELTAE, _parse_deltae},
    {OSH_BEAM_KEY_DEMIN, _parse_demin},
    {OSH_BEAM_KEY_DUMPEVERY, _parse_dumpevery},
    {OSH_BEAM_KEY_EMTRANS, _parse_emtrans},
    {OSH_BEAM_KEY_EXTSPEC, _parse_extspec},
    {OSH_BEAM_KEY_MAKELN, _parse_makeln},
    {OSH_BEAM_KEY_MAXTIME, _parse_maxtime},
    {OSH_BEAM_KEY_MSCAT, _parse_mscat},
    {OSH_BEAM_KEY_NEUTRFAST, _parse_neutrfast},
    {OSH_BEAM_KEY_NEUTRLCUT, _parse_neutrlcut},
    {OSH_BEAM_KEY_NSTAT, _parse_nstat},
    {OSH_BEAM_KEY_NUCRE, _parse_nucre},
    {OSH_BEAM_KEY_PRIMARY, _parse_primary},
    {OSH_BEAM_KEY_RNDSEED, _parse_rndseed},
    {OSH_BEAM_KEY_STRAGG, _parse_stragg},
    {OSH_BEAM_KEY_TMAX0, _parse_tmax0},
    {OSH_BEAM_KEY_TCUT0, _parse_tcut0},
    {OSH_BEAM_KEY_USEBMOD, _parse_usebmod},
    {OSH_BEAM_KEY_USECBEAM, _parse_usecbeam},
    {OSH_BEAM_KEY_USEPARLEV, _parse_useparlev},
    {NULL, NULL} /* sentinel */
};
// clang-format on

/* ---- Main parser entry point --------------------------------------------- */

enum osh_status osh_beam_parse(struct oshfile *oshf,
                               struct osh_diag_sink const *diag,
                               struct osh_beam_workspace *beam,
                               struct osh_beam_spot *spot_out,
                               char **spotlist_path_out) {
    char *lline = NULL;
    char *key = NULL;
    char *args = NULL;
    struct beam_parse_state state = {diag, spot_out, spotlist_path_out};
    int lineno;
    int i;

    if (spotlist_path_out) {
        *spotlist_path_out = NULL;
    }

    while (osh_readline_key(oshf, &lline, &key, &args, &lineno) != -1) {
        int found = 0;
        for (i = 0; key[i] != '\0'; i++) {
            key[i] = (char) tolower((unsigned char) key[i]);
        }
        for (i = 0; _dispatch_table[i].key != NULL; i++) {
            if (strcmp(_dispatch_table[i].key, key) == 0) {
                int rc = _dispatch_table[i].handler(beam, oshf, args, &state);
                if (rc != OSH_OK) {
                    free(lline);
                    return rc;
                }
                found = 1;
                break;
            }
        }
        if (!found) {
            OSH_DIAG_ERRORF(diag, "in %s line %d: unknown key '%s'", oshf->filename, lineno, key);
            free(lline);
            return OSH_EPARSE;
        }
        free(lline);
        lline = NULL;
    }

    /* Warn when USECBEAM is present but BEAMSAD was never set.
     * In that case spot x/y positions from the file are used as-is (no
     * back-projection), which is physically correct only for parallel beam
     * delivery (SAD = infinity).  This is a valid research use-case, but
     * easy to overlook — hence the warning.
     * Set BEAMSAD INF explicitly to suppress it. */
    if (spotlist_path_out && *spotlist_path_out && !beam->shared.sad_was_set) {
        OSH_DIAG_WARNF(diag,
                       "USECBEAM loaded without BEAMSAD: assuming parallel beam delivery "
                       "(SAD = infinity, no fan-out correction). Spot x/y positions are "
                       "used as-is. If your nozzle has a finite source distance, add "
                       "BEAMSAD <sadX> [sadY]. To suppress this warning for intentional "
                       "parallel delivery, set BEAMSAD INF.");
    }

    return OSH_OK;
}

static int _resolve_input_relative_path(struct oshfile *oshf, char const *input_path, char **resolved_path_out) {
    char *wdir;
    int rc;

    if (!oshf || !input_path || !resolved_path_out) {
        return OSH_EINVAL;
    }

    *resolved_path_out = NULL;
    wdir = osh_path_dirname(oshf->filename);
    /* NULL means the beam file has no directory separator (bare filename);
     * treat it as relative to the current working directory ("."). */

    rc = osh_relative_path_to_file(resolved_path_out, wdir ? wdir : ".", input_path);
    free(wdir);
    if (rc != 0) {
        return OSH_ENOMEM;
    }

    return OSH_OK;
}

static struct osh_beam_spot *spot_from_state(struct beam_parse_state *state) {
    if (!state) {
        return NULL;
    }
    return state->spot_out;
}

/* ---- Handler implementations --------------------------------------------- */

/**
 * @brief Enable aperture correction.
 *
 * APCORR is a flag card: its presence alone enables correction.
 * No argument is expected or read.
 *
 * @param[in,out] beam  Sets beam->apcorr = 1.
 * @param[in]     oshf  Unused.
 * @param[in]     args  Unused.
 *
 * @returns OSH_OK.
 */
static int _parse_apcorr(PARSE_HANDLER_ARGS) {
    (void) oshf;
    (void) args;
    beam->apcorr = 1;
    return OSH_OK;
}

/**
 * @brief Parse BEAMDIR: beam direction in spherical coordinates.
 *
 * @details
 * Syntax: BEAMDIR \<theta\> \<phi\>
 *
 * theta is the polar angle from the +Z axis [degrees, 0–180].
 * phi is the azimuthal angle in the XY plane [degrees, 0–360].
 * Both are stored internally in radians.
 *
 * @param[in,out] beam  Writes beam->shared.theta and beam->shared.phi [rad].
 * @param[in]     oshf  Used for error diagnostics.
 * @param[in]     args  Two whitespace-separated floats: theta phi.
 *
 * @returns OSH_OK on success.
 */
static int _parse_beamdir(PARSE_HANDLER_ARGS) {
    float _f[2];
    if (sscanf(args, "%f %f", &_f[0], &_f[1]) != 2) {
        OSH_DIAG_ERRORF(state->diag, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
        return OSH_EPARSE;
    }
    if (_f[0] < 0.0f || _f[0] > 180.0f) {
        OSH_DIAG_ERRORF(state->diag, "in %s line %i: theta must be within [0:180] deg", oshf->filename, oshf->lineno);
        return OSH_EPARSE;
    }
    if (_f[1] < 0.0f || _f[1] > 360.0f) {
        OSH_DIAG_ERRORF(state->diag, "in %s line %i: phi must be within [0:360] deg", oshf->filename, oshf->lineno);
        return OSH_EPARSE;
    }
    beam->shared.theta = (double) _f[0] * OSH_M_PI_180;
    beam->shared.phi = (double) _f[1] * OSH_M_PI_180;
    return OSH_OK;
}

/**
 * @brief Parse BEAMDIV: beam angular divergence and focus distance.
 *
 * @details
 * Syntax: BEAMDIV \<divX\> \<divY\> [\<focus\>]
 *
 * divX and divY are the half-angle divergences in the X and Y planes [mrad],
 * stored internally in radians. focus is the optional distance from the beam
 * source to the beam waist [cm], stored as-is; when omitted it defaults to 0.
 *
 * use_div is derived later by osh_beam_spots_set() from the final spot array.
 *
 * @param[in,out] beam  Writes shared.focus.
 * @param[in]     oshf  Used for error diagnostics.
 * @param[in]     args  Two or three whitespace-separated floats: divX divY [focus].
 *
 * @returns OSH_OK on success.
 */
static int _parse_beamdiv(PARSE_HANDLER_ARGS) {
    struct osh_beam_spot *spot;
    float _f[3] = {0.0f, 0.0f, 0.0f};
    int nread;

    spot = spot_from_state(state);
    if (!spot) {
        return OSH_ESTATE;
    }
    nread = sscanf(args, "%f %f %f", &_f[0], &_f[1], &_f[2]);
    if (nread < 2 || nread > 3) {
        OSH_DIAG_ERRORF(state->diag, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
        return OSH_EPARSE;
    }
    /* divergence is given in mrad in the file — convert to rad */
    spot->div[0] = (double) _f[0] * 0.001;
    spot->div[1] = (double) _f[1] * 0.001;
    beam->shared.focus = (double) _f[2];
    return OSH_OK;
}

/**
 * @brief Parse BEAMPOS: beam spot starting position.
 *
 * @details
 * Syntax: BEAMPOS \<X\> \<Y\> \<Z\>  [cm]
 *
 * @param[in,out] beam  Unused.
 * @param[in]     oshf  Used for error diagnostics.
 * @param[in]     args  Three whitespace-separated floats: X Y Z.
 *
 * @returns OSH_OK on success.
 */
static int _parse_beampos(PARSE_HANDLER_ARGS) {
    struct osh_beam_spot *spot;
    float _f[3];
    int i;

    (void) beam;
    spot = spot_from_state(state);
    if (!spot) {
        return OSH_ESTATE;
    }
    if (sscanf(args, "%f %f %f", &_f[0], &_f[1], &_f[2]) != 3) {
        OSH_DIAG_ERRORF(state->diag, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
        return OSH_EPARSE;
    }
    for (i = 0; i < 3; i++) {
        spot->p[i] = (double) _f[i];
    }
    return OSH_OK;
}

/**
 * @brief Parse BEAMSAD: Source-to-Axis Distance for fan-out correction.
 *
 * @details
 * Syntax: BEAMSAD \<sadX\> [\<sadY\>]  [cm]
 *
 * SAD is the positive distance from the virtual upstream point source to the
 * isocenter.  It is always a positive number regardless of beam direction,
 * because it is a physical distance along the beam axis, not a signed
 * coordinate.  It is used to tilt each primary so that it converges toward
 * the isocenter (scanning-magnet fan-out correction).
 *
 * One value gives a symmetric nozzle (same SAD for X and Y).
 * Two values allow an asymmetric nozzle (common in clinical systems).
 *
 * The special value INF (or Inf, inf) explicitly selects parallel beam
 * delivery (SAD = infinity, no fan-out correction).  This suppresses the
 * warning that would otherwise appear when USECBEAM is used without BEAMSAD.
 *
 * @param[in,out] beam  Writes shared.sad[0,1], sets shared.use_sad and
 *                      shared.sad_was_set.
 * @param[in]     oshf  Used for error diagnostics.
 * @param[in]     args  One or two whitespace-separated positive floats or INF.
 *
 * @returns OSH_OK on success, OSH_EPARSE on invalid input.
 */
static int _parse_beamsad(PARSE_HANDLER_ARGS) {
    float _f[2];
    int n = sscanf(args, "%f %f", &_f[0], &_f[1]);
    switch (n) {
    case 1:
        beam->shared.sad[0] = _f[0];
        beam->shared.sad[1] = _f[0];
        break;
    case 2:
        beam->shared.sad[0] = _f[0];
        beam->shared.sad[1] = _f[1];
        break;
    default:
        OSH_DIAG_ERRORF(state->diag, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
        return OSH_EPARSE;
    }

    beam->shared.sad_was_set = 1;

    if (isinf((float) beam->shared.sad[0]) || isinf((float) beam->shared.sad[1])) {
        /* INF means parallel beam (SAD = infinity). No fan-out correction is
         * applied and back-projection of spotlist coordinates is skipped.
         * Suppress the "USECBEAM without BEAMSAD" warning since the user
         * has explicitly declared their intent. */
        OSH_DIAG_INFOF(state->diag, "BEAMSAD INF: parallel beam delivery (no fan-out correction)");
        beam->shared.sad[0] = 0.0;
        beam->shared.sad[1] = 0.0;
        beam->shared.use_sad = 0;
    } else if (beam->shared.sad[0] > 0.0 && beam->shared.sad[1] > 0.0) {
        beam->shared.use_sad = 1;
    } else {
        OSH_DIAG_ERRORF(state->diag,
                        "in %s line %i: BEAMSAD must be a positive value or INF (got '%s')",
                        oshf->filename,
                        oshf->lineno,
                        args);
        return OSH_EPARSE;
    }
    return OSH_OK;
}

/**
 * @brief Parse BEAMSIGMA: beam spot size and shape.
 *
 * @details
 * Syntax: BEAMSIGMA \<sx\> [\<sy\>]  [cm]
 *
 * The sign of the arguments encodes the spot shape:
 *
 * If only one value is given, it is applied symmetrically: sy = sx.
 *
 *   sx > 0, sy > 0  —  Gaussian,   sigma_x = sx, sigma_y = sy
 *   sx = 0, sy = 0  —  pencil beam (point source, no lateral spread)
 *   sx >= 0, sy < 0 —  circular uniform disk, radius = |sy|  (sx ignored)
 *   sx < 0, sy < 0  —  square uniform, half-width_x = |sx|, half-width_y = |sy|
 *
 * This sign convention is inherited from SHIELD-HIT12A and avoids a
 * separate shape keyword.  The size is stored as a positive magnitude;
 * the shape enum carries the interpretation used by the sampler.
 *
 * @param[in,out] beam  Unused.
 * @param[in]     oshf  Used for error diagnostics.
 * @param[in]     args  One or two whitespace-separated floats: sx [sy].
 *
 * @returns OSH_OK on success.
 */
static int _parse_beamsigma(PARSE_HANDLER_ARGS) {
    struct osh_beam_spot *spot;
    float _f[2] = {0.0f, 0.0f};
    int nread;

    (void) beam;
    spot = spot_from_state(state);
    if (!spot) {
        return OSH_ESTATE;
    }
    nread = sscanf(args, "%f %f", &_f[0], &_f[1]);
    if (nread < 1 || nread > 2) {
        OSH_DIAG_ERRORF(state->diag, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
        return OSH_EPARSE;
    }
    if (nread == 1) {
        _f[1] = _f[0];
    }
    if (_f[0] < 0.0f && _f[1] < 0.0f) {
        spot->shape = OSH_BEAM_SHAPE_SQUARE;
        spot->size[0] = fabs(_f[0]);
        spot->size[1] = fabs(_f[1]);
    } else if (_f[0] >= 0.0f && _f[1] < 0.0f) {
        spot->shape = OSH_BEAM_SHAPE_CIRCULAR;
        spot->size[0] = fabs(_f[1]);
        spot->size[1] = 0.0;
    } else if (_f[0] > 0.0f || _f[1] > 0.0f) {
        spot->shape = OSH_BEAM_SHAPE_GAUSSIAN;
        spot->size[0] = _f[0];
        spot->size[1] = _f[1];
    } else {
        spot->shape = OSH_BEAM_SHAPE_PENCIL;
        spot->size[0] = 0.0;
        spot->size[1] = 0.0;
    }
    return OSH_OK;
}

/**
 * @brief Parse BMODMC: ripple filter Monte Carlo modifier (not yet implemented).
 *
 * @details
 * TODO: struct ripple_filter is not yet fully defined; wire this once
 * osh_beam_rifi.h exists and beam->rifi is allocated by USEBMOD.
 *
 * @param[in,out] beam  Unused.
 * @param[in]     oshf  Used for warning diagnostics.
 * @param[in]     args  Unused.
 *
 * @returns OSH_OK.
 */
static int _parse_bmodmc(PARSE_HANDLER_ARGS) {
    (void) beam;
    (void) args;
    OSH_DIAG_WARNF(
        state->diag, "in %s line %i: BMODMC parsed but rifi not yet implemented", oshf->filename, oshf->lineno);
    return OSH_OK;
}

/**
 * @brief Parse BMODTRANS: deprecated beam modifier card — warn and ignore.
 *
 * @param[in,out] beam  Unused.
 * @param[in]     oshf  Used for warning diagnostics.
 * @param[in]     args  Unused.
 *
 * @returns OSH_OK.
 */
static int _parse_bmodtrans(PARSE_HANDLER_ARGS) {
    (void) beam;
    (void) args;
    OSH_DIAG_WARNF(
        state->diag, "in %s line %i: BMODTRANS is deprecated and will be ignored", oshf->filename, oshf->lineno);
    return OSH_OK;
}

/**
 * @brief Parse DELTAE: maximum relative energy loss per transport step.
 *
 * @details
 * Syntax: DELTAE \<f\>  (dimensionless fraction, e.g. 0.005)
 *
 * Controls the step-size algorithm: a particle loses at most deltae * E
 * per step.  Smaller values improve accuracy at the cost of more steps.
 *
 * @param[in,out] beam  Writes beam->deltae.
 * @param[in]     oshf  Used for error diagnostics.
 * @param[in]     args  Single float.
 *
 * @returns OSH_OK on success.
 */
static int _parse_deltae(PARSE_HANDLER_ARGS) {
    float _f;
    if (sscanf(args, "%f", &_f) != 1) {
        OSH_DIAG_ERRORF(state->diag, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
        return OSH_EPARSE;
    }
    beam->deltae = _f;
    return OSH_OK;
}

/**
 * @brief Parse DEMIN: minimum energy step size.
 *
 * @details
 * Syntax: DEMIN \<e\>  [MeV/nucleon]
 *
 * Prevents the step-size algorithm from producing infinitely small steps
 * near the Bragg peak where dE/dx diverges.  Typical value: 0.025 MeV/nucleon.
 *
 * @param[in,out] beam  Writes beam->demin.
 * @param[in]     oshf  Used for error diagnostics.
 * @param[in]     args  Single float.
 *
 * @returns OSH_OK on success.
 */
static int _parse_demin(PARSE_HANDLER_ARGS) {
    float _f;
    if (sscanf(args, "%f", &_f) != 1) {
        OSH_DIAG_ERRORF(state->diag, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
        return OSH_EPARSE;
    }
    beam->demin = _f;
    return OSH_OK;
}

/**
 * @brief Parse EMTRANS: electromagnetic transport mode flag.
 *
 * @details
 * Syntax: EMTRANS \<mode\>
 *
 * Stored as a single character; interpretation is transport-engine specific.
 *
 * @param[in,out] beam  Writes beam->emtrans.
 * @param[in]     oshf  Used for error diagnostics.
 * @param[in]     args  Single character.
 *
 * @returns OSH_OK on success.
 */
static int _parse_emtrans(PARSE_HANDLER_ARGS) {
    if (sscanf(args, "%c", &(beam->emtrans)) != 1) {
        OSH_DIAG_ERRORF(state->diag, "in %s line %i: unknown EMTRANS mode '%s'", oshf->filename, oshf->lineno, args);
        return OSH_EPARSE;
    }
    return OSH_OK;
}

/**
 * @brief Parse EXTSPEC: external energy spectrum (not yet implemented).
 *
 * @details
 * TODO: EXTSPEC (external spectrum file) is not yet implemented.
 * Return an explicit parse error instead of exiting, so callers can
 * decide how to report configuration failures.
 *
 * @param[in,out] beam  Unused.
 * @param[in]     oshf  Used for error diagnostics.
 * @param[in]     args  Unused.
 *
 * @returns OSH_EPARSE always — EXTSPEC is not yet implemented.
 */
static int _parse_extspec(PARSE_HANDLER_ARGS) {
    (void) beam;
    (void) args;
    OSH_DIAG_ERRORF(state->diag, "in %s line %i: EXTSPEC is not yet implemented", oshf->filename, oshf->lineno);
    return OSH_EPARSE;
}

/**
 * @brief Parse PRIMARY: identify the beam's primary particle.
 *
 * @details
 * Accepts three unambiguous forms:
 *
 *   PRIMARY proton   — canonical name (case-insensitive)
 *   PRIMARY 2212     — PDG Monte Carlo numbering scheme particle code
 *   PRIMARY 6 12     — atomic number Z followed by mass number A
 *
 * The two-integer form covers all ions, including light ones such as
 * proton (Z=1 A=1), deuteron (Z=1 A=2), or helium-4 (Z=2 A=4).
 * For exotic or future particles not yet in the name table, the PDG
 * form provides a stable, registry-based fallback.
 *
 * On success, beam->primary is filled with species identity only and
 * beam->has_primary is set to 1.
 * On failure, beam->primary is left untouched and the caller must not
 * assume any particle has been set.
 *
 * @param[in,out] beam  Writes beam->primary and beam->has_primary.
 * @param[in]     oshf  Used for error diagnostics.
 * @param[in]     args  Particle name, PDG code, or two integers Z A.
 *
 * @returns OSH_OK on success, OSH_EINVAL if the particle is unrecognised
 *          or the argument format is invalid.
 */
static int _parse_primary(PARSE_HANDLER_ARGS) {
    char tok1[64];
    char tok2[64];
    char extra;
    int n;
    int pdg;
    struct particle part;

    n = sscanf(args, "%63s %63s", tok1, tok2);
    if (n < 1) {
        OSH_DIAG_WARNF(
            state->diag, "in %s line %i: PRIMARY expects a name, PDG code, or '<Z> <A>'", oshf->filename, oshf->lineno);
        return OSH_EINVAL;
    }

    if (n == 2) {
        /* Two integers: Z A */
        unsigned int z;
        unsigned int a;
        struct isotope iso;
        char ex2;
        if (sscanf(tok1, "%u%c", &z, &extra) != 1 || sscanf(tok2, "%u%c", &a, &ex2) != 1) {
            OSH_DIAG_WARNF(
                state->diag, "in %s line %i: PRIMARY with two values expects '<Z> <A>'", oshf->filename, oshf->lineno);
            return OSH_EINVAL;
        }
        if (!osh_isotope_from_za(&iso, z, a)) {
            OSH_DIAG_WARNF(
                state->diag, "in %s line %i: unknown ion with Z=%u A=%u", oshf->filename, oshf->lineno, z, a);
            return OSH_EINVAL;
        }
        pdg = OSH_PART_PDG_HIBASE + (int) z * 10000 + (int) a * 10;
        if (!osh_particle_from_pdg(&part, pdg)) {
            OSH_DIAG_WARNF(state->diag,
                           "in %s line %i: failed to construct ion with Z=%u A=%u",
                           oshf->filename,
                           oshf->lineno,
                           z,
                           a);
            return OSH_EINVAL;
        }
    } else if (sscanf(tok1, "%d%c", &pdg, &extra) == 1) {
        /* Single pure integer: PDG code */
        if (!osh_particle_from_pdg(&part, pdg)) {
            OSH_DIAG_WARNF(state->diag, "in %s line %i: unknown PDG code %d", oshf->filename, oshf->lineno, pdg);
            return OSH_EINVAL;
        }
    } else {
        /* Name lookup */
        if (!osh_particle_from_name(&part, tok1)) {
            OSH_DIAG_WARNF(state->diag, "in %s line %i: unknown particle '%s'", oshf->filename, oshf->lineno, tok1);
            return OSH_EINVAL;
        }
    }

    beam->primary.pdg = part.pdg;
    beam->primary.z = part.z;
    beam->primary.a = part.a;
    beam->has_primary = 1;
    return OSH_OK;
}

/**
 * @brief Parse MAKELN: enable source particle file output.
 *
 * @details
 * Syntax: MAKELN \<0|1\>
 *
 * 1 = write a "sour" phase-space file of all primary particles at birth;
 * 0 = skip.  Primarily used for debugging or phase-space export.
 *
 * @param[in,out] beam  Writes beam->makeln.
 * @param[in]     oshf  Used for error diagnostics.
 * @param[in]     args  Single integer: 0 or 1.
 *
 * @returns OSH_OK on success.
 */
static int _parse_makeln(PARSE_HANDLER_ARGS) {
    int _i;
    if (sscanf(args, "%i", &_i) != 1) {
        OSH_DIAG_ERRORF(state->diag, "in %s line %i: unknown MAKELN mode '%s'", oshf->filename, oshf->lineno, args);
        return OSH_EPARSE;
    }
    beam->makeln = (char) _i;
    return OSH_OK;
}

/**
 * @brief Parse MAXTIME: wall-time budget for a graceful, partial-result stop.
 *
 * @details
 * Syntax: MAXTIME \<duration\>
 *
 * The argument is a duration token accepted by @ref osh_parse_duration: a
 * non-negative number with an optional unit suffix (`s`, `m`, `h`; bare number
 * = seconds).  When the budget elapses the run stops cleanly at the next safe
 * point — in-flight histories finish, all secondary families drain, and the
 * partial result is saved normalised by the true completed-primary count.
 * 0 means unlimited.  A CLI `--max-time` overrides this card.
 *
 * @param[in,out] beam  Writes beam->wall_budget_s [s].
 * @param[in]     oshf  Used for error diagnostics.
 * @param[in]     args  A single duration token, e.g. "30m" or "3600".
 *
 * @returns OSH_OK on success, OSH_EPARSE on a malformed duration.
 */
static int _parse_maxtime(PARSE_HANDLER_ARGS) {
    char tok[64];
    double seconds;

    if (sscanf(args, "%63s", tok) != 1 || !osh_parse_duration(tok, &seconds)) {
        OSH_DIAG_ERRORF(state->diag,
                        "in %s line %i: invalid MAXTIME duration '%s' (use e.g. 30s, 30m, 1h, or 3600)",
                        oshf->filename,
                        oshf->lineno,
                        args);
        return OSH_EPARSE;
    }
    beam->wall_budget_s = seconds;
    return OSH_OK;
}

/**
 * @brief Parse DUMPEVERY: wall-time cadence for periodic partial-result dumps.
 *
 * @details
 * Syntax: DUMPEVERY \<duration\>
 *
 * The argument is a duration token accepted by @ref osh_parse_duration (a
 * non-negative number with an optional `s`/`m`/`h` suffix; bare number = seconds).
 * Every such interval of wall time the run overwrites its output files with the
 * current partial result, taken at a family-complete checkpoint so the snapshot is
 * physically exact.  0 disables periodic dumps.  A CLI `--dump-every` overrides
 * this card; the count-cadence equivalent is the NSTAT save step (@ref nsave) /
 * `--dump-every-primaries`.
 *
 * @param[in,out] beam  Writes beam->dump_every_s [s].
 * @param[in]     oshf  Used for error diagnostics.
 * @param[in]     args  A single duration token, e.g. "10m" or "600".
 *
 * @returns OSH_OK on success, OSH_EPARSE on a malformed duration.
 */
static int _parse_dumpevery(PARSE_HANDLER_ARGS) {
    char tok[64];
    double seconds;

    if (sscanf(args, "%63s", tok) != 1 || !osh_parse_duration(tok, &seconds)) {
        OSH_DIAG_ERRORF(state->diag,
                        "in %s line %i: invalid DUMPEVERY duration '%s' (use e.g. 30s, 10m, 1h, or 600)",
                        oshf->filename,
                        oshf->lineno,
                        args);
        return OSH_EPARSE;
    }
    beam->dump_every_s = seconds;
    return OSH_OK;
}

/**
 * @brief Parse MSCAT: multiple Coulomb scattering model.
 *
 * @details
 * Syntax: MSCAT \<mode\>
 *
 *   0 — off
 *   1 — Gaussian (Rossi-Greisen, fast)
 *   2 — Molière  (more accurate, standard choice for clinical use)
 *   3 — Wentzel-VI/Urban — reserved: accepted here but not yet implemented;
 *       transport setup rejects it with OSH_ENOTSUP (cf. Vavilov straggling).
 *
 * @param[in,out] beam  Writes beam->scatter.
 * @param[in]     oshf  Used for error diagnostics.
 * @param[in]     args  Single integer mode.
 *
 * @returns OSH_OK on success.
 */
static int _parse_mscat(PARSE_HANDLER_ARGS) {
    int _i;
    if (sscanf(args, "%i", &_i) != 1) {
        OSH_DIAG_ERRORF(state->diag, "in %s line %i: unknown MSCAT mode '%s'", oshf->filename, oshf->lineno, args);
        return OSH_EPARSE;
    }
    beam->scatter = (char) _i;
    if (beam->scatter > OSH_BEAM_MSCAT_WENTZEL || beam->scatter < OSH_BEAM_MSCAT_OFF) {
        OSH_DIAG_ERRORF(
            state->diag, "in %s line %i: invalid MSCAT mode '%i'", oshf->filename, oshf->lineno, beam->scatter);
        return OSH_EPARSE;
    }
    return OSH_OK;
}

/**
 * @brief Parse NEUTRFAST: fast neutron transport mode.
 *
 * @details
 * Syntax: NEUTRFAST \<mode\>
 *
 * Controls whether neutrons above a threshold energy use a simplified
 * (faster) transport algorithm instead of full analog tracking.
 *
 * @param[in,out] beam  Writes beam->neutrfast.
 * @param[in]     oshf  Used for error diagnostics.
 * @param[in]     args  Single integer mode.
 *
 * @returns OSH_OK on success.
 */
static int _parse_neutrfast(PARSE_HANDLER_ARGS) {
    int _i;
    if (sscanf(args, "%i", &_i) != 1) {
        OSH_DIAG_ERRORF(state->diag, "in %s line %i: unknown NEUTRFAST mode '%s'", oshf->filename, oshf->lineno, args);
        return OSH_EPARSE;
    }
    beam->neutrfast = (char) _i;
    return OSH_OK;
}

/**
 * @brief Parse NEUTRLCUT: lower energy cutoff for neutron transport.
 *
 * @details
 * Syntax: NEUTRLCUT \<e\>  [MeV]
 *
 * Neutrons at or below this threshold are killed; in positive-density material
 * their remaining kinetic energy is point-deposited.  0.0 keeps the transport
 * default.
 *
 * @param[in,out] beam  Writes beam->ncut.
 * @param[in]     oshf  Used for error diagnostics.
 * @param[in]     args  Single float.
 *
 * @returns OSH_OK on success.
 */
static int _parse_neutrlcut(PARSE_HANDLER_ARGS) {
    float _f;
    if (sscanf(args, "%f", &_f) != 1) {
        OSH_DIAG_ERRORF(state->diag, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
        return OSH_EPARSE;
    }
    beam->ncut = _f;
    return OSH_OK;
}

/**
 * @brief Parse NSTAT: number of primary histories and save interval.
 *
 * @details
 * Syntax: NSTAT \<n\> [\<step\>]
 *
 * n    — total number of primary particle histories to simulate.
 * step — interval at which intermediate results are written to disk;
 *        -1 disables intermediate saves (write only at the end).
 *
 * @param[in,out] beam  Writes beam->nstat and beam->nsave.
 * @param[in]     oshf  Used for error diagnostics.
 * @param[in]     args  One or two integers: n [step].
 *
 * @returns OSH_OK on success.
 */
static int _parse_nstat(PARSE_HANDLER_ARGS) {
    int _i[2] = {0, 0};
    int nread;

    nread = sscanf(args, "%i %i", &_i[0], &_i[1]);
    if (nread < 1 || nread > 2) {
        OSH_DIAG_ERRORF(state->diag, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
        return OSH_EPARSE;
    }
    if (_i[0] < 0) {
        OSH_DIAG_ERRORF(state->diag, "in %s line %i: NSTAT must be >= 0", oshf->filename, oshf->lineno);
        return OSH_EPARSE;
    }

    beam->nstat = (size_t) _i[0];
    beam->nsave = (_i[1] > 0) ? (size_t) _i[1] : 0u;
    return OSH_OK;
}

/**
 * @brief Parse NUCRE: nuclear reaction switch.
 *
 * @details
 * Syntax: NUCRE \<0|1|2|3\>
 *
 *   0 — off (electromagnetic transport only, no hadronic interactions)
 *   1 — all nuclear reactions (Tripathi absorption + pp elastic scattering)
 *   2 — elastic only (pp elastic scattering, no inelastic absorption)
 *   3 — inelastic only (Tripathi absorption; diagnostic mode)
 *
 * Disabling nuclear reactions is useful for pure range/dose validation
 * where hadronic secondaries would complicate comparison with analytic models.
 *
 * @param[in,out] beam  Writes beam->nuclear_inelastic and beam->nuclear_elastic.
 * @param[in]     oshf  Used for error diagnostics.
 * @param[in]     args  Single integer: 0, 1, 2, or 3.
 *
 * @returns OSH_OK on success.
 */
static int _parse_nucre(PARSE_HANDLER_ARGS) {
    int _i;
    if (sscanf(args, "%i", &_i) != 1) {
        OSH_DIAG_ERRORF(state->diag, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
        return OSH_EPARSE;
    }
    if (_i < 0 || _i > 3) {
        OSH_DIAG_ERRORF(state->diag,
                        "in %s line %i: invalid NUCRE mode '%i' (expected 0, 1, 2, or 3)",
                        oshf->filename,
                        oshf->lineno,
                        _i);
        return OSH_EPARSE;
    }
    beam->nuclear_inelastic = (_i == 1 || _i == 3) ? 1 : 0;
    beam->nuclear_elastic = (_i == 1 || _i == 2) ? 1 : 0;
    return OSH_OK;
}

/**
 * @brief Parse RNDSEED: random number generator seed.
 *
 * @details
 * Syntax: RNDSEED \<seed\>
 *
 * Fixing the seed makes a run bit-for-bit reproducible.  Using different
 * seeds produces statistically independent runs, which is the standard
 * approach for Monte Carlo uncertainty estimation.
 *
 * @param[in,out] beam  Writes beam->rndseed.
 * @param[in]     oshf  Used for error diagnostics.
 * @param[in]     args  Single integer.
 *
 * @returns OSH_OK on success.
 */
static int _parse_rndseed(PARSE_HANDLER_ARGS) {
    int _i;
    if (sscanf(args, "%i", &_i) != 1) {
        OSH_DIAG_ERRORF(state->diag, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
        return OSH_EPARSE;
    }
    beam->rndseed = _i;
    return OSH_OK;
}

/**
 * @brief Parse STRAGG: energy straggling model.
 *
 * @details
 * Syntax: STRAGG \<mode\>
 *
 *   0 — off
 *   1 — Gaussian (Bohr straggling, fast)
 *   2 — Vavilov  (statistically correct for thin absorbers)
 *   3 — Urbán    (reserved; not yet implemented — rejected at transport setup)
 *
 * Vavilov converges to Gaussian for thick absorbers, so mode 2 is the
 * safe default.  Gaussian may be preferred when speed is critical and
 * the absorber is thick compared to the range.
 *
 * @param[in,out] beam  Writes beam->straggl.
 * @param[in]     oshf  Used for error diagnostics.
 * @param[in]     args  Single integer mode.
 *
 * @returns OSH_OK on success.
 */
static int _parse_stragg(PARSE_HANDLER_ARGS) {
    int _i;
    if (sscanf(args, "%i", &_i) != 1) {
        OSH_DIAG_ERRORF(state->diag, "in %s line %i: unknown STRAGG mode '%s'", oshf->filename, oshf->lineno, args);
        return OSH_EPARSE;
    }
    beam->straggl = (char) _i;
    if (beam->straggl > OSH_BEAM_STRAGG_URBAN || beam->straggl < OSH_BEAM_STRAGG_OFF) {
        OSH_DIAG_ERRORF(
            state->diag, "in %s line %i: invalid STRAGG mode '%i'", oshf->filename, oshf->lineno, beam->straggl);
        return OSH_EPARSE;
    }
    return OSH_OK;
}

/**
 * @brief Parse TMAX0: initial kinetic energy or momentum of the primary particle.
 *
 * @details
 * Syntax: TMAX0 \<value\> [\<spread\>]
 *
 * Sign convention (inherited from SHIELD-HIT12A):
 *   positive value — kinetic energy T0 [MeV/nucleon] for ions, [MeV] otherwise
 *   negative value — total momentum p0 [MeV/c] (stored as |value|)
 *
 * The same sign convention applies to the optional spread parameter.
 * Only one of t0/p0 (and tsigma/psigma) is set here — whichever the file
 * specifies; the other is left at zero.  The post-parse step derives the
 * complementary quantity using the particle mass, and checks which was given
 * by testing which field is non-zero.
 *
 * @param[in,out] beam  Unused.
 * @param[in]     oshf  Used for error diagnostics.
 * @param[in]     args  One or two floats: value [spread].
 *
 * @returns OSH_OK on success.
 */
static int _parse_tmax0(PARSE_HANDLER_ARGS) {
    struct osh_beam_spot *spot;
    float _f[2];

    (void) beam;
    spot = spot_from_state(state);
    if (!spot) {
        return OSH_ESTATE;
    }
    _f[0] = 0.0f;
    _f[1] = 0.0f;
    if (sscanf(args, "%f %f", &_f[0], &_f[1]) < 1) {
        OSH_DIAG_ERRORF(state->diag, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
        return OSH_EPARSE;
    }

    if (_f[0] < 0.0f) {
        spot->p0 = fabs((double) _f[0]); /* momentum given */
        spot->t0 = 0.0;
        spot->t0_per_nucleon = 0;
    } else {
        spot->t0 = (double) _f[0]; /* energy per nucleon until post-parse for ions */
        spot->p0 = 0.0;
        spot->t0_per_nucleon = 1;
        if (_f[0] < OSH_BEAM_TMIN) {
            OSH_DIAG_ERRORF(state->diag,
                            "in %s line %i: TMAX0 %.4f MeV/nucleon is below the minimum beam energy %.4f MeV/nucleon",
                            oshf->filename,
                            oshf->lineno,
                            (double) _f[0],
                            OSH_BEAM_TMIN);
            return OSH_EPARSE;
        }
    }

    if (_f[1] < 0.0f) {
        spot->psigma = fabs((double) _f[1]); /* momentum spread given */
        spot->tsigma = 0.0;
        spot->tsigma_per_nucleon = 0;
    } else {
        spot->tsigma = (double) _f[1]; /* energy spread per nucleon until post-parse for ions */
        spot->psigma = 0.0;
        spot->tsigma_per_nucleon = 1;
    }

    return OSH_OK;
}

/**
 * @brief Parse TCUT0: primary energy sampling window (truncated Gaussian).
 *
 * @details
 * Syntax: TCUT0 \<lower\> \<upper\>  [MeV/nucleon]
 *
 * Confines the TMAX0 energy spread to [lower, upper]: primary energies are
 * drawn from N(t0, tsigma^2) **truncated** to the window rather than from the
 * plain untruncated Gaussian used when TCUT0 is absent.  The draw is an exact
 * inverse-CDF transform (see osh_rng_gauss_trunc.h), so the window may sit
 * anywhere relative to t0.
 *
 * This bounds only the energy primaries are born with.  It is unrelated to the
 * ion transport cutoff (beam->tcut_transport), which kills particles that slow
 * down past a threshold in flight and is not settable from beam.dat.
 *
 * Both bounds are converted from MeV/nucleon to absolute MeV for the primary
 * species during the post-parse step (see _wb_postparse() in osh_beam.c),
 * mirroring how TMAX0's t0/tsigma are scaled.
 *
 * Rejected input: a missing argument, a negative bound (kinetic energy cannot
 * be negative), an upper bound below the lower one, and a non-positive upper
 * bound — the last would force every primary to 0 MeV, which is never what the
 * user meant and would otherwise be indistinguishable from "no TCUT0".
 *
 * @param[in,out] beam  Writes beam->tcut_lower and beam->tcut_upper.
 * @param[in]     oshf  Used for error diagnostics.
 * @param[in]     args  Two floats: lower upper.
 *
 * @returns OSH_OK on success, OSH_EPARSE on any rejected input above.
 */
static int _parse_tcut0(PARSE_HANDLER_ARGS) {
    float _f[2];
    _f[0] = 0.0f;
    _f[1] = 0.0f;
    if (sscanf(args, "%f %f", &_f[0], &_f[1]) < 2) {
        OSH_DIAG_ERRORF(state->diag, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
        return OSH_EPARSE;
    }
    if (_f[0] < 0.0f || _f[1] < 0.0f) {
        OSH_DIAG_ERRORF(state->diag, "in %s line %i: TCUT0 bounds must not be negative", oshf->filename, oshf->lineno);
        return OSH_EPARSE;
    }
    if (_f[1] <= 0.0f) {
        OSH_DIAG_ERRORF(state->diag, "in %s line %i: TCUT0 upper bound must be > 0", oshf->filename, oshf->lineno);
        return OSH_EPARSE;
    }
    if (_f[0] > _f[1]) {
        OSH_DIAG_ERRORF(
            state->diag, "in %s line %i: TCUT0 upper bound must be >= lower bound", oshf->filename, oshf->lineno);
        return OSH_EPARSE;
    }
    beam->tcut_lower = _f[0];
    beam->tcut_upper = _f[1];
    return OSH_OK;
}

/**
 * @brief Parse USEBMOD: load a ripple filter beam modifier (not yet implemented).
 *
 * @details
 * Syntax: USEBMOD \<scale\> \<filename\>
 *
 * TODO: osh_beam_rifi_load is not yet declared; wire once osh_beam_rifi.h
 * exists and beam->rifi is allocated.
 *
 * @param[in,out] beam  Unused (rifi not yet allocated).
 * @param[in]     oshf  Used for warning diagnostics.
 * @param[in]     args  Float scale factor and path to ripple filter file.
 *
 * @returns OSH_OK.
 */
static int _parse_usebmod(PARSE_HANDLER_ARGS) {
    float _f;
    char tmpstr[256];
    char *_path = NULL;
    (void) beam;
    (void) state;
    if (sscanf(args, "%f %255s", &_f, tmpstr) > 2) {
        OSH_DIAG_ERRORF(state->diag, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
        return OSH_EPARSE;
    }
    if (_resolve_input_relative_path(oshf, tmpstr, &_path) != OSH_OK) {
        return OSH_ENOMEM;
    }
    OSH_DIAG_WARNF(
        state->diag, "in %s line %i: USEBMOD parsed but rifi loader not yet implemented", oshf->filename, oshf->lineno);
    free(_path);
    return OSH_OK;
}

/**
 * @brief Parse USECBEAM: load an external spot list for SOBP delivery.
 *
 * @details
 * Syntax: USECBEAM \<filename\>
 *
 * The spot list file defines one beam spot per line (energy, position,
 * size, weight). When present, app-layer import later replaces the
 * single parsed template spot, which then serves as the fallback default
 * source for columns not specified in the spot file.
 *
 * The path is resolved relative to the directory containing beam.dat,
 * so clinical plans can reference files by relative path without
 * hard-coding absolute directories.
 *
 * @param[in,out] beam  Sets OSH_BEAM_MODE_SOBP.
 * @param[in]     oshf  Used for error diagnostics.
 * @param[in]     args  Path to the spot list file (relative or absolute).
 *
 * @returns OSH_OK on success.
 */
static int _parse_usecbeam(PARSE_HANDLER_ARGS) {
    char tmpstr[256];
    char *_path = NULL;
    if (sscanf(args, "%255s", tmpstr) != 1) {
        OSH_DIAG_ERRORF(state->diag, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
        return OSH_EPARSE;
    }
    if (_resolve_input_relative_path(oshf, tmpstr, &_path) != OSH_OK) {
        return OSH_ENOMEM;
    }
    if (state && state->spotlist_path_out) {
        free(*state->spotlist_path_out);
        *state->spotlist_path_out = _path;
        OSH_DIAG_INFOF(state->diag, "USECBEAM enabled: queued external spotlist %s", *state->spotlist_path_out);
        _path = NULL;
    }
    free(_path);
    beam->beam_mode = OSH_BEAM_MODE_SOBP;
    return OSH_OK;
}

/**
 * @brief Parse USEPARLEV: load a parallel lever beam optics file (not yet implemented).
 *
 * @details
 * Syntax: USEPARLEV \<filename\>
 *
 * TODO: osh_beam_parlev_load is not yet declared; wire once the parlev
 * loader header exists.
 *
 * @param[in,out] beam  Unused (parlev not yet loaded).
 * @param[in]     oshf  Used for warning diagnostics.
 * @param[in]     args  Path to the parlev file.
 *
 * @returns OSH_OK.
 */
static int _parse_useparlev(PARSE_HANDLER_ARGS) {
    char tmpstr[256];
    char *_path = NULL;
    (void) beam;
    (void) state;
    if (sscanf(args, "%255s", tmpstr) != 1) {
        OSH_DIAG_ERRORF(state->diag, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
        return OSH_EPARSE;
    }
    if (_resolve_input_relative_path(oshf, tmpstr, &_path) != OSH_OK) {
        return OSH_ENOMEM;
    }
    OSH_DIAG_WARNF(state->diag,
                   "in %s line %i: USEPARLEV parsed but parlev loader not yet implemented",
                   oshf->filename,
                   oshf->lineno);
    free(_path);
    return OSH_OK;
}
