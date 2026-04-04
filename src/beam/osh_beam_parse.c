#include "beam/osh_beam_parse.h"

#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "beam/osh_beam.h"
#include "beam/osh_beam_parse_keys.h"
#include "common/osh_const.h"
#include "common/osh_exit.h"
#include "common/osh_file.h"
#include "common/osh_logger.h"
#include "common/osh_rc.h"
#include "common/osh_readline.h"

/* ---- Parse-phase handlers ------------------------------------------------
 *
 * Each handler fills raw values from the file into the beam workspace.
 * No derived quantities (p0 from t0, _tm matrix, emax) are computed here —
 * that is the job of the post-parse step in osh_beam_setup_from_path().
 *
 * Design: keeping parse and post-parse as separate phases makes each
 * independently testable and clearly separates "what the file says" from
 * "what the simulation needs". */

static int _parse_apcorr(struct beam_workspace *beam, struct oshfile *oshf, char const *args);
static int _parse_beamdir(struct beam_workspace *beam, struct oshfile *oshf, char const *args);
static int _parse_beamdiv(struct beam_workspace *beam, struct oshfile *oshf, char const *args);
static int _parse_beampos(struct beam_workspace *beam, struct oshfile *oshf, char const *args);
static int _parse_beamsad(struct beam_workspace *beam, struct oshfile *oshf, char const *args);
static int _parse_beamsigma(struct beam_workspace *beam, struct oshfile *oshf, char const *args);
static int _parse_bmodmc(struct beam_workspace *beam, struct oshfile *oshf, char const *args);
static int _parse_bmodtrans(struct beam_workspace *beam, struct oshfile *oshf, char const *args);
static int _parse_deltae(struct beam_workspace *beam, struct oshfile *oshf, char const *args);
static int _parse_demin(struct beam_workspace *beam, struct oshfile *oshf, char const *args);
static int _parse_emtrans(struct beam_workspace *beam, struct oshfile *oshf, char const *args);
static int _parse_extspec(struct beam_workspace *beam, struct oshfile *oshf, char const *args);
static int _parse_hiproj(struct beam_workspace *beam, struct oshfile *oshf, char const *args);
static int _parse_jpart0(struct beam_workspace *beam, struct oshfile *oshf, char const *args);
static int _parse_makeln(struct beam_workspace *beam, struct oshfile *oshf, char const *args);
static int _parse_mscat(struct beam_workspace *beam, struct oshfile *oshf, char const *args);
static int _parse_neutrfast(struct beam_workspace *beam, struct oshfile *oshf, char const *args);
static int _parse_neutrlcut(struct beam_workspace *beam, struct oshfile *oshf, char const *args);
static int _parse_nstat(struct beam_workspace *beam, struct oshfile *oshf, char const *args);
static int _parse_nucre(struct beam_workspace *beam, struct oshfile *oshf, char const *args);
static int _parse_rndseed(struct beam_workspace *beam, struct oshfile *oshf, char const *args);
static int _parse_stragg(struct beam_workspace *beam, struct oshfile *oshf, char const *args);
static int _parse_tmax0(struct beam_workspace *beam, struct oshfile *oshf, char const *args);
static int _parse_tcut0(struct beam_workspace *beam, struct oshfile *oshf, char const *args);
static int _parse_usebmod(struct beam_workspace *beam, struct oshfile *oshf, char const *args);
static int _parse_usecbeam(struct beam_workspace *beam, struct oshfile *oshf, char const *args);
static int _parse_useparlev(struct beam_workspace *beam, struct oshfile *oshf, char const *args);

/* ---- Dispatch table ------------------------------------------------------ */

struct _beam_dispatch_entry {
    char const *key;
    int (*handler)(struct beam_workspace *beam, struct oshfile *oshf, char const *args);
};

static struct _beam_dispatch_entry _dispatch_table[] = {
    {OSH_BEAM_KEY_APCORR, _parse_apcorr},       {OSH_BEAM_KEY_BEAMDIR, _parse_beamdir},
    {OSH_BEAM_KEY_BEAMDIV, _parse_beamdiv},     {OSH_BEAM_KEY_BEAMPOS, _parse_beampos},
    {OSH_BEAM_KEY_BEAMSAD, _parse_beamsad},     {OSH_BEAM_KEY_BEAMSIGMA, _parse_beamsigma},
    {OSH_BEAM_KEY_BMODMC, _parse_bmodmc},       {OSH_BEAM_KEY_BMODTRANS, _parse_bmodtrans},
    {OSH_BEAM_KEY_DELTAE, _parse_deltae},       {OSH_BEAM_KEY_DEMIN, _parse_demin},
    {OSH_BEAM_KEY_EMTRANS, _parse_emtrans},     {OSH_BEAM_KEY_EXTSPEC, _parse_extspec},
    {OSH_BEAM_KEY_HIPROJ, _parse_hiproj},       {OSH_BEAM_KEY_JPART0, _parse_jpart0},
    {OSH_BEAM_KEY_MAKELN, _parse_makeln},       {OSH_BEAM_KEY_MSCAT, _parse_mscat},
    {OSH_BEAM_KEY_NEUTRFAST, _parse_neutrfast}, {OSH_BEAM_KEY_NEUTRLCUT, _parse_neutrlcut},
    {OSH_BEAM_KEY_NSTAT, _parse_nstat},         {OSH_BEAM_KEY_NUCRE, _parse_nucre},
    {OSH_BEAM_KEY_RNDSEED, _parse_rndseed},     {OSH_BEAM_KEY_STRAGG, _parse_stragg},
    {OSH_BEAM_KEY_TMAX0, _parse_tmax0},         {OSH_BEAM_KEY_TCUT0, _parse_tcut0},
    {OSH_BEAM_KEY_USEBMOD, _parse_usebmod},     {OSH_BEAM_KEY_USECBEAM, _parse_usecbeam},
    {OSH_BEAM_KEY_USEPARLEV, _parse_useparlev}, {NULL, NULL} /* sentinel */
};

/* ---- Main parser entry point --------------------------------------------- */

int osh_beam_parse(struct oshfile *oshf, struct beam_workspace *beam) {
    char *lline = NULL;
    char *key = NULL;
    char *args = NULL;
    int lineno;
    int i;

    while (osh_readline_key(oshf, &lline, &key, &args, &lineno) != -1) {
        int found = 0;
        for (i = 0; key[i] != '\0'; i++) {
            key[i] = (char) tolower((unsigned char) key[i]);
        }
        for (i = 0; _dispatch_table[i].key != NULL; i++) {
            if (strcmp(_dispatch_table[i].key, key) == 0) {
                _dispatch_table[i].handler(beam, oshf, args);
                found = 1;
                break;
            }
        }
        if (!found) {
            osh_warn("Line %d: Unknown key '%s'", lineno, key);
        }
        free(lline);
        lline = NULL;
    }
    return OSH_OK;
}

/* ---- Handler implementations --------------------------------------------- */

static int _parse_apcorr(struct beam_workspace *beam, struct oshfile *oshf, char const *args) {
    (void) oshf;
    (void) args;
    beam->apcorr = 1;
    return OSH_OK;
}

static int _parse_beamdir(struct beam_workspace *beam, struct oshfile *oshf, char const *args) {
    float _f[2];
    if (sscanf(args, "%f %f", &_f[0], &_f[1]) != 2) {
        osh_error(EX_CONFIG, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
    }
    if (_f[0] < 0.0f || _f[0] > 180.0f) {
        osh_error(EX_CONFIG, "in %s line %i: theta must be within [0:180] deg", oshf->filename, oshf->lineno);
    }
    if (_f[1] < 0.0f || _f[1] > 360.0f) {
        osh_error(EX_CONFIG, "in %s line %i: phi must be within [0:360] deg", oshf->filename, oshf->lineno);
    }
    beam->shared.theta = (double) _f[0] * OSH_M_PI_180;
    beam->shared.phi = (double) _f[1] * OSH_M_PI_180;
    return OSH_OK;
}

static int _parse_beamdiv(struct beam_workspace *beam, struct oshfile *oshf, char const *args) {
    float _f[3];
    if (sscanf(args, "%f %f %f", &_f[0], &_f[1], &_f[2]) > 3) {
        osh_error(EX_CONFIG, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
    }
    /* divergence is given in mrad in the file — convert to rad */
    beam->spots[0].div[0] = (double) _f[0] * 0.001;
    beam->spots[0].div[1] = (double) _f[1] * 0.001;
    beam->shared.focus = (double) _f[2];
    if (fabs(_f[0]) > 0.0 || fabs(_f[1]) > 0.0) {
        beam->shared.use_div = 1;
    }
    return OSH_OK;
}

static int _parse_beampos(struct beam_workspace *beam, struct oshfile *oshf, char const *args) {
    float _f[3];
    int i;
    if (sscanf(args, "%f %f %f", &_f[0], &_f[1], &_f[2]) != 3) {
        osh_error(EX_CONFIG, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
    }
    for (i = 0; i < 3; i++) {
        beam->spots[0].p[i] = (double) _f[i];
    }
    return OSH_OK;
}

static int _parse_beamsad(struct beam_workspace *beam, struct oshfile *oshf, char const *args) {
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
        osh_error(EX_CONFIG, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
    }
    if (beam->shared.sad[0] > 0.0 && beam->shared.sad[1] > 0.0) {
        beam->shared.use_sad = 1;
    } else {
        osh_error(EX_CONFIG, "in %s line %i: SAD must be > 0.0", oshf->filename, oshf->lineno);
    }
    return OSH_OK;
}

static int _parse_beamsigma(struct beam_workspace *beam, struct oshfile *oshf, char const *args) {
    float _f[2];
    if (sscanf(args, "%f %f", &_f[0], &_f[1]) > 2) {
        osh_error(EX_CONFIG, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
    }
    if (_f[0] < 0.0f && _f[1] < 0.0f) {
        beam->spots[0].shape = OSH_BEAM_SHAPE_SQUARE;
        beam->spots[0].size[0] = fabs(_f[0]);
        beam->spots[0].size[1] = fabs(_f[1]);
    } else if (_f[0] >= 0.0f && _f[1] < 0.0f) {
        beam->spots[0].shape = OSH_BEAM_SHAPE_CIRCULAR;
        beam->spots[0].size[0] = fabs(_f[1]);
        beam->spots[0].size[1] = 0.0;
    } else if (_f[0] > 0.0f || _f[1] > 0.0f) {
        beam->spots[0].shape = OSH_BEAM_SHAPE_GAUSSIAN;
        beam->spots[0].size[0] = _f[0];
        beam->spots[0].size[1] = _f[1];
    } else {
        beam->spots[0].shape = OSH_BEAM_SHAPE_PENCIL;
        beam->spots[0].size[0] = 0.0;
        beam->spots[0].size[1] = 0.0;
    }
    return OSH_OK;
}

static int _parse_bmodmc(struct beam_workspace *beam, struct oshfile *oshf, char const *args) {
    /* TODO: struct ripple_filter is not yet fully defined; wire this once
     * osh_beam_rifi.h exists and beam->rifi is allocated by USEBMOD. */
    (void) beam;
    (void) args;
    osh_warn("in %s line %i: BMODMC parsed but rifi not yet implemented", oshf->filename, oshf->lineno);
    return OSH_OK;
}

static int _parse_bmodtrans(struct beam_workspace *beam, struct oshfile *oshf, char const *args) {
    /* BMODTRANS is deprecated — warn and ignore. */
    (void) beam;
    (void) args;
    osh_warn("in %s line %i: BMODTRANS is deprecated and will be ignored", oshf->filename, oshf->lineno);
    return OSH_OK;
}

static int _parse_deltae(struct beam_workspace *beam, struct oshfile *oshf, char const *args) {
    float _f;
    if (sscanf(args, "%f", &_f) != 1) {
        osh_error(EX_CONFIG, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
    }
    beam->deltae = _f;
    return OSH_OK;
}

static int _parse_demin(struct beam_workspace *beam, struct oshfile *oshf, char const *args) {
    float _f;
    if (sscanf(args, "%f", &_f) != 1) {
        osh_error(EX_CONFIG, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
    }
    beam->demin = _f;
    return OSH_OK;
}

static int _parse_emtrans(struct beam_workspace *beam, struct oshfile *oshf, char const *args) {
    if (sscanf(args, "%c", &(beam->emtrans)) != 1) {
        osh_error(EX_CONFIG, "in %s line %i: unknown EMTRANS mode '%s'", oshf->filename, oshf->lineno, args);
    }
    return OSH_OK;
}

static int _parse_extspec(struct beam_workspace *beam, struct oshfile *oshf, char const *args) {
    /* TODO: EXTSPEC (external spectrum) is not yet implemented. */
    (void) beam;
    (void) args;
    osh_error(EX_CONFIG, "in %s line %i: EXTSPEC is not yet implemented", oshf->filename, oshf->lineno);
    return OSH_OK;
}

static int _parse_hiproj(struct beam_workspace *beam, struct oshfile *oshf, char const *args) {
    /* TODO: particle lookup by A/Z not yet wired — requires particle table. */
    (void) beam;
    (void) args;
    osh_warn("in %s line %i: HIPROJ parsed but particle lookup not yet implemented", oshf->filename, oshf->lineno);
    return OSH_OK;
}

static int _parse_jpart0(struct beam_workspace *beam, struct oshfile *oshf, char const *args) {
    /* TODO: particle lookup by JPART0 index not yet wired — requires particle table. */
    (void) beam;
    (void) args;
    osh_warn("in %s line %i: JPART0 parsed but particle lookup not yet implemented", oshf->filename, oshf->lineno);
    return OSH_OK;
}

static int _parse_makeln(struct beam_workspace *beam, struct oshfile *oshf, char const *args) {
    int _i;
    if (sscanf(args, "%i", &_i) != 1) {
        osh_error(EX_CONFIG, "in %s line %i: unknown MAKELN mode '%s'", oshf->filename, oshf->lineno, args);
    }
    beam->makeln = (char) _i;
    return OSH_OK;
}

static int _parse_mscat(struct beam_workspace *beam, struct oshfile *oshf, char const *args) {
    int _i;
    if (sscanf(args, "%i", &_i) != 1) {
        osh_error(EX_CONFIG, "in %s line %i: unknown MSCAT mode '%s'", oshf->filename, oshf->lineno, args);
    }
    beam->scatter = (char) _i;
    if (beam->scatter > OSH_BEAM_MSCAT_MOLIERE || beam->scatter < OSH_BEAM_MSCAT_OFF) {
        osh_error(EX_CONFIG, "in %s line %i: invalid MSCAT mode '%i'", oshf->filename, oshf->lineno, beam->scatter);
    }
    return OSH_OK;
}

static int _parse_neutrfast(struct beam_workspace *beam, struct oshfile *oshf, char const *args) {
    int _i;
    if (sscanf(args, "%i", &_i) != 1) {
        osh_error(EX_CONFIG, "in %s line %i: unknown NEUTRFAST mode '%s'", oshf->filename, oshf->lineno, args);
    }
    beam->neutrfast = (char) _i;
    return OSH_OK;
}

static int _parse_neutrlcut(struct beam_workspace *beam, struct oshfile *oshf, char const *args) {
    float _f;
    if (sscanf(args, "%f", &_f) != 1) {
        osh_error(EX_CONFIG, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
    }
    beam->ncut = _f;
    return OSH_OK;
}

static int _parse_nstat(struct beam_workspace *beam, struct oshfile *oshf, char const *args) {
    int _i[2];
    if (sscanf(args, "%i %i", &_i[0], &_i[1]) > 2) {
        osh_error(EX_CONFIG, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
    }
    beam->nstat = (size_t) _i[0];
    beam->nsave = (size_t) _i[1];
    return OSH_OK;
}

static int _parse_nucre(struct beam_workspace *beam, struct oshfile *oshf, char const *args) {
    int _i;
    if (sscanf(args, "%i", &_i) != 1) {
        osh_error(EX_CONFIG, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
    }
    beam->nuclear = (char) _i;
    if (beam->nuclear > 1 || beam->nuclear < 0) {
        osh_error(EX_CONFIG, "in %s line %i: invalid NUCRE mode '%i'", oshf->filename, oshf->lineno, beam->nuclear);
    }
    return OSH_OK;
}

static int _parse_rndseed(struct beam_workspace *beam, struct oshfile *oshf, char const *args) {
    int _i;
    if (sscanf(args, "%i", &_i) != 1) {
        osh_error(EX_CONFIG, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
    }
    beam->rndseed = _i;
    return OSH_OK;
}

static int _parse_stragg(struct beam_workspace *beam, struct oshfile *oshf, char const *args) {
    int _i;
    if (sscanf(args, "%i", &_i) != 1) {
        osh_error(EX_CONFIG, "in %s line %i: unknown STRAGG mode '%s'", oshf->filename, oshf->lineno, args);
    }
    beam->straggl = (char) _i;
    if (beam->straggl > OSH_BEAM_STRAGG_VAVILOV || beam->straggl < OSH_BEAM_STRAGG_OFF) {
        osh_error(EX_CONFIG, "in %s line %i: invalid STRAGG mode '%i'", oshf->filename, oshf->lineno, beam->straggl);
    }
    return OSH_OK;
}

static int _parse_tmax0(struct beam_workspace *beam, struct oshfile *oshf, char const *args) {
    float _f[2];

    /* TMAX0 takes up to two values: primary energy/momentum and spread.
     *
     * Sign convention (inherited from SHIELD-HIT12A):
     *   positive value → total kinetic energy [MeV]
     *   negative value → total momentum [MeV/c] (stored as fabs)
     *
     * Only one of t0/p0 (and tsigma/psigma) is set here — whichever the file
     * specifies. The post-parse step derives the other using the particle mass.
     * Distinguishing which was given is done by checking which field is non-zero. */
    _f[0] = 0.0f;
    _f[1] = 0.0f;
    if (sscanf(args, "%f %f", &_f[0], &_f[1]) < 1) {
        osh_error(EX_CONFIG, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
    }

    if (_f[0] < 0.0f) {
        beam->spots[0].p0 = fabs((double) _f[0]); /* momentum given */
        beam->spots[0].t0 = 0.0;
    } else {
        beam->spots[0].t0 = (double) _f[0]; /* energy given */
        beam->spots[0].p0 = 0.0;
        if (_f[0] < OSH_BEAM_TMIN) {
            osh_error(EX_CONFIG,
                      "in %s line %i: TMAX0 %.4f MeV is below transport threshold %.4f MeV",
                      oshf->filename,
                      oshf->lineno,
                      (double) _f[0],
                      OSH_BEAM_TMIN);
        }
    }

    if (_f[1] < 0.0f) {
        beam->spots[0].psigma = fabs((double) _f[1]); /* momentum spread given */
        beam->spots[0].tsigma = 0.0;
    } else {
        beam->spots[0].tsigma = (double) _f[1]; /* energy spread given */
        beam->spots[0].psigma = 0.0;
    }

    return OSH_OK;
}

static int _parse_tcut0(struct beam_workspace *beam, struct oshfile *oshf, char const *args) {
    float _f[2];

    /* TCUT0 sets lower and upper transport energy cutoffs [MeV/nucleon].
     * Only the lower cutoff is stored; the upper is not kept because it is
     * effectively emax (the maximum beam energy), derived in the post-parse step. */
    _f[0] = 0.0f;
    _f[1] = 0.0f;
    if (sscanf(args, "%f %f", &_f[0], &_f[1]) > 2) {
        osh_error(EX_CONFIG, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
    }
    if (_f[0] > _f[1]) {
        osh_error(EX_CONFIG, "in %s line %i: TCUT0 upper bound must be >= lower bound", oshf->filename, oshf->lineno);
    }
    beam->tcut = fabs(_f[0]);
    return OSH_OK;
}

static int _parse_usebmod(struct beam_workspace *beam, struct oshfile *oshf, char const *args) {
    float _f;
    char tmpstr[256];
    char *_path = NULL;
    if (sscanf(args, "%f %255s", &_f, tmpstr) > 2) {
        osh_error(EX_CONFIG, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
    }
    osh_relative_path_to_file(&_path, beam->wdir, tmpstr);
    /* TODO: osh_beam_rifi_load not yet declared in a header */
    osh_warn("in %s line %i: USEBMOD parsed but rifi loader not yet implemented", oshf->filename, oshf->lineno);
    free(_path);
    return OSH_OK;
}

static int _parse_usecbeam(struct beam_workspace *beam, struct oshfile *oshf, char const *args) {
    char tmpstr[256];
    char *_path = NULL;
    size_t len;
    if (sscanf(args, "%255s", tmpstr) != 1) {
        osh_error(EX_CONFIG, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
    }
    osh_relative_path_to_file(&_path, beam->wdir, tmpstr);
    len = strlen(_path);
    free(beam->fname_spotlist);
    beam->fname_spotlist = (char *) malloc(len + 1);
    if (!beam->fname_spotlist) {
        osh_alloc_failed("beam->fname_spotlist");
    }
    memcpy(beam->fname_spotlist, _path, len + 1);
    osh_info("USECBEAM enabled: queued external spotlist %s", beam->fname_spotlist);
    free(_path);
    beam->beam_mode = OSH_BEAM_MODE_SOBP;
    return OSH_OK;
}

static int _parse_useparlev(struct beam_workspace *beam, struct oshfile *oshf, char const *args) {
    char tmpstr[256];
    char *_path = NULL;
    if (sscanf(args, "%255s", tmpstr) != 1) {
        osh_error(EX_CONFIG, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
    }
    osh_relative_path_to_file(&_path, beam->wdir, tmpstr);
    /* TODO: osh_beam_parlev_load not yet declared in a header */
    osh_warn("in %s line %i: USEPARLEV parsed but parlev loader not yet implemented", oshf->filename, oshf->lineno);
    free(_path);
    return OSH_OK;
}
