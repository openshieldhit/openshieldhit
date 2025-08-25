#include "osh_beam.h"
#include "osh_beam_parse.h"
#include "osh_beam_parse_keys.h"

#include "osh_logger.h"

/* Forward declarations of handler functions */
static int _parse_apcorr(struct beam_workspace *beam, struct oshfile *oshf, const char *args);
static int _parse_beamdir(struct beam_workspace *beam, struct oshfile *oshf, const char *args);
static int _parse_beamdiv(struct beam_workspace *beam, struct oshfile *oshf, const char *args);
static int _parse_beampos(struct beam_workspace *beam, struct oshfile *oshf, const char *args);
static int _parse_beamsad(struct beam_workspace *beam, struct oshfile *oshf, const char *args);
static int _parse_beamsigma(struct beam_workspace *beam, struct oshfile *oshf, const char *args);
static int _parse_bmodmc(struct beam_workspace *beam, struct oshfile *oshf, const char *args);
static int _parse_bmodtrans(struct beam_workspace *beam, struct oshfile *oshf, const char *args);
static int _parse_deltae(struct beam_workspace *beam, struct oshfile *oshf, const char *args);
static int _parse_demin(struct beam_workspace *beam, struct oshfile *oshf, const char *args);
static int _parse_emtrans(struct beam_workspace *beam, struct oshfile *oshf, const char *args);
static int _parse_extspec(struct beam_workspace *beam, struct oshfile *oshf, const char *args);
static int _parse_hiproj(struct beam_workspace *beam, struct oshfile *oshf, const char *args);
static int _parse_jpart0(struct beam_workspace *beam, struct oshfile *oshf, const char *args);
static int _parse_makeln(struct beam_workspace *beam, struct oshfile *oshf, const char *args);
static int _parse_mscat(struct beam_workspace *beam, struct oshfile *oshf, const char *args);
static int _parse_neutrfast(struct beam_workspace *beam, struct oshfile *oshf, const char *args);
static int _parse_neutrlcut(struct beam_workspace *beam, struct oshfile *oshf, const char *args);
static int _parse_nstat(struct beam_workspace *beam, struct oshfile *oshf, const char *args);
static int _parse_nucre(struct beam_workspace *beam, struct oshfile *oshf, const char *args);
static int _parse_rndseed(struct beam_workspace *beam, struct oshfile *oshf, const char *args);
static int _parse_stragg(struct beam_workspace *beam, struct oshfile *oshf, const char *args);
static int _parse_tmax0(struct beam_workspace *beam, struct oshfile *oshf, const char *args);
static int _parse_tcut0(struct beam_workspace *beam, struct oshfile *oshf, const char *args);
static int _parse_usebmod(struct beam_workspace *beam, struct oshfile *oshf, const char *args);
static int _parse_usecbeam(struct beam_workspace *beam, struct oshfile *oshf, const char *args);
static int _parse_useparlev(struct beam_workspace *beam, struct oshfile *oshf, const char *args);

/* ... */

/* Dispatch entry */
struct beam_dispatch_entry {
    const char *key;
    int (*handler)(struct beam_workspace *beam, struct oshfile *oshf, const char *args);
};

static struct beam_dispatch_entry beam_dispatch_table[] = {
    { OSH_BEAM_KEY_APCORR,     _parse_apcorr },
    { OSH_BEAM_KEY_BEAMDIR,    _parse_beamdir },
    { OSH_BEAM_KEY_BEAMDIV,    _parse_beamdiv },
    { OSH_BEAM_KEY_BEAMPOS,    _parse_beampos },
    { OSH_BEAM_KEY_BEAMSAD,    _parse_beamsad },
    { OSH_BEAM_KEY_BEAMSIGMA,  _parse_beamsigma },
    { OSH_BEAM_KEY_BMODMC,     _parse_bmodmc },
    { OSH_BEAM_KEY_BMODTRANS,  _parse_bmodtrans },
    { OSH_BEAM_KEY_DELTAE,     _parse_deltae },
    { OSH_BEAM_KEY_DEMIN,      _parse_demin },
    { OSH_BEAM_KEY_EMTRANS,    _parse_emtrans },
    { OSH_BEAM_KEY_EXTSPEC,    _parse_extspec },
    { OSH_BEAM_KEY_HIPROJ,     _parse_hiproj },
    { OSH_BEAM_KEY_JPART0,     _parse_jpart0 },
    { OSH_BEAM_KEY_MAKELN,     _parse_makeln },
    { OSH_BEAM_KEY_MSCAT,      _parse_mscat },
    { OSH_BEAM_KEY_NEUTRFAST,  _parse_neutrfast },
    { OSH_BEAM_KEY_NEUTRLCUT,  _parse_neutrlcut },
    { OSH_BEAM_KEY_NSTAT,      _parse_nstat },
    { OSH_BEAM_KEY_NUCRE,      _parse_nucre },
    { OSH_BEAM_KEY_RNDSEED,    _parse_rndseed },
    { OSH_BEAM_KEY_STRAGG,     _parse_stragg },
    { OSH_BEAM_KEY_TMAX0,      _parse_tmax0 },
    { OSH_BEAM_KEY_TCUT0,      _parse_tcut0 },
    { OSH_BEAM_KEY_USEBMOD,    _parse_usebmod },
    { OSH_BEAM_KEY_USECBEAM,   _parse_usecbeam },
    { OSH_BEAM_KEY_USEPARLEV,  _parse_useparlev },
    { NULL,                    NULL }  /* Sentinel */
};


/* Main parser entry point */
int osh_beam_parse(struct oshfile *oshf, struct beam_workspace *beam) {
    char *lline = NULL, *key = NULL, *args = NULL;
    int lineno;
    int i;

    while (osh_readline_key(oshf, &lline, &key, &args, &lineno) != -1) {
        int found = 0;
        for (i = 0; beam_dispatch_table[i].key != NULL; i++) {
            if (strcmp(beam_dispatch_table[i].key, key) == 0) {
                beam_dispatch_table[i].handler(beam, args, lineno);
                found = 1;
                break;
            }
        }
        if (!found) {
            osh_warn("Line %d: Unknown key '%s'\n", lineno, key);
        }
        free(lline);
    }
    return 0;
}


static int _parse_apcorr(struct beam_workspace *beam, struct oshfile *oshf, const char *args) {
    return 0;
}


static int _parse_beamdir(struct beam_workspace *beam, struct oshfile *oshf, const char *args) {
    float _f[2];
    if (sscanf(args, "%f %f", &_f[0], &_f[1]) != 2) {
        osh_err(EX_CONFIG, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
    }
    if (_f[0] <
        0.0) osh_err(EX_CONFIG, "in %s line %i: theta must be within [0:180] deg", oshf->filename, oshf->lineno);
    if (_f[0] > 180.0) osh_err(EX_CONFIG, "in %s line %i: theta must be within [0:180] deg", oshf->filename,
                               oshf->lineno);
    if (_f[1] < 0.0) osh_err(EX_CONFIG, "in %s line %i: phi must be within [0:360] deg", oshf->filename, oshf->lineno);
    if (_f[1] >
        360.0) osh_err(EX_CONFIG, "in %s line %i: phi must be within [0:360] deg", oshf->filename, oshf->lineno);
    //beam->spot0->theta = (double) _f[0] * M_PI / 180.0;
    //beam->spot0->phi = (double) _f[1] * M_PI / 180.0;
    return 0;
}


static int _parse_beamdiv(struct beam_workspace *beam, struct oshfile *oshf, const char *args) {
    float _f[3];
    if (sscanf(args, "%f %f %f", &_f[0], &_f[1], &_f[2]) > 3) {
        osh_err(EX_CONFIG, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
    }
    beam->spot0->div[0] = (double) _f[0] * 0.001;
    beam->spot0->div[1] = (double) _f[1] * 0.001;
    beam->spot0->focus = (double) _f[2];
    if (fabs(_f[0]) > 0.0 || fabs(_f[1]) > 0.0) beam->spot0->use_div = 1;
    return 0;
}


static int _parse_beampos(struct beam_workspace *beam, struct oshfile *oshf, const char *args) {
    float _f[3];
    int i;
    if (sscanf(args, "%f %f %f", &_f[0], &_f[1], &_f[2]) != 3) {
        osh_err(EX_CONFIG, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
    }
    for (i = 0; i < 3; i++) {
        beam->spot0->p[i] = (double) _f[i];
    }
    return 0;
}


static int _parse_beamsad(struct beam_workspace *beam, struct oshfile *oshf, const char *args) {
    float _f[2];
    int i = sscanf(args, "%f %f", &_f[0], &_f[1]);
    switch (i) {
    case 1:
        beam->spot0->sad[0] = _f[0];
        beam->spot0->sad[1] = _f[0];
        break;
    case 2:
        beam->spot0->sad[0] = _f[0];
        beam->spot0->sad[1] = _f[1];
        break;
    default:
        osh_err(EX_CONFIG, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
    }
    if ((beam->spot0->sad[0] > 0.0) && (beam->spot0->sad[1] > 0.0)) beam->spot0->use_sad = 1;
    else {
        osh_err(EX_CONFIG, "In %s line %i: SAD must be > 0.0.", oshf->filename, oshf->lineno);
    }
    return 0;
}


static int _parse_beamsigma(struct beam_workspace *beam, struct oshfile *oshf, const char *args) {
    float _f[2];
    if (sscanf(args, "%f %f", &_f[0], &_f[1]) > 2) {
        osh_err(EX_CONFIG, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
    }
    if (_f[0] < 0.0 && _f[1] < 0.0) {
        beam->spot0->shape = OSH_BEAM_SHAPE_SQUARE;
        beam->spot0->size[0] = fabs(_f[0]);
        beam->spot0->size[1] = fabs(_f[1]);
    }
    else if (_f[0] >= 0.0 && _f[1] < 0.0) {
        beam->spot0->shape = OSH_BEAM_SHAPE_CIRCULAR;
        beam->spot0->size[0] = fabs(_f[1]);
        beam->spot0->size[1] = 0.0;
    }
    else if ((_f[0] >= 0.0 && _f[1] > 0.0) ||  (_f[0] > 0.0 && _f[1] >= 0.0)) {
        beam->spot0->shape = OSH_BEAM_SHAPE_GAUSSIAN;
        beam->spot0->size[0] = _f[0];
        beam->spot0->size[1] = _f[1];
    }
    else {
        beam->spot0->shape = OSH_BEAM_SHAPE_PENCIL;
        beam->spot0->size[0] = 0.0;
        beam->spot0->size[1] = 0.0;
    }
    return 0;
}


static int _parse_bmodmc(struct beam_workspace *beam, struct oshfile *oshf, const char *args) {
    if (sscanf(args, "%c", &(beam->rifi->mode_mc)) != 1) {
        osh_err(EX_CONFIG, "in %s line %i: unknown BMODMC mode '%s'", oshf->filename, oshf->lineno, args);
    }
    return 0;
}


static int _parse_bmodtrans(struct beam_workspace *beam, struct oshfile *oshf, const char *args) {
    if (sscanf(args, "%c", &(beam->rifi->mode_trans)) != 1) {
        osh_err(EX_CONFIG, "in %s line %i: unknown BMODMC mode '%s'", oshf->filename, oshf->lineno, args);
    }
    osh_warn("in %s line %i: BMODTRANS deprecated, will be ignored.\n", oshf->filename, oshf->lineno);
    return 0;
}


static int _parse_deltae(struct beam_workspace *beam, struct oshfile *oshf, const char *args) {
    float _f[1];
    if (sscanf(args, "%f", &_f[0]) != 1) {
        osh_err(EX_CONFIG, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
    }
    beam->deltae = _f[0];
    return 0;
}


static int _parse_demin(struct beam_workspace *beam, struct oshfile *oshf, const char *args) {
    float _f[1];
    if (sscanf(args, "%f", &_f[0]) != 1) {
        osh_err(EX_CONFIG, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
    }
    beam->demin = _f[0];
    return 0;
}


static int _parse_emtrans(struct beam_workspace *beam, struct oshfile *oshf, const char *args) {
    if (sscanf(args, "%c", &(beam->emtrans)) != 1) {
        osh_err(EX_CONFIG, "in %s line %i: unknown EMTRANS mode '%s'", oshf->filename, oshf->lineno, args);
    }
    return 0;
}


static int _parse_extspec(struct beam_workspace *beam, struct oshfile *oshf, const char *args) {
    osh_err(EX_CONFIG, "in %s line %i: EXTSPEC not implemented", oshf->filename, oshf->lineno);
    return 0;
}


static int _parse_hiproj(struct beam_workspace *beam, struct oshfile *oshf, const char *args) {
    float _f[2];
    if (sscanf(args, "%f %f", &_f[0], &_f[1]) != 2) {
        osh_err(EX_CONFIG, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
    }
    beam->spot0->part->a = (int) _f[0];
    beam->spot0->part->z = (int) _f[1];
    return 0;
}


static int _parse_jpart0(struct beam_workspace *beam, struct oshfile *oshf, const char *args) {
    int _i[1];
    if (sscanf(args, "%i", &_i[0]) != 1) {
        osh_err(EX_CONFIG, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
    }
    if ((_i[0] > 25) || (_i[0] < 1)) {
        osh_err(EX_CONFIG, "in %s line %i: invalid JPART0 number '%i'", oshf->filename, oshf->lineno, _i[0]);
    }
    beam->spot0->part->_jpart = _i[0];
    return 0;
}


static int _parse_makeln(struct beam_workspace *beam, struct oshfile *oshf, const char *args) {
    int _i[1];
    if (sscanf(args, "%i", &_i[0]) != 1) {
        osh_err(EX_CONFIG, "in %s line %i: unknown MAKELN mode '%s'", oshf->filename, oshf->lineno, args);
    }
    beam->makeln = (char) _i[0];
    return 0;
}


static int _parse_mscat(struct beam_workspace *beam, struct oshfile *oshf, const char *args) {
    int _i[1];
    if (sscanf(args, "%i", &_i[0]) != 1) {
        osh_err(EX_CONFIG, "in %s line %i: unknown MSCAT mode '%s'", oshf->filename, oshf->lineno, args);
    }
    beam->scatter = (char) _i[0];
    if (((beam->scatter) > OSH_BEAM_MSCAT_MOLIERE) || ((beam->scatter) < OSH_BEAM_MSCAT_OFF)) {
        osh_err(EX_CONFIG, "in %s line %i: invalid MSCAT mode '%i'", oshf->filename, oshf->lineno, beam->scatter);
    }
    return 0;
}


static int _parse_neutrfast(struct beam_workspace *beam, struct oshfile *oshf, const char *args) {
    int _i[1];
    if (sscanf(args, "%i", &_i[0]) != 1) {
        osh_err(EX_CONFIG, "in %s line %i: unknown NEUTRFAST mode '%s'", oshf->filename, oshf->lineno, args);
    }
    beam->neutrfast =  (char) _i[0];
    return 0;
}


static int _parse_neutrlcut(struct beam_workspace *beam, struct oshfile *oshf, const char *args) {
    float _f[1];
    if (sscanf(args, "%f", &_f[0]) != 1) {
        osh_err(EX_CONFIG, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
    }
    beam->oln = _f[0];
    return 0;
}


static int _parse_nstat(struct beam_workspace *beam, struct oshfile *oshf, const char *args) {
    int _i[2];
    if (sscanf(args, "%i %i", &_i[0], &_i[1]) > 2) {
        osh_err(EX_CONFIG, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
    }
    beam->nstat = (size_t) _i[0];
    beam->nsave = (size_t) _i[1];
    return 0;
}


static int _parse_nucre(struct beam_workspace *beam, struct oshfile *oshf, const char *args) {
    int _i[1];
    if (sscanf(args, "%i", &_i[0]) != 1) {
        osh_err(EX_CONFIG, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
    }
    beam->nuclear = (char) _i[0];
    if (((beam->nuclear) > 1) || ((beam->nuclear) < 0)) {
        osh_err(EX_CONFIG, "in %s line %i: invalid NUCRE mode '%i'", oshf->filename, oshf->lineno, beam->nuclear);
    }
    return 0;
}


static int _parse_rndseed(struct beam_workspace *beam, struct oshfile *oshf, const char *args) {
    int _i[1];
    if (sscanf(args, "%i", &_i[0]) != 1) {
        osh_err(EX_CONFIG, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
    }
    beam->rndseed = _i[0];
    return 0;
}


static int _parse_stragg(struct beam_workspace *beam, struct oshfile *oshf, const char *args) {
    int _i[1];
    if (sscanf(args, "%i", &_i[0]) != 1) {
        osh_err(EX_CONFIG, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
    }
    beam->straggl = (char) _i[0];
    if (((beam->straggl) > OSH_BEAM_STRAGG_VAVILOV) || ((beam->straggl) < OSH_BEAM_STRAGG_OFF)) {
        osh_err(EX_CONFIG, "in %s line %i: invalid STRAGG mode '%i'", oshf->filename, oshf->lineno, beam->straggl);
    }
    return 0;
}


static int _parse_tmax0(struct beam_workspace *beam, struct oshfile *oshf, const char *args) {
    float _f[2];
    if (sscanf(args, "%f %f", &_f[0], &_f[1]) > 2) {
        osh_err(EX_CONFIG, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
    }
    if (_f[0] < 0.0) {
        beam->spot0->pmax0 = fabs(_f[0]);
        beam->spot0->use_pmax0 = 1;
    }
    else {
        beam->spot0->tmax0 = _f[0];
        beam->spot0->use_pmax0 = 0;
        if (_f[0] < OSH_BEAM_TMIN0)
            osh_err(EX_CONFIG, "in %s line %i: TMAX0 is below transport threshold '%s'",
                    oshf->filename, oshf->lineno, args);
    }
    if (_f[1] < 0.0) {
        beam->spot0->psigma0 = fabs(_f[1]);
        beam->spot0->use_psigma0 = 1;
    }
    else {
        beam->spot0->tsigma0 = _f[1];
        beam->spot0->use_psigma0 = 0;
    }
    return 0;
}


static int _parse_tcut0(struct beam_workspace *beam, struct oshfile *oshf, const char *args) {
    float _f[2];
    if (sscanf(args, "%f %f", &_f[0], &_f[1]) > 2) {
        osh_err(EX_CONFIG, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
    }
    if (_f[0] > _f[1])
        osh_err(EX_CONFIG, "in %s line %i: TCUT0 upper bound must be larger than lower bound.",
                oshf->filename, oshf->lineno, args);
    beam->spot0->tcut0[0] = fabs(_f[0]);
    beam->spot0->tcut0[1] = fabs(_f[1]);
    return 0;
}


static int _parse_usebmod(struct beam_workspace *beam, struct oshfile *oshf, const char *args) {
    float _f[1];
    char tmpstr[256];
    char *_path = NULL;
    if (sscanf(args, "%f %s", &_f[0], tmpstr) > 2) {
        osh_err(EX_CONFIG, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
    }
    osh_relative_path_to_file(&_path, beam->wdir, tmpstr);
    osh_beam_rifi_load(_path, _f[0], beam->rifi);
    free(_path);
    return 0;
}


static int _parse_usecbeam(struct beam_workspace *beam, struct oshfile *oshf, const char *args) {
    char tmpstr[256];
    char *_path = NULL;
    if (sscanf(args, "%s", tmpstr) > 1) {
        osh_err(EX_CONFIG, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
    }
    osh_relative_path_to_file(&_path, beam->wdir, tmpstr);
    beam->_sobp_fname = calloc(strlen(_path) + 1, sizeof(char));
    if (!beam->_sobp_fname) osh_malloc_err("beam->sobp_fname");
    sprintf(beam->_sobp_fname, "%s", _path);
    free(_path);
    beam->beam_mode = OSH_BEAM_MODE_SOBP;
    return 0;
}


static int _parse_useparlev(struct beam_workspace *beam, struct oshfile *oshf, const char *args) {
    char tmpstr[256];
    char *_path = NULL;
    if (sscanf(args, "%s", tmpstr) > 1) {
        osh_err(EX_CONFIG, "in %s line %i: parse error '%s'", oshf->filename, oshf->lineno, args);
    }
    osh_relative_path_to_file(&_path, beam->wdir, tmpstr);
    osh_beam_parlev_load(_path, beam->parlev);
    free(_path);
    return 0;
}

