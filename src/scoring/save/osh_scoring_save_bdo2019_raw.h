#ifndef OSH_SCORING_SAVE_BDO2019_RAW_H
#define OSH_SCORING_SAVE_BDO2019_RAW_H

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include "openshieldhit/status.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---- BDO 2019 constants ------------------------------------------------- */

#define OSH_SCORING_BDO2019_MAGIC_NUMBER "xSH12A"

/* NumPy-style payload type strings, as used by the existing BDO readers. */
#define OSH_SCORING_BDO2019_PL_TYPE_NONE "V"
#define OSH_SCORING_BDO2019_PL_TYPE_CHAR "S"
#define OSH_SCORING_BDO2019_PL_TYPE_LLUINT "u8"
#define OSH_SCORING_BDO2019_PL_TYPE_LLSINT "i8"
#define OSH_SCORING_BDO2019_PL_TYPE_DOUBLE "f8"

/* Preamble endian field (2 bytes): SH12A convention ("II" = Intel/LE, "MM" = Motorola/BE). */
#define OSH_SCORING_BDO2019_ENDIAN_BIG "MM"
#define OSH_SCORING_BDO2019_ENDIAN_LITTLE "II"

/* NumPy dtype endian prefix (1 char): prepended to pltype strings in token headers.
 * Must use NumPy byte-order characters so np.fromfile(dtype=...) can parse them. */
#define OSH_SCORING_BDO2019_DTYPE_ENDIAN_BIG ">"
#define OSH_SCORING_BDO2019_DTYPE_ENDIAN_LITTLE "<"

#define OSH_SCORING_BDO2019_FORMAT_ID 2

/* ---- BDO 2019 tag IDs --------------------------------------------------- */

enum osh_scoring_bdo2019_tag_id {
    /* Group 0x0000 - 0x00FF: miscellaneous file info. */
    OSHBDO_SHVERSION = 0x0000, /* [char*] OpenShieldHIT / SHIELD-HIT12A version string. */
    OSHBDO_SHBUILDDATE,        /* [char*] build date. */
    OSHBDO_FILEDATE,           /* [char*] BDO file creation date, RFC 2822 compliant. */
    OSHBDO_USER,               /* [char*] optional login name. */
    OSHBDO_HOST,               /* [char*] optional host where this file was created. */
    OSHBDO_FORMAT,             /* [int] optional flavour/version ID to help BDO parsers. */

    /* Group 0xCB00 - 0xCBFF: beam configuration. */
    OSHBDO_JPART0 = 0xCB00, /* [int] primary particle ID in SH12A JPART terminology (32768 = INVALID). */
    OSHBDO_APRO0,           /* [int] projectile nucleon number A; only written if A > 0. */
    OSHBDO_ZPRO0,           /* [int] projectile charge Z; may be negative (32768 = INVALID). */
    OSHBDO_BEAMX,           /* [float] beam start position, X coordinate. */
    OSHBDO_BEAMY,           /* [float] beam start position, Y coordinate. */
    OSHBDO_BEAMZ,           /* [float] beam start position, Z coordinate. */
    OSHBDO_SIGMAX,          /* [float] lateral beam extension in X direction. */
    OSHBDO_SIGMAY,          /* [float] lateral beam extension in Y direction. */
    OSHBDO_TMAX0,           /* [float] initial projectile energy (unit depends on projectile type). */
    OSHBDO_SIGMAT0,         /* [float] primary-particle energy spread. */
    OSHBDO_BEAMTHETA,       /* [float] beam polar angle. */
    OSHBDO_BEAMPHI,         /* [float] beam azimuth angle. */
    OSHBDO_BEAMDIVX,        /* [float] beam divergence, X coordinate. */
    OSHBDO_BEAMDIVY,        /* [float] beam divergence, Y coordinate. */
    OSHBDO_BEAMDIVK,        /* [float] beam divergence focus. */
    OSHBDO_TMAX0MEV,        /* [double] initial projectile energy, always in MeV. */
    OSHBDO_TMAX0AMU,        /* [double] initial projectile energy in MeV/amu; only written if mass > 1e-6 u. */
    OSHBDO_TMAX0NUC,        /* [double] initial projectile energy in MeV/nucleon; only written if A > 0. */

    /* Group 0xCC00 - 0xCCFF: run/physics configuration. */
    OSHBDO_DELE = 0xCC00,
    OSHBDO_DEMIN,
    OSHBDO_ITYPST,
    OSHBDO_ITYPMS,
    OSHBDO_OLN,
    OSHBDO_INUCRE,
    OSHBDO_IEMTRANS,
    OSHBDO_IEXTSPEC,
    OSHBDO_INTRFAST,
    OSHBDO_INTRSLOW,
    OSHBDO_APZLSCL,
    OSHBDO_IOFFSET,
    OSHBDO_IRIFIMC,
    OSHBDO_IRIFITRANS,
    OSHBDO_IRIFIZONE,
    OSHBDO_EXT_NPROJ, /* [int] requested number of projectiles, not necessarily the simulated count. */
    OSHBDO_EXT_PTVDOSE,
    OSHBDO_IXFIRS,

    /* Group 0xCE00 - 0xCEFF: CT-specific tags. */
    OSHBDO_CT_ANG = 0xCE00, /* [double] couch and gantry angles. */
    OSHBDO_CT_ICNT,         /* [int/double] CT grid counts, three values. */
    OSHBDO_CT_LEN,          /* [int/double] CT grid lengths, three values. */

    /* Group 0xE000 - 0xE0FF: scoring geometry. */
    OSHBDO_GEO_TYPE = 0xE000, /* [char*] geometry type name, e.g. MSH, CYL, ZONE, VOXSCORE. */
    OSHBDO_GEO_NAME,          /* [char*] user-given geometry name. */
    OSHBDO_GEO_P,             /* [double] geometry start values, e.g. xmin, ymin, zmin. */
    OSHBDO_GEO_Q,             /* [double] geometry stop values, e.g. xmax, ymax, zmax. */
    OSHBDO_GEO_N,             /* [int] number of bins along each geometry axis. */
    OSHBDO_GEO_ROT,           /* [future] geometry rotation. */
    OSHBDO_GEO_VOL,           /* [double] bin volume(s) in cm3; may become a list for heterogeneous bins. */
    OSHBDO_GEO_ZONES,         /* [int] single GEMCA zone or list of zone IDs. */
    OSHBDO_GEO_NEQGRID,       /* [double] non-equidistant grid coordinates; tag only used if set. */
    OSHBDO_GEO_UNITS,         /* [char*] ASCII string of semicolon-separated units along each geometry axis. */
    OSHBDO_GEO_UNITIDS,       /* [int] unit IDs, one unit along each geometry axis. */

    /* Group 0xEE00 - 0xEEFF: estimator/output metadata. */
    OSHBDO_EST_FILENAME = 0xEE00, /* [char*] output filename for this estimator. */
    OSHBDO_EST_COUNT,             /* [int] unique estimator/output number, starting at 0. */
    OSHBDO_EST_NPAGES,            /* [int] number of pages in this estimator/output. */
    OSHBDO_EST_RESCALE_NSTAT,     /* [double] estimator per-particle rescale; absent or 1 if no rescaling. */
    /* Data are not multiplied by this tag; readers apply it when page normalisation says so. */

    /* Group 0xDD30 - 0xDDFF: page-specific metadata and data. */
    OSHBDO_PAG_TYPE = 0xDD30, /* [int] score/detector type; starts a new page block. */
    OSHBDO_PAG_COUNT,         /* [int] page number, starting at 0; total page count is OSHBDO_EST_NPAGES. */
    OSHBDO_PAG_NORMALIZE,     /* [int] page postprocess/normalisation mode. */
                              /* Result X, with units OSHBDO_PAG_DATA_UNIT, is interpreted as:
                               *   0: X = x_0                              GEOMAP/raw page
                               *   1: X = sum_j x_j                        SUM/COUNT-like page
                               *   2: X = (sum_j x_j) / (sum_j I_j)        NORM page
                               *   3: X = (sum_j x_j * I_j) / (sum_j I_j)  AVER page
                               *   4: X = [x_0, x_1, ... x_j]              APPEND/MCPL-like page
                               * where j indexes independent runs/files and I_j is each file's nstat. */
    OSHBDO_PAG_RESCALE,       /* [double] if present and != 1, data were multiplied by this factor. */
    OSHBDO_PAG_OFFSET,        /* [double] if present and != 0, data were offset by this value. */
    OSHBDO_PAG_MEDIUM_TRANSP, /* [future] ASCII string for transport medium. */
    OSHBDO_PAG_MEDIUM_SCORE,  /* [future] ASCII string for scoring medium. */
    OSHBDO_PAG_UNITIDS,       /* [int] unit IDs for requested scored quantity, Diff1, and Diff2. */

    /* Page data. */
    OSHBDO_PAG_DATA = 0xDDBB, /* [double] page data block; terminates a page in legacy BDO readers. */
    OSHBDO_PAG_DATA_UNIT, /* [char*] unit of OSHBDO_PAG_DATA values, including differential divisions after postprocess.
                           */

    /* Page differential-axis metadata. */
    OSHBDO_PAG_DIF_SET = 0xDDD0, /* [int] differential binning flag: 1 = linear, -1 = log10/log axis, 0 = none. */
    OSHBDO_PAG_DIF_TYPE,         /* [int] one or two differential quantity type IDs. */
    OSHBDO_PAG_DIF_START,        /* [double] one or two lower axis bounds, for Diff1 and optional Diff2. */
    OSHBDO_PAG_DIF_STOP,         /* [double] one or two upper axis bounds, for Diff1 and optional Diff2. */
    OSHBDO_PAG_DIF_SIZE,         /* [int] one or two bin counts, for Diff1 and optional Diff2. */
    OSHBDO_PAG_DIF_UNITS,        /* [char*] semicolon-separated component units: value;Diff1;Diff2. */

    /* Settings data attached to page. */
    OSHBDO_PAG_SETTINGS_NAME = 0xDDE0, /* [char*] space-delimited list of attached Settings names. */

    /* Filter data attached to page. */
    OSHBDO_PAG_FILTER_NAME = 0xDDF0, /* [char*] space-delimited list of attached Filter names. */
    OSHBDO_PAG_FILTER_NRULES,        /* [int] number of filter rules applied. */
    OSHBDO_PAG_FILTER_EMIN,          /* [double] lower filter energy threshold, Emin. */
    OSHBDO_PAG_FILTER_EMAX,          /* [double] upper filter energy threshold, Emax. */

    /* Group 0xAA00 - 0xAAFF: runtime variables. */
    OSHBDO_RT_NSTAT = 0xAA00, /* [int] number of actually simulated primary particles. */
    OSHBDO_RT_TIME,           /* [double] optional total runtime in seconds. */
    OSHBDO_RT_TIMESIM,        /* [double] optional simulation time in seconds, excluding initialisation/finalisation. */
    OSHBDO_RT_COMPLETENESS,   /**< Partial-result honesty label (issue #193/#195): "exact" for a
                                   family-complete result, "families_pending=…" for a mid-run
                                   snapshot taken before every secondary family was drained. */

    /* Group 0xFF00 - 0xFFFF: diagnostics / informational payloads. */
    OSHBDO_COMMENT = 0xFFCC, /* [char*] optional human-readable comment. */
    OSHBDO_DEBUG,            /* [char*] optional debug payload. */
    OSHBDO_ERROR             /* [char*] optional error payload. */
};

/* ---- Raw tag layout ----------------------------------------------------- */

/**
 * @brief One BDO 2019 token header.
 *
 * @details
 * This layout is fixed-width on disk: 8-byte tag ID, 8-byte payload type
 * string, 8-byte payload element count. Existing readers expect this exact
 * 24-byte header shape.
 */
struct osh_scoring_bdo2019_tag {
    uint64_t tag;
    char pltype[8];
    uint64_t len;
};

/* ---- Raw writer API ----------------------------------------------------- */

/*
 * Each write_token_* function emits one BDO 2019 tag: a 24-byte header
 * (osh_scoring_bdo2019_tag) followed immediately by the payload bytes.
 * Callers should use these primitives rather than writing headers manually.
 */

/** @brief Write the fixed BDO 2019 file preamble (magic number + version tag). */
enum osh_status osh_scoring_bdo2019_write_preamble(FILE *fp, char const *version_string);

/** @brief Write a string-payload token (payload type "S"). */
enum osh_status osh_scoring_bdo2019_write_token_str(FILE *fp, uint64_t tag_id, char const *str);

/** @brief Write a 64-bit signed integer array token (payload type "i8"). */
enum osh_status
osh_scoring_bdo2019_write_token_llint(FILE *fp, uint64_t tag_id, long long int const *values, size_t nvalues);

/** @brief Write a 32-bit signed integer array token (payload type "i4"). */
enum osh_status osh_scoring_bdo2019_write_token_int(FILE *fp, uint64_t tag_id, int const *values, size_t nvalues);

/** @brief Write a double-precision float array token (payload type "f8"). */
enum osh_status osh_scoring_bdo2019_write_token_double(FILE *fp, uint64_t tag_id, double const *values, size_t nvalues);

/** @brief Write a single-precision float array token (payload type "f4"). */
enum osh_status osh_scoring_bdo2019_write_token_float(FILE *fp, uint64_t tag_id, float const *values, size_t nvalues);

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_SAVE_BDO2019_RAW_H */
