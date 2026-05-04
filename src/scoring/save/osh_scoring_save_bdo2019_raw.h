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
#define OSH_SCORING_BDO2019_ENDIAN_BIG    "MM"
#define OSH_SCORING_BDO2019_ENDIAN_LITTLE "II"

/* NumPy dtype endian prefix (1 char): prepended to pltype strings in token headers.
 * Must use NumPy byte-order characters so np.fromfile(dtype=...) can parse them. */
#define OSH_SCORING_BDO2019_DTYPE_ENDIAN_BIG    ">"
#define OSH_SCORING_BDO2019_DTYPE_ENDIAN_LITTLE "<"

#define OSH_SCORING_BDO2019_FORMAT_ID 2

/* ---- BDO 2019 tag IDs --------------------------------------------------- */

enum osh_scoring_bdo2019_tag_id {
    OSHBDO_SHVERSION = 0x0000,
    OSHBDO_SHBUILDDATE,
    OSHBDO_FILEDATE,
    OSHBDO_USER,
    OSHBDO_HOST,
    OSHBDO_FORMAT,

    OSHBDO_JPART0 = 0xCB00,
    OSHBDO_APRO0,
    OSHBDO_ZPRO0,
    OSHBDO_BEAMX,
    OSHBDO_BEAMY,
    OSHBDO_BEAMZ,
    OSHBDO_SIGMAX,
    OSHBDO_SIGMAY,
    OSHBDO_TMAX0,
    OSHBDO_SIGMAT0,
    OSHBDO_BEAMTHETA,
    OSHBDO_BEAMPHI,
    OSHBDO_BEAMDIVX,
    OSHBDO_BEAMDIVY,
    OSHBDO_BEAMDIVK,
    OSHBDO_TMAX0MEV,
    OSHBDO_TMAX0AMU,
    OSHBDO_TMAX0NUC,

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
    OSHBDO_EXT_NPROJ,
    OSHBDO_EXT_PTVDOSE,
    OSHBDO_IXFIRS,

    OSHBDO_CT_ANG = 0xCE00,
    OSHBDO_CT_ICNT,
    OSHBDO_CT_LEN,

    OSHBDO_GEO_TYPE = 0xE000,
    OSHBDO_GEO_NAME,
    OSHBDO_GEO_P,
    OSHBDO_GEO_Q,
    OSHBDO_GEO_N,
    OSHBDO_GEO_ROT,
    OSHBDO_GEO_VOL,
    OSHBDO_GEO_ZONES,
    OSHBDO_GEO_NEQGRID,
    OSHBDO_GEO_UNITS,
    OSHBDO_GEO_UNITIDS,

    OSHBDO_EST_FILENAME = 0xEE00,
    OSHBDO_EST_COUNT,
    OSHBDO_EST_NPAGES,
    OSHBDO_EST_RESCALE_NSTAT,

    OSHBDO_PAG_TYPE = 0xDD30,
    OSHBDO_PAG_COUNT,
    OSHBDO_PAG_NORMALIZE,
    OSHBDO_PAG_RESCALE,
    OSHBDO_PAG_OFFSET,
    OSHBDO_PAG_MEDIUM_TRANSP,
    OSHBDO_PAG_MEDIUM_SCORE,
    OSHBDO_PAG_UNITIDS,

    OSHBDO_PAG_DATA = 0xDDBB,
    OSHBDO_PAG_DATA_UNIT,

    OSHBDO_PAG_DIF_SET = 0xDDD0,
    OSHBDO_PAG_DIF_TYPE,
    OSHBDO_PAG_DIF_START,
    OSHBDO_PAG_DIF_STOP,
    OSHBDO_PAG_DIF_SIZE,
    OSHBDO_PAG_DIF_UNITS,

    OSHBDO_PAG_SETTINGS_NAME = 0xDDE0,

    OSHBDO_PAG_FILTER_NAME = 0xDDF0,
    OSHBDO_PAG_FILTER_NRULES,
    OSHBDO_PAG_FILTER_EMIN,
    OSHBDO_PAG_FILTER_EMAX,

    OSHBDO_RT_NSTAT = 0xAA00,
    OSHBDO_RT_TIME,
    OSHBDO_RT_TIMESIM,

    OSHBDO_COMMENT = 0xFFCC,
    OSHBDO_DEBUG,
    OSHBDO_ERROR
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
