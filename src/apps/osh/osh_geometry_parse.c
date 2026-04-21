#include "apps/osh/osh_geometry_parse.h"

#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_geometry_parse_keys.h"
#include "common/osh_diag.h"
#include "common/osh_readline.h"
#include "gemca/osh_gemca2_defines.h"
#include "openshieldhit/geometry.h"
#include "openshieldhit/geometry_defs.h"
#include "openshieldhit/status.h"

/* ---- Internal helpers ---------------------------------------------------- */

static int _rewind_oshfile(struct oshfile *shf, struct osh_diag_sink const *diag);

/**
 * @brief Map a body-type key string to an OSH_GEOMETRY_BODY_* code.
 *
 * @details Values are identical to the internal OSH_GEMCA_BODY_* constants so
 * that osh_geometry_workspace_prepare() can use them without conversion.
 *
 * @param[in] key  Lowercase or mixed-case body type token (e.g. "sph", "rpp").
 *
 * @returns The matching OSH_GEOMETRY_BODY_* constant, or OSH_GEOMETRY_BODY_NONE
 *          if @p key is unrecognized.
 */
static int _body_type_from_key(char const *key) {
    if (strcasecmp(key, "sph") == 0) {
        return OSH_GEOMETRY_BODY_SPH;
    }
    if (strcasecmp(key, "wed") == 0) {
        return OSH_GEOMETRY_BODY_WED;
    }
    if (strcasecmp(key, "arb") == 0) {
        return OSH_GEOMETRY_BODY_ARB;
    }
    if (strcasecmp(key, "box") == 0) {
        return OSH_GEOMETRY_BODY_BOX;
    }
    if (strcasecmp(key, "vox") == 0) {
        return OSH_GEOMETRY_BODY_VOX;
    }
    if (strcasecmp(key, "rpp") == 0) {
        return OSH_GEOMETRY_BODY_RPP;
    }
    if (strcasecmp(key, "rcc") == 0) {
        return OSH_GEOMETRY_BODY_RCC;
    }
    if (strcasecmp(key, "rec") == 0) {
        return OSH_GEOMETRY_BODY_REC;
    }
    if (strcasecmp(key, "trc") == 0) {
        return OSH_GEOMETRY_BODY_TRC;
    }
    if (strcasecmp(key, "ell") == 0) {
        return OSH_GEOMETRY_BODY_ELL;
    }
    if (strcasecmp(key, "yzp") == 0) {
        return OSH_GEOMETRY_BODY_YZP;
    }
    if (strcasecmp(key, "xzp") == 0) {
        return OSH_GEOMETRY_BODY_XZP;
    }
    if (strcasecmp(key, "xyp") == 0) {
        return OSH_GEOMETRY_BODY_XYP;
    }
    if (strcasecmp(key, "pla") == 0) {
        return OSH_GEOMETRY_BODY_PLA;
    }
    if (strcasecmp(key, "rot") == 0) {
        return OSH_GEOMETRY_BODY_ROT;
    }
    if (strcasecmp(key, "cpy") == 0) {
        return OSH_GEOMETRY_BODY_CPY;
    }
    if (strcasecmp(key, "mov") == 0) {
        return OSH_GEOMETRY_BODY_MOV;
    }
    return OSH_GEOMETRY_BODY_NONE;
}

/**
 * @brief Check whether @p key is a zone boolean-expression continuation token.
 *
 * @param[in] key  Token to test.
 *
 * @returns 1 if @p key starts with +, -, |, ( or ), 0 otherwise.
 */
static int _is_zone_continuation(char const *key) {
    return key && (key[0] == '+' || key[0] == '-' || key[0] == '|' || key[0] == '(' || key[0] == ')');
}

/**
 * @brief Count body cards in the geometry file without retaining state.
 *
 * @param[in,out] shf  Open geometry file; rewound internally.
 *
 * @returns Number of body cards found, or 0 on I/O error.
 */
static size_t _count_bodies(struct oshfile *shf, struct osh_diag_sink const *diag) {
    int lineno;
    char *line = NULL;
    char *key = NULL;
    char *args = NULL;
    size_t nbody = 0u;

    if (!_rewind_oshfile(shf, diag)) {
        return 0u;
    }

    while (osh_readline_key(shf, &line, &key, &args, &lineno) > 0) {
        if (strcasecmp(key, OSH_GEO_KEY_END) != 0 && _body_type_from_key(key) != OSH_GEOMETRY_BODY_NONE) {
            ++nbody;
        }
        free(line);
        line = NULL;
    }

    free(line);
    return nbody;
}

/**
 * @brief Count zones in the geometry file without building any cold objects.
 *
 * @param[in,out] shf  Open geometry file; rewound internally.
 *
 * @returns Number of zone header cards found, or 0 on I/O error.
 */
static size_t _count_zones(struct oshfile *shf, struct osh_diag_sink const *diag) {
    int lineno;
    char *line = NULL;
    char *key = NULL;
    char *args = NULL;
    size_t nzone = 0u;
    int in_block = 0;

    if (!_rewind_oshfile(shf, diag)) {
        return 0u;
    }

    while (osh_readline_key(shf, &line, &key, &args, &lineno) > 0) {
        if (strcasecmp(key, OSH_GEO_KEY_END) != 0) {
            if (in_block == 2 && !_is_zone_continuation(key)) {
                ++nzone;
            }
        } else {
            in_block++;
        }

        if (in_block == 0) {
            in_block = 1;
        }

        free(line);
        line = NULL;
    }

    free(line);
    OSH_DIAG_INFOF(diag, "Found %zu zones in geo.dat file", nzone);
    return nzone;
}

/**
 * @brief Append string @p b to the heap-allocated string at @p *a.
 *
 * @param[in,out] a  Pointer to the destination heap string; grown as needed.
 * @param[in]     b  NUL-terminated string to append.
 *
 * @returns OSH_OK on success, OSH_ENOMEM on allocation failure.
 */
static enum osh_status _str_append(char **a, char const *b) {
    size_t la = strlen(*a);
    size_t lb = strlen(b);
    char *p = (char *) realloc(*a, la + lb + 1u);
    if (!p) {
        return OSH_ENOMEM;
    }
    memcpy(p + la, b, lb + 1u);
    *a = p;
    return OSH_OK;
}

/**
 * @brief Find a zone by exact name in the flat cold-zone array.
 *
 * @param[in]  name     Zone name to look up.
 * @param[in]  zones    Array of cold zone definitions.
 * @param[in]  nzones   Length of @p zones.
 * @param[out] idx_out  Receives the index of the matching zone.
 *
 * @returns 1 if the name is found (fills @p *idx_out), 0 otherwise.
 */
static int
_zone_index_from_name(char const *name, struct osh_geometry_zone const *zones, size_t nzones, size_t *idx_out) {
    size_t i;
    for (i = 0u; i < nzones; ++i) {
        if (zones[i].name && strcmp(zones[i].name, name) == 0) {
            *idx_out = i;
            return 1;
        }
    }
    return 0;
}

/**
 * @brief Return the next whitespace-delimited token from @p *cursor.
 *
 * @details Advances @p *cursor past the returned token. Modifies the buffer
 * in place by inserting a NUL terminator after each token.
 *
 * @param[in,out] cursor  Pointer into the buffer; updated to point past the token.
 *
 * @returns Pointer to the start of the token, or NULL when no more tokens remain.
 */
static char *_next_token(char **cursor) {
    char *start;
    char *end;

    if (!cursor || !*cursor) {
        return NULL;
    }
    start = *cursor;
    while (*start != '\0' && isspace((unsigned char) *start)) {
        start++;
    }
    if (*start == '\0') {
        *cursor = start;
        return NULL;
    }
    end = start;
    while (*end != '\0' && !isspace((unsigned char) *end)) {
        end++;
    }
    if (*end != '\0') {
        *end = '\0';
        end++;
    }
    *cursor = end;
    return start;
}

/**
 * @brief Normalize legacy numeric material names.
 *
 * @details Maps "0" to "blackhole" and "1000" to "vacuum" with a warning;
 * all other names pass through unchanged.
 *
 * @param[in] raw       Raw material name token from the file.
 * @param[in] filename  Source file name used in warning messages.
 * @param[in] lineno    Source line number used in warning messages.
 *
 * @returns Normalized name string (either @p raw or a literal constant).
 */
static char const *
_normalize_material_name(char const *raw, char const *filename, int lineno, struct osh_diag_sink const *diag) {
    if (strcmp(raw, "0") == 0) {
        OSH_DIAG_WARNF(diag,
                       "%s line %d: legacy material '0' mapped to 'blackhole'; use 'blackhole' explicitly",
                       filename,
                       lineno);
        return "blackhole";
    }
    if (strcmp(raw, "1000") == 0) {
        OSH_DIAG_WARNF(
            diag, "%s line %d: legacy material '1000' mapped to 'vacuum'; use 'vacuum' explicitly", filename, lineno);
        return "vacuum";
    }
    return raw;
}

/**
 * @brief Duplicate a normalized material name into @p z->material_name.
 *
 * @param[in,out] z         Zone to update.
 * @param[in]     raw_name  Raw material name token (will be normalized).
 * @param[in]     filename  Source file name used in warning messages.
 * @param[in]     lineno    Source line number used in warning messages.
 *
 * @returns OSH_OK on success, OSH_ENOMEM on allocation failure.
 */
static enum osh_status _assign_material(struct osh_geometry_zone *z,
                                        char const *raw_name,
                                        char const *filename,
                                        int lineno,
                                        struct osh_diag_sink const *diag) {
    char const *name = _normalize_material_name(raw_name, filename, lineno, diag);
    char *copy = strdup(name);
    if (!copy) {
        return OSH_ENOMEM;
    }
    free(z->material_name);
    z->material_name = copy;
    return OSH_OK;
}

/**
 * @brief Rewind @p shf to the beginning and reset its line counter.
 *
 * @param[in,out] shf  Open geometry file to rewind.
 *
 * @returns 1 on success, 0 on I/O error.
 */
static int _rewind_oshfile(struct oshfile *shf, struct osh_diag_sink const *diag) {
    if (fseek(shf->fp, 0L, SEEK_SET) != 0) {
        OSH_DIAG_ERRORF(diag, "Failed to rewind geometry file '%s'", shf->filename);
        return 0;
    }
    shf->lineno = 0;
    return 1;
}

/* ---- Phase 1: body section ----------------------------------------------- */

/**
 * @brief Parse the body section from file start up to the first END card.
 *
 * @details Fills ws->bodies[0..ws->nbodies-1] with body type, name, and raw
 * argument arrays. Leaves the file positioned just after the first END card.
 *
 * @param[in,out] shf  Open geometry file; rewound internally.
 * @param[in,out] ws   Geometry workspace with pre-allocated bodies array.
 *
 * @returns OSH_OK on success, OSH_EPARSE or OSH_ENOMEM on failure.
 */
static enum osh_status
_parse_bodies(struct oshfile *shf, struct osh_diag_sink const *diag, struct osh_geometry_workspace *ws) {
    char *line = NULL;
    char *key = NULL;
    char *args = NULL;
    int lineno;
    int btype_new;
    int btype = OSH_GEOMETRY_BODY_NONE;
    int nt;
    int npar = 0;
    int off = 0;
    char nstr[OSH_GEMCA_BODY_NAME_MAXLEN];
    double par[OSH_GEMCA_NARGS_MAX];
    size_t ibody = 0;
    int body_active = 0;
    enum osh_status rc = OSH_OK;

    if (!_rewind_oshfile(shf, diag)) {
        return OSH_EIO;
    }

    while (osh_readline_key(shf, &line, &key, &args, &lineno) > 0) {

        if (strcasecmp(key, OSH_GEO_KEY_END) == 0) {
            /* Finalize last body. */
            if (body_active) {
                ws->bodies[ibody].type = btype;
                ws->bodies[ibody].name = strdup(nstr);
                if (!ws->bodies[ibody].name) {
                    rc = OSH_ENOMEM;
                    goto done;
                }
                if (npar > 0) {
                    ws->bodies[ibody].a = (double *) calloc((size_t) npar, sizeof(double));
                    if (!ws->bodies[ibody].a) {
                        rc = OSH_ENOMEM;
                        goto done;
                    }
                    memcpy(ws->bodies[ibody].a, par, (size_t) npar * sizeof(double));
                }
                ws->bodies[ibody].na = npar;
            }
            free(line);
            return OSH_OK; /* File is now positioned just after the first END. */
        }

        btype_new = _body_type_from_key(key);

        if (btype_new != OSH_GEOMETRY_BODY_NONE) {
            /* Finalize the previous body, if any. */
            if (body_active) {
                ws->bodies[ibody].type = btype;
                ws->bodies[ibody].name = strdup(nstr);
                if (!ws->bodies[ibody].name) {
                    rc = OSH_ENOMEM;
                    goto done;
                }
                if (npar > 0) {
                    ws->bodies[ibody].a = (double *) calloc((size_t) npar, sizeof(double));
                    if (!ws->bodies[ibody].a) {
                        rc = OSH_ENOMEM;
                        goto done;
                    }
                    memcpy(ws->bodies[ibody].a, par, (size_t) npar * sizeof(double));
                }
                ws->bodies[ibody].na = npar;
                ibody++;
            }

            if (ibody >= ws->nbodies) {
                OSH_DIAG_ERRORF(diag, "%s line %d: too many bodies (max=%zu)", shf->filename, lineno, ws->nbodies);
                rc = OSH_EPARSE;
                goto done;
            }

            btype = btype_new;
            body_active = 1;

            if (!args) {
                OSH_DIAG_ERRORF(diag, "%s line %d: missing body name or parameters", shf->filename, lineno);
                rc = OSH_EPARSE;
                goto done;
            }
            /* First six values: name + up to 5 floats on the same line. */
            nt = sscanf(args, "%s %lf %lf %lf %lf %lf %lf", nstr, &par[0], &par[1], &par[2], &par[3], &par[4], &par[5]);
            npar = nt - 1;
            off = 6;

        } else {
            /* Continuation line: accumulate more float arguments. */
            if (!body_active) {
                /* Allow preamble/header lines before the first body card. */
                free(line);
                line = NULL;
                continue;
            }

            if ((off + 5) >= OSH_GEMCA_NARGS_MAX) {
                OSH_DIAG_ERRORF(diag, "%s line %d: too many body arguments", shf->filename, lineno);
                rc = OSH_EPARSE;
                goto done;
            }

            nt = sscanf(key, "%lf", &par[off]);
            if (nt > 0) {
                npar += nt;
            }
            nt = 0;
            if (args) {
                nt = sscanf(args,
                            "%lf %lf %lf %lf %lf",
                            &par[1 + off],
                            &par[2 + off],
                            &par[3 + off],
                            &par[4 + off],
                            &par[5 + off]);
            }
            if (nt > 0) {
                npar += nt;
            }
            off += 6;
        }

        free(line);
        line = NULL;
    }

done:
    free(line);
    return rc;
}

/* ---- Phase 2: zone section ----------------------------------------------- */

/**
 * @brief Parse the zone section between the first and second END cards.
 *
 * @details Fills ws->zones[].name and ws->zones[].expr. Precondition: file is
 * positioned just after the first END card as left by _parse_bodies().
 *
 * @param[in,out] shf  Open geometry file; read from current position.
 * @param[in,out] ws   Geometry workspace with pre-allocated zones array.
 *
 * @returns OSH_OK on success, OSH_EPARSE or OSH_ENOMEM on failure.
 */
static enum osh_status
_parse_zones(struct oshfile *shf, struct osh_diag_sink const *diag, struct osh_geometry_workspace *ws) {
    char *line = NULL;
    char *key = NULL;
    char *args = NULL;
    int lineno;
    size_t izone = 0u;
    int zone_active = 0;
    char *expr = NULL;
    enum osh_status rc = OSH_OK;

    /* Accumulator for the current zone's raw boolean expression.
     * Starts as an empty NUL-terminated string; grows via _str_append(). */
    expr = (char *) calloc(1u, sizeof(char));
    if (!expr) {
        return OSH_ENOMEM;
    }

    while (osh_readline_key(shf, &line, &key, &args, &lineno) > 0) {

        if (_is_zone_continuation(key)) {
            if (!zone_active) {
                OSH_DIAG_ERRORF(diag, "%s line %d: zone continuation before any zone header", shf->filename, lineno);
                rc = OSH_EPARSE;
                goto done;
            }
            rc = _str_append(&expr, key);
            if (rc != OSH_OK) {
                goto done;
            }

        } else {
            /* Non-continuation key: either a new zone name or the END card.
             * In both cases, finalize the expression for the current zone first. */
            if (zone_active && expr[0] != '\0') {
                if (izone >= ws->nzones) {
                    OSH_DIAG_ERRORF(diag, "%s line %d: too many zones (max=%zu)", shf->filename, lineno, ws->nzones);
                    rc = OSH_EPARSE;
                    goto done;
                }
                ws->zones[izone].expr = strdup(expr);
                if (!ws->zones[izone].expr) {
                    rc = OSH_ENOMEM;
                    goto done;
                }
                expr[0] = '\0'; /* reset accumulator for next zone */
            }

            if (strcasecmp(key, OSH_GEO_KEY_END) == 0) {
                free(line);
                line = NULL;
                goto done; /* File now positioned just after the second END. */
            }

            /* New zone: advance the index (except for the very first zone). */
            if (!zone_active) {
                zone_active = 1;
                /* izone stays at 0 */
            } else {
                izone++;
            }

            if (izone >= ws->nzones) {
                OSH_DIAG_ERRORF(diag, "%s line %d: too many zones (max=%zu)", shf->filename, lineno, ws->nzones);
                rc = OSH_EPARSE;
                goto done;
            }

            ws->zones[izone].name = strdup(key);
            if (!ws->zones[izone].name) {
                rc = OSH_ENOMEM;
                goto done;
            }
        }

        /* Append args (if any) to the current zone's expression accumulator. */
        if (args) {
            rc = _str_append(&expr, args);
            if (rc != OSH_OK) {
                goto done;
            }
        }

        free(line);
        line = NULL;
    }

done:
    free(line);
    free(expr);
    return rc;
}

/* ---- Phase 3: material section ------------------------------------------- */

/**
 * @brief Handle one ASSIGNMAT / ASSIGNMA card.
 *
 * @details Assigns a material to a range of zones specified by name, with an
 * optional stride. Updates ws->zones[].material_name for each matching zone.
 *
 * @param[in,out] ws        Geometry workspace whose zones are updated.
 * @param[in,out] args      Argument string after the ASSIGNMAT keyword; tokenized in place.
 * @param[in]     filename  Source file name used in error messages.
 * @param[in]     lineno    Source line number used in error messages.
 *
 * @returns OSH_OK on success, OSH_EPARSE on invalid arguments, OSH_ENOMEM on failure.
 */
static enum osh_status _parse_assignmat(
    struct osh_geometry_workspace *ws, char *args, char const *filename, int lineno, struct osh_diag_sink const *diag) {
    char *cursor = args;
    char *mat_name;
    char *zname_start;
    char *zname_end;
    char *stride_str;
    size_t iz_start = 0u;
    size_t iz_end = 0u;
    size_t stride;
    size_t iz;
    enum osh_status rc;

    mat_name = _next_token(&cursor);
    if (!mat_name) {
        OSH_DIAG_ERRORF(diag, "%s line %d: ASSIGNMAT: missing material name", filename, lineno);
        return OSH_EPARSE;
    }

    zname_start = _next_token(&cursor);
    if (!zname_start) {
        OSH_DIAG_ERRORF(diag, "%s line %d: ASSIGNMAT: missing zone name", filename, lineno);
        return OSH_EPARSE;
    }
    if (!_zone_index_from_name(zname_start, ws->zones, ws->nzones, &iz_start)) {
        OSH_DIAG_ERRORF(diag, "%s line %d: ASSIGNMAT: unknown zone '%s'", filename, lineno, zname_start);
        return OSH_EPARSE;
    }

    zname_end = _next_token(&cursor);
    if (zname_end) {
        if (!_zone_index_from_name(zname_end, ws->zones, ws->nzones, &iz_end)) {
            OSH_DIAG_ERRORF(diag, "%s line %d: ASSIGNMAT: unknown zone '%s'", filename, lineno, zname_end);
            return OSH_EPARSE;
        }
    } else {
        iz_end = iz_start;
    }

    stride_str = _next_token(&cursor);
    stride = stride_str ? (size_t) atoi(stride_str) : 1u;
    if (stride == 0u) {
        OSH_DIAG_ERRORF(diag, "%s line %d: ASSIGNMAT: invalid stride 0", filename, lineno);
        return OSH_EPARSE;
    }

    if (iz_end < iz_start) {
        OSH_DIAG_ERRORF(diag, "%s line %d: ASSIGNMAT: zone range ends before it starts", filename, lineno);
        return OSH_EPARSE;
    }

    for (iz = iz_start; iz <= iz_end; iz += stride) {
        rc = _assign_material(&ws->zones[iz], mat_name, filename, lineno, diag);
        if (rc != OSH_OK) {
            return rc;
        }
    }
    return OSH_OK;
}

/**
 * @brief Parse the material section after the second END card.
 *
 * @details Handles both ASSIGNMAT cards and the legacy positional zone-list
 * format. Rewinds the file internally and skips past both END cards.
 *
 * @param[in,out] shf  Open geometry file; rewound internally.
 * @param[in,out] ws   Geometry workspace whose zones receive material names.
 *
 * @returns OSH_OK on success, OSH_EPARSE or OSH_ENOMEM on failure.
 */
static enum osh_status
_parse_media(struct oshfile *shf, struct osh_diag_sink const *diag, struct osh_geometry_workspace *ws) {
    char *line = NULL;
    char *key = NULL;
    char *args = NULL;
    int lineno;
    int nend = 0;
    int in_media = 0;
    int warned_legacy = 0;
    size_t izone = 0u;
    char *cursor;
    char *tok;
    enum osh_status rc = OSH_OK;

    if (!_rewind_oshfile(shf, diag)) {
        return OSH_EIO;
    }

    /* Skip past the second END card. */
    while (osh_readline_key(shf, &line, &key, &args, &lineno) > 0) {
        if (strcasecmp(key, OSH_GEO_KEY_END) == 0) {
            ++nend;
            if (nend == 2) {
                free(line);
                line = NULL;
                break;
            }
        }
        free(line);
        line = NULL;
    }

    if (nend < 2) {
        /* File has no material section — nothing to parse. */
        free(line);
        return OSH_OK;
    }

    /* Parse material assignments.
     *
     * Legacy positional format:
     *   - First nzones tokens form the zone-number list (ignored; zones are
     *     already named).
     *   - Next nzones tokens are material names assigned positionally.
     *
     * Preferred ASSIGNMAT format:
     *   ASSIGNMAT <material> <zone-start> [<zone-end> [<stride>]]
     */
    while (osh_readline_key(shf, &line, &key, &args, &lineno) > 0) {

        if ((strcasecmp(key, OSH_GEO_KEY_ASSIGNMAT) == 0) || (strcasecmp(key, OSH_GEO_KEY_ASSIGNMA) == 0)) {
            rc = _parse_assignmat(ws, args, shf->filename, lineno, diag);
            if (rc != OSH_OK) {
                free(line);
                return rc;
            }
            free(line);
            line = NULL;
            continue;
        }

        /* Legacy: count tokens and assign materials when in the media half. */
        izone++;
        if (izone > ws->nzones) {
            OSH_DIAG_ERRORF(
                diag, "%s line %d: too many entries in material section (max=%zu)", shf->filename, lineno, ws->nzones);
            free(line);
            return OSH_EPARSE;
        }

        if (in_media) {
            if (!warned_legacy) {
                OSH_DIAG_WARNF(diag, "Implicit geo.dat material lists are legacy; use ASSIGNMAT instead");
                warned_legacy = 1;
            }
            rc = _assign_material(&ws->zones[izone - 1u], key, shf->filename, lineno, diag);
            if (rc != OSH_OK) {
                free(line);
                return rc;
            }
        }

        cursor = args;
        tok = _next_token(&cursor);
        while (tok != NULL) {
            izone++;
            if (izone > ws->nzones) {
                OSH_DIAG_ERRORF(diag,
                                "%s line %d: too many entries in material section (max=%zu)",
                                shf->filename,
                                lineno,
                                ws->nzones);
                free(line);
                return OSH_EPARSE;
            }
            if (in_media) {
                rc = _assign_material(&ws->zones[izone - 1u], tok, shf->filename, lineno, diag);
                if (rc != OSH_OK) {
                    free(line);
                    return rc;
                }
            }
            tok = _next_token(&cursor);
        }

        /* After counting nzones entries we switch from the zone-number list
         * to the material-name list. */
        if (izone == ws->nzones) {
            in_media = 1;
            izone = 0u;
        }

        free(line);
        line = NULL;
    }

    free(line);
    return OSH_OK;
}

/* ---- Public entry point -------------------------------------------------- */

enum osh_status
osh_geometry_parse(struct oshfile *oshf, struct osh_diag_sink const *diag, struct osh_geometry_workspace *ws) {
    size_t nbodies;
    size_t nzones;
    enum osh_status rc;

    if (!oshf || !ws) {
        return OSH_EINVAL;
    }

    /* Count bodies and zones without retaining any state. */
    nbodies = _count_bodies(oshf, diag);
    nzones = _count_zones(oshf, diag);

    if (nbodies == 0u || nzones == 0u) {
        OSH_DIAG_ERRORF(diag, "geometry: '%s' has no bodies or no zones", oshf->filename);
        return OSH_EPARSE;
    }

    ws->nbodies = nbodies;
    ws->bodies = (struct osh_geometry_body *) calloc(nbodies, sizeof(struct osh_geometry_body));
    if (!ws->bodies) {
        return OSH_ENOMEM;
    }

    ws->nzones = nzones;
    ws->zones = (struct osh_geometry_zone *) calloc(nzones, sizeof(struct osh_geometry_zone));
    if (!ws->zones) {
        return OSH_ENOMEM;
    }

    /* Phase 1: parse bodies.
     * Rewinds internally; leaves file positioned after the first END card. */
    rc = _parse_bodies(oshf, diag, ws);
    if (rc != OSH_OK) {
        return rc;
    }

    /* Phase 2: parse zones.
     * Continues from where _parse_bodies stopped;
     * leaves file positioned after the second END card. */
    rc = _parse_zones(oshf, diag, ws);
    if (rc != OSH_OK) {
        return rc;
    }

    /* Phase 3: parse material assignments.
     * Rewinds internally; skips past both END cards. */
    rc = _parse_media(oshf, diag, ws);
    if (rc != OSH_OK) {
        return rc;
    }

    return OSH_OK;
}
