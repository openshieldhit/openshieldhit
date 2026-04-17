#include "apps/osh/osh_geometry_parse.h"

#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_geometry_parse_keys.h"
#include "gemca/osh_gemca2_defines.h"
#include "openshieldhit/geometry.h"
#include "openshieldhit/geometry_defs.h"
#include "openshieldhit/logger.h"
#include "openshieldhit/readline.h"
#include "openshieldhit/status.h"

/* ---- Internal helpers ---------------------------------------------------- */

static int _rewind_oshfile(struct oshfile *shf);

/* Map a body-type key string to an OSH_GEOMETRY_BODY_* code.
 * Values are identical to the internal OSH_GEMCA_BODY_* constants so that
 * osh_geometry_workspace_prepare() can use them without conversion. */
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

/* Return non-zero if key begins a zone continuation line (operator token). */
static int _is_zone_continuation(char const *key) {
    return key && (key[0] == '+' || key[0] == '-' || key[0] == '|' || key[0] == '(' || key[0] == ')');
}

/* Count body cards in the geometry file without retaining state.
 * This stays app-side because it is part of geo.dat parsing, not geometry
 * preparation. */
static size_t _count_bodies(struct oshfile *shf) {
    int lineno;
    char *line = NULL;
    char *key = NULL;
    char *args = NULL;
    size_t nbody = 0u;

    if (!_rewind_oshfile(shf)) {
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

/* Count zones in the geometry file without building any cold objects.
 * Like _count_bodies(), this belongs to the app parser because it depends on
 * geo.dat section structure. */
static size_t _count_zones(struct oshfile *shf) {
    int lineno;
    char *line = NULL;
    char *key = NULL;
    char *args = NULL;
    size_t nzone = 0u;
    int in_block = 0;

    if (!_rewind_oshfile(shf)) {
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
    osh_info("Found %zu zones in geo.dat file", nzone);
    return nzone;
}

/* Append NUL-terminated string b to *a, growing the allocation as needed.
 * *a must point to a heap-allocated, NUL-terminated string. */
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

/* Find a zone by exact name in the flat cold-zone array.
 * Returns 1 on success (fills *idx_out), 0 if not found. */
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

/* Return the next whitespace-delimited token from *cursor, NULL when done.
 * Modifies the buffer in place by inserting NUL bytes. */
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

/* Normalize legacy numeric material names.
 * "0" → "blackhole", "1000" → "vacuum", all others pass through unchanged. */
static char const *_normalize_material_name(char const *raw, char const *filename, int lineno) {
    if (strcmp(raw, "0") == 0) {
        osh_warn("%s line %d: legacy material '0' mapped to 'blackhole'; use 'blackhole' explicitly", filename, lineno);
        return "blackhole";
    }
    if (strcmp(raw, "1000") == 0) {
        osh_warn("%s line %d: legacy material '1000' mapped to 'vacuum'; use 'vacuum' explicitly", filename, lineno);
        return "vacuum";
    }
    return raw;
}

/* Duplicate a normalized material name into z->material_name. */
static enum osh_status
_assign_material(struct osh_geometry_zone *z, char const *raw_name, char const *filename, int lineno) {
    char const *name = _normalize_material_name(raw_name, filename, lineno);
    char *copy = strdup(name);
    if (!copy) {
        return OSH_ENOMEM;
    }
    free(z->material_name);
    z->material_name = copy;
    return OSH_OK;
}

static int _rewind_oshfile(struct oshfile *shf) {
    if (fseek(shf->fp, 0L, SEEK_SET) != 0) {
        osh_error("Failed to rewind geometry file '%s'", shf->filename);
        return 0;
    }
    shf->lineno = 0;
    return 1;
}

/* ---- Phase 1: body section ----------------------------------------------- */

/* Read the body section (from file start up to and including the first END).
 * Fills ws->bodies[0..ws->nbodies-1] with type, name, and raw argument arrays.
 * Leaves the file position just after the first END card. */
static enum osh_status _parse_bodies(struct oshfile *shf, struct osh_geometry_workspace *ws) {
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

    if (!_rewind_oshfile(shf)) {
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
                osh_error("%s line %d: too many bodies (max=%zu)", shf->filename, lineno, ws->nbodies);
                rc = OSH_EPARSE;
                goto done;
            }

            btype = btype_new;
            body_active = 1;

            if (!args) {
                osh_error("%s line %d: missing body name or parameters", shf->filename, lineno);
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
                osh_error("%s line %d: too many body arguments", shf->filename, lineno);
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

/* Read the zone section (from just after the first END to the second END).
 * Fills ws->zones[].name and ws->zones[].expr.
 * Precondition: file is positioned just after the first END (left by _parse_bodies). */
static enum osh_status _parse_zones(struct oshfile *shf, struct osh_geometry_workspace *ws) {
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
                osh_error("%s line %d: zone continuation before any zone header", shf->filename, lineno);
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
                    osh_error("%s line %d: too many zones (max=%zu)", shf->filename, lineno, ws->nzones);
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
                osh_error("%s line %d: too many zones (max=%zu)", shf->filename, lineno, ws->nzones);
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

/* Handle one ASSIGNMAT / ASSIGNMA card: assign material to a zone range. */
static enum osh_status
_parse_assignmat(struct osh_geometry_workspace *ws, char *args, char const *filename, int lineno) {
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
        osh_error("%s line %d: ASSIGNMAT: missing material name", filename, lineno);
        return OSH_EPARSE;
    }

    zname_start = _next_token(&cursor);
    if (!zname_start) {
        osh_error("%s line %d: ASSIGNMAT: missing zone name", filename, lineno);
        return OSH_EPARSE;
    }
    if (!_zone_index_from_name(zname_start, ws->zones, ws->nzones, &iz_start)) {
        osh_error("%s line %d: ASSIGNMAT: unknown zone '%s'", filename, lineno, zname_start);
        return OSH_EPARSE;
    }

    zname_end = _next_token(&cursor);
    if (zname_end) {
        if (!_zone_index_from_name(zname_end, ws->zones, ws->nzones, &iz_end)) {
            osh_error("%s line %d: ASSIGNMAT: unknown zone '%s'", filename, lineno, zname_end);
            return OSH_EPARSE;
        }
    } else {
        iz_end = iz_start;
    }

    stride_str = _next_token(&cursor);
    stride = stride_str ? (size_t) atoi(stride_str) : 1u;
    if (stride == 0u) {
        osh_error("%s line %d: ASSIGNMAT: invalid stride 0", filename, lineno);
        return OSH_EPARSE;
    }

    if (iz_end < iz_start) {
        osh_error("%s line %d: ASSIGNMAT: zone range ends before it starts", filename, lineno);
        return OSH_EPARSE;
    }

    for (iz = iz_start; iz <= iz_end; iz += stride) {
        rc = _assign_material(&ws->zones[iz], mat_name, filename, lineno);
        if (rc != OSH_OK) {
            return rc;
        }
    }
    return OSH_OK;
}

/* Read the material section (after the second END).
 * Handles both ASSIGNMAT cards and the legacy positional zone-list format.
 * Rewinds the file and skips past both END cards before reading. */
static enum osh_status _parse_media(struct oshfile *shf, struct osh_geometry_workspace *ws) {
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

    if (!_rewind_oshfile(shf)) {
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
            rc = _parse_assignmat(ws, args, shf->filename, lineno);
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
            osh_error("%s line %d: too many entries in material section (max=%zu)", shf->filename, lineno, ws->nzones);
            free(line);
            return OSH_EPARSE;
        }

        if (in_media) {
            if (!warned_legacy) {
                osh_warn("Implicit geo.dat material lists are legacy; use ASSIGNMAT instead");
                warned_legacy = 1;
            }
            rc = _assign_material(&ws->zones[izone - 1u], key, shf->filename, lineno);
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
                osh_error(
                    "%s line %d: too many entries in material section (max=%zu)", shf->filename, lineno, ws->nzones);
                free(line);
                return OSH_EPARSE;
            }
            if (in_media) {
                rc = _assign_material(&ws->zones[izone - 1u], tok, shf->filename, lineno);
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

enum osh_status osh_geometry_parse(struct oshfile *oshf, struct osh_geometry_workspace *ws) {
    size_t nbodies;
    size_t nzones;
    enum osh_status rc;

    if (!oshf || !ws) {
        return OSH_EINVAL;
    }

    /* Count bodies and zones without retaining any state. */
    nbodies = _count_bodies(oshf);
    nzones = _count_zones(oshf);

    if (nbodies == 0u || nzones == 0u) {
        osh_error("geometry: '%s' has no bodies or no zones", oshf->filename);
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
    rc = _parse_bodies(oshf, ws);
    if (rc != OSH_OK) {
        return rc;
    }

    /* Phase 2: parse zones.
     * Continues from where _parse_bodies stopped;
     * leaves file positioned after the second END card. */
    rc = _parse_zones(oshf, ws);
    if (rc != OSH_OK) {
        return rc;
    }

    /* Phase 3: parse material assignments.
     * Rewinds internally; skips past both END cards. */
    rc = _parse_media(oshf, ws);
    if (rc != OSH_OK) {
        return rc;
    }

    return OSH_OK;
}
