#include "gemca/parse/osh_gemca2_parse_medium.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_logger.h"
#include "common/osh_readline.h"
#include "gemca/osh_gemca2.h"
#include "gemca/osh_gemca2_defines.h"
#include "gemca/parse/osh_gemca2_parse_keys.h"

static enum osh_status _assign_material(struct gemca_workspace *g, char *args, int lineno);
static enum osh_status
_assign_material_name_to_zone(struct gemca_workspace *g, struct zone *z, char const *raw_name, int lineno);
static int _get_zone_index_from_name(char const *zname, struct gemca_workspace const *g, size_t *index_out);
static char const *_normalize_material_name(char const *raw_name, struct gemca_workspace const *g, int lineno);
static char *_next_token(char **cursor);
static int _rewind_oshfile(struct oshfile *shf);

static char const GEMCA_MATERIAL_NAME_BLACKHOLE[] = "blackhole";
static char const GEMCA_MATERIAL_NAME_VACUUM[] = "vacuum";

/**
 * @brief Parse zone information
 *
 * @details This function parses the material part of the geo.dat file.
 * (1st part is body description, 2nd is zone description, 3rd is material description)
 * Material names and zone names are stored as strings here. The parser does not look up material definitions and does
 * not populate zone::material_idx; that cross-reference is resolved later by the pre-simulation assembly layer after
 * both geometry and material workspaces exist.
 *
 * @param[in] fp - file pointer to file open for reading.
 * @param[in,out] g - gemca workspace pointer
 *
 * @returns OSH_OK on success, OSH_E* on failure.
 *
 * @author Niels Bassler
 */
enum osh_status osh_gemca_parse_media(struct oshfile *shf, struct gemca_workspace *g) {
    char *key = NULL;
    char *args = NULL;
    char *line = NULL;
    char *arg = NULL;
    char *cursor = NULL;
    enum osh_status rc;

    /* The legacy implicit material assignment consists of two sets of data.
     *   - First set is a list of body numbers. This will only be counted, but not parsed.
     *   - Second set contains the material names
     * These two sets are matched one by one.
     * This legacy format is discouraged because it relies on positional zone ordering. Prefer ASSIGNMAT with explicit
     * material and zone names.
     */
    int in_media = 0; /* flag whether we are in the media set */
    int warned_legacy_media_list = 0;
    size_t izone = 0;

    int lineno;

    /* move to the second end statement */
    if (!_rewind_oshfile(shf)) {
        return OSH_EIO;
    }

    while (osh_readline_key(shf, &line, &key, &args, &lineno) > 0) {
        if (strcasecmp(OSH_GEMCA_KEY_END, key) == 0) {
            in_media++;
            if (in_media == 2) { /* jump out after second END statement was read */
                break;
            }
        }
        free(line);
    }
    free(line);
    in_media = 0;

    /* we know how many zones there are, so we will read exactly this number of zones into the medium list. */
    /* readline is maybe not so optimal, since there is no key */
    while (osh_readline_key(shf, &line, &key, &args, &lineno) > 0) {

        /* optionally, materials can also be assigned to zones by the (paritally) FLUKA compatible ASSIGNMAT key */
        if ((strcasecmp(OSH_GEMCA_KEY_ASSIGNMAT, key) == 0) || (strcasecmp(OSH_GEMCA_KEY_ASSIGNMA, key) == 0)) {
            rc = _assign_material(g, args, lineno);
            if (rc != OSH_OK) {
                free(line);
                return rc;
            }
            free(line);
            continue; /* next line */
        }

        izone++; /* first zone in the line just read, will always be in the key, so count it. */
        if (izone > g->nzones) {
            osh_error("Too many zones found: %llu (expected %llu) in %s line %i",
                      (long long unsigned int) izone,
                      (long long unsigned int) g->nzones,
                      g->filename,
                      lineno);
            free(line);
            return OSH_EPARSE;
        }
        /* if we are in the media block, assign the value */
        if (in_media) {
            if (!warned_legacy_media_list) {
                osh_warn("Implicit geo.dat material lists are legacy; use ASSIGNMAT <material-name> <zone> instead");
                warned_legacy_media_list = 1;
            }
            rc = _assign_material_name_to_zone(g, g->zones[izone - 1], key, lineno);
            if (rc != OSH_OK) {
                free(line);
                return rc;
            }
        }

        /* next run trough the remaining args on the line */

        cursor = args;
        arg = _next_token(&cursor);
        while (arg != NULL) {
            izone++;

            if (izone > g->nzones) {
                osh_error("Too many zones found: %llu (expected %llu) in %s line %i",
                          (long long unsigned int) izone,
                          (long long unsigned int) g->nzones,
                          g->filename,
                          lineno);
                free(line);
                return OSH_EPARSE;
            }
            /* again, if we are in the media block, assign it */
            if (in_media) {
                rc = _assign_material_name_to_zone(g, g->zones[izone - 1], arg, lineno);
                if (rc != OSH_OK) {
                    free(line);
                    return rc;
                }
            }
            arg = _next_token(&cursor); /* next arguments */
        }

        /* check if we have read all the expected zones */
        if (izone == g->nzones) {
            /* we switched to from the zone-block to the media block */
            in_media = 1;
            izone = 0;
        }
        free(line);
    } /* end of while loop */
    free(line);
    return OSH_OK;
}

/**
 * @brief Assign material to zones based on ASSIGNMA(T) key
 *
 * @param[in,out] g Gemca workspace pointer
 * @param[in] args Arguments string containing material and zone information
 * @param[in] lineno Line number in the input file for error reporting
 *
 * @returns OSH_OK on success, OSH_E* on failure.
 *
 * @author Niels Bassler
 */
static enum osh_status _assign_material(struct gemca_workspace *g, char *args, int lineno) {

    char *arg;
    char *cursor;
    enum osh_status rc;

    size_t zone_start_index;
    size_t zone_end_index;
    size_t stride;
    size_t iz;

    char const *material_name;

    if (args == NULL) {
        osh_error("No arguments found in ASSIGNMA(T) %s line number %i", g->filename, lineno);
        return OSH_EPARSE;
    }

    /* first argument is the material name */
    cursor = args;
    arg = _next_token(&cursor);
    if (arg != NULL) {
        material_name = arg;
    } else {
        osh_error("No material name found in ASSIGNMA(T) %s line number %i", g->filename, lineno);
        return OSH_EPARSE;
    }

    /* next two arguments are exact zone names, or a range delimited by two exact zone names.
     * There is intentionally no numeric fallback: "001" and "1" are different zone names. */
    arg = _next_token(&cursor);
    if (arg == NULL) {
        osh_error("No zone name found in ASSIGNMA(T) %s line number %i", g->filename, lineno);
        return OSH_EPARSE;
    }
    if (!_get_zone_index_from_name(arg, g, &zone_start_index)) {
        osh_error("Unknown zone name '%s' in ASSIGNMA(T) %s line number %i", arg, g->filename, lineno);
        return OSH_EPARSE;
    }

    arg = _next_token(&cursor);
    if (arg != NULL) {
        if (!_get_zone_index_from_name(arg, g, &zone_end_index)) {
            osh_error("Unknown zone name '%s' in ASSIGNMA(T) %s line number %i", arg, g->filename, lineno);
            return OSH_EPARSE;
        }
    } else {
        zone_end_index = zone_start_index;
    }

    /* last supported argument is stride of the range */
    arg = _next_token(&cursor);
    stride = 1; /* default stride */
    if (arg != NULL) {
        stride = atoi(arg);
        if (stride == 0) {
            osh_error("Invalid stride 0 in ASSIGNMA(T) %s line number %i", g->filename, lineno);
            return OSH_EPARSE;
        }
    }

    if (zone_end_index < zone_start_index) {
        osh_error("Zone range '%s' to '%s' ends before it starts in ASSIGNMA(T) %s line number %i",
                  g->zones[zone_start_index]->name,
                  g->zones[zone_end_index]->name,
                  g->filename,
                  lineno);
        return OSH_EPARSE;
    }

    /* assign the material to the zones in the specified range */
    for (iz = zone_start_index; iz <= zone_end_index; iz += stride) {
        rc = _assign_material_name_to_zone(g, g->zones[iz], material_name, lineno);
        if (rc != OSH_OK) {
            return rc;
        }
        osh_debug("    Assigned material '%s' to zone index %llu named '%s'",
                  g->zones[iz]->material_name,
                  (long long unsigned int) iz,
                  g->zones[iz]->name);
    }
    return OSH_OK;
}

/**
 * @brief Assign a material name to a zone.
 *
 * @details Normalizes only the explicitly supported legacy material names
 * "0" and "1000". All other material names are stored exactly as strings and
 * are resolved against mat.dat later by the pre-simulation assembly layer.
 *
 * @param[in] g         GEMCA workspace for diagnostics.
 * @param[in,out] z     Zone to update.
 * @param[in] raw_name  Material name token from geo.dat.
 * @param[in] lineno    geo.dat line number for diagnostics.
 *
 * @returns OSH_OK on success, OSH_E* on failure.
 */
static enum osh_status
_assign_material_name_to_zone(struct gemca_workspace *g, struct zone *z, char const *raw_name, int lineno) {
    char const *material_name;
    char *copy;
    size_t len;

    if (!z || !raw_name) {
        return OSH_EINVAL;
    }

    material_name = _normalize_material_name(raw_name, g, lineno);
    len = strlen(material_name);
    copy = (char *) malloc(len + 1u);
    if (!copy) {
        return OSH_ENOMEM;
    }
    memcpy(copy, material_name, len + 1u);

    free(z->material_name);
    z->material_name = copy;

    return OSH_OK;
}

/**
 * @brief Normalize supported legacy material-name tokens.
 *
 * @details "0" maps to "blackhole" and "1000" maps to "vacuum" for
 * compatibility with old geo.dat files. Numeric-looking names otherwise remain
 * ordinary strings.
 *
 * @param[in] raw_name  Material name token from geo.dat.
 * @param[in] g         GEMCA workspace for diagnostics.
 * @param[in] lineno    geo.dat line number for diagnostics.
 *
 * @returns Normalized material name pointer. The return value is either a
 *          static string or @p raw_name; callers must copy it if they need
 *          ownership.
 */
static char const *_normalize_material_name(char const *raw_name, struct gemca_workspace const *g, int lineno) {
    if (strcmp(raw_name, "0") == 0) {
        osh_warn("%s line %i: legacy material name '0' maps to '%s'; use '%s' explicitly",
                 g->filename,
                 lineno,
                 GEMCA_MATERIAL_NAME_BLACKHOLE,
                 GEMCA_MATERIAL_NAME_BLACKHOLE);
        return GEMCA_MATERIAL_NAME_BLACKHOLE;
    }
    if (strcmp(raw_name, "1000") == 0) {
        osh_warn("%s line %i: legacy material name '1000' maps to '%s'; use '%s' explicitly",
                 g->filename,
                 lineno,
                 GEMCA_MATERIAL_NAME_VACUUM,
                 GEMCA_MATERIAL_NAME_VACUUM);
        return GEMCA_MATERIAL_NAME_VACUUM;
    }

    return raw_name;
}

/**
 * @brief Find a zone by exact user-facing name.
 *
 * @details Zone names are strings; "001" and "1" are intentionally different.
 *
 * @param[in]  zname      Zone name to find.
 * @param[in]  g          GEMCA workspace to search.
 * @param[out] index_out  Receives the dense internal zone index on success.
 *
 * @returns 1 if found, 0 if no matching zone exists.
 */
static int _get_zone_index_from_name(char const *zname, struct gemca_workspace const *g, size_t *index_out) {
    size_t iz;

    for (iz = 0; iz < g->nzones; iz++) {
        if (strcmp(zname, g->zones[iz]->name) == 0) {
            *index_out = iz;
            return 1;
        }
    }
    return 0; /* not found */
}

/**
 * @brief Return the next whitespace-delimited token from a mutable string.
 *
 * @details This helper intentionally replaces `strtok()`: MSVC warns about
 * `strtok()` as unsafe, while `strtok_s()` is not portable across the C
 * libraries we support. The parser only needs simple in-place whitespace
 * tokenization, so a local helper is clearer than platform-specific wrappers.
 *
 * @param[in,out] cursor  On input, current parse position. On output, position
 *                        after the returned token.
 *
 * @returns Pointer to the next token inside the input buffer, or NULL when no
 *          token remains. The input buffer is modified by inserting NUL
 *          terminators.
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
 * @brief Rewind an `oshfile` stream and reset its tracked line number.
 *
 * @param[in,out] shf Open geometry file wrapper to rewind.
 *
 * @returns Nothing. Exits via `osh_error()` if rewinding fails.
 */
static int _rewind_oshfile(struct oshfile *shf) {
    if (fseek(shf->fp, 0L, SEEK_SET) != 0) {
        osh_error("Failed to rewind geometry file '%s'", shf->filename);
        return 0;
    }
    shf->lineno = 0;
    return 1;
}
