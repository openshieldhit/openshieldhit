#include "gemca/parse/osh_gemca2_parse_medium.h"

#include <ctype.h> /* isalpha() */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_logger.h"
#include "common/osh_readline.h"
#include "gemca/osh_gemca2.h"
#include "gemca/osh_gemca2_defines.h"
#include "gemca/parse/osh_gemca2_parse_keys.h"

int _assign_material(struct gemca_workspace *g, char *args, int lineno);
int _get_zoneid_from_name(char const *zname, struct gemca_workspace const *g);

/**
 * @brief Parse zone information
 *
 * @details This function parses the material part of the geo.dat file.
 * (1st part is body description, 2nd is zone description, 3rd is material description)
 *
 * @param[in] fp - file pointer to file open for reading.
 * @param[in,out] g - gemca workspace pointer
 *
 * @returns 1 if format is OK, 0 otherwise
 *
 * @author Niels Bassler
 */
int osh_gemca_parse_media(struct oshfile *shf, struct gemca_workspace *g) {
    char *key = NULL;
    char *args = NULL;
    char *line = NULL;
    char *arg = NULL;

    /* The material assignment consists of two sets of data.
     *   - First set is a list of body numbers. This will only be counted, but not parsed.
     *   - Second set contains the media data
     * These two sets are matched one by one.
     */
    int in_media = 0; /* flag whether we are in the media set */
    size_t izone = 0;

    int lineno;

    /* move to the second end statement */
    rewind(shf->fp);

    while (osh_readline_key(shf, &line, &key, &args, &lineno) > 0) {
        if (strcasecmp(OSH_GEMCA_KEY_END, key) == 0) {
            in_media++;
            if (in_media == 2) /* jump out after second END statement was read */
                break;
        }
        free(line);
    }
    free(line);

    /* we know how many zones there are, so we will read exactly this number of zones into the medium list. */
    /* readline is maybe not so optimal, since there is no key */
    while (osh_readline_key(shf, &line, &key, &args, &lineno) > 0) {

        /* optionally, materials can also be assigned to zones by the (paritally) FLUKA compatible ASSIGNMAT key */
        if ((strcasecmp(OSH_GEMCA_KEY_ASSIGNMAT, key) == 0) || (strcasecmp(OSH_GEMCA_KEY_ASSIGNMA, key) == 0)) {
            _assign_material(g, args, lineno);
            free(line);
            continue; /* next line */
        }

        izone++; /* first zone in the line just read, will always be in the key, so count it. */
        /* if we are in the media block, assign the value */
        if (in_media) {
            g->zones[izone - 1]->medium = atoi(key);
        }

        /* next run trough the remaining args on the line */

        arg = strtok(args, " ");
        while (arg != NULL) {
            izone++;

            if (izone > g->nzones) {
                fprintf(stderr,
                        "    Found zones %llu, but expected %llu\n",
                        (long long unsigned int) izone,
                        (long long unsigned int) g->nzones);
                osh_error(EX_CONFIG, "Too many zones found %s line number %i\n", g->filename, lineno);
            }
            /* again, if we are in the media block, assign it */
            if (in_media) {
                g->zones[izone - 1]->medium = atoi(arg);
            }
            arg = strtok(NULL, " "); /* next arguments */
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
    return 1;
}

/**
 * @brief Assign material to zones based on ASSIGNMA(T) key
 *
 * @param[in,out] g Gemca workspace pointer
 * @param[in] args Arguments string containing material and zone information
 * @param[in] lineno Line number in the input file for error reporting
 *
 * @return int Returns 1 on success
 *
 * @author Niels Bassler
 */
int _assign_material(struct gemca_workspace *g, char *args, int lineno) {

    char *arg;

    size_t zone_start;
    size_t zone_end;
    size_t stride;
    size_t iz;

    int medium = 0;

    /* first argument is the material index (names will be supported later)*/
    arg = strtok(args, " ");
    if (arg != NULL) {
        medium = atoi(arg);
    } else {
        osh_error(EX_CONFIG, "No medium index found in ASSIGNMA(T) %s line number %i\n", g->filename, lineno);
    }

    /* next two arguments are the zone indices, or a range of zone indices */
    arg = strtok(NULL, " ");
    if (arg == NULL) {
        osh_error(EX_CONFIG, "No zone index found in ASSIGNMA(T) %s line number %i\n", g->filename, lineno);
    }
    /* check if zone name was given, if so, get the index for that zone name */
    iz = _get_zoneid_from_name(arg, g);
    if (iz > 0 && iz < g->nzones) {
        /* zone name found, set both start and end to this index +1 (since zones are 1-based) */
        zone_start = iz;
    } else {
        /* not a zone name, assume its a zone index or range */
        zone_start = atoi(arg);
    }

    arg = strtok(NULL, " ");
    if (arg != NULL) {
        /* check if zone name was given, if so, get the index for that zone name */
        iz = _get_zoneid_from_name(arg, g);
        if (iz > 0 && iz < g->nzones) {
            zone_end = iz;
        } else {
            zone_end = atoi(arg);
        }
    } else {
        /* only a single zone index given, set end to start */
        zone_end = zone_start;
    }

    /* last supported argument is stride of the range */
    arg = strtok(NULL, " ");
    stride = 1; /* default stride */
    if (arg != NULL) {
        stride = atoi(arg);
        if (stride == 0) {
            osh_error(EX_CONFIG, "Invalid stride 0 in ASSIGNMA(T) %s line number %i\n", g->filename, lineno);
        }
    }

    /* assign the medium to the zones in the specified range */
    for (iz = zone_start; iz <= zone_end; iz += stride) {
        if (iz == 0 || iz > g->nzones) {
            osh_error(EX_CONFIG,
                      "Zone index %llu out of range in ASSIGNMA(T) %s line number %i\n",
                      (long long unsigned int) iz,
                      g->filename,
                      lineno);
        }
        g->zones[iz - 1]->medium = medium;
        printf("    Assigned medium %i to zoneID %llu named '%s'\n",
               medium,
               (long long unsigned int) iz,
               g->zones[iz - 1]->name);
    }
    return 1;
}

int _get_zoneid_from_name(char const *zname, struct gemca_workspace const *g) {
    size_t iz;

    for (iz = 0; iz < g->nzones; iz++) {
        if (strcmp(zname, g->zones[iz]->name) == 0) {
            return g->zones[iz]->id;
        }
    }
    return 0; /* not found */
}