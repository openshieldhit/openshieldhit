#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_logger.h"
#include "common/osh_readline.h"
#include "gemca/osh_gemca2.h"
#include "gemca/parse/osh_gemca2_parse_body.h"
#include "gemca/parse/osh_gemca2_parse_keys.h"
#include "gemca/parse/osh_gemca2_parse_medium.h"
#include "gemca/parse/osh_gemca2_parse_zone.h"

static int _test_format(struct oshfile *shf, size_t *nbody, size_t *nzone);
static int _rewind_oshfile(struct oshfile *shf);

/**
 * @brief loads and parses the geometry geo.dat file (or whatever filename was specified.

 * @param[in] filename - relatove path to filename to be opened (or absolute)s
 *
 * @param[in,out] *g - pointer to gemca workspace, may be initialized with NULL, when invoking this function.
 *
 * @returns 1
 *
 * @author Niels Bassler
 */
int osh_gemca_parse(char const *filename, struct gemca_workspace *g) {

    size_t i;
    size_t nbody; /* number of bodies */
    size_t nzone; /* number of zones */

    struct body **body = NULL;
    struct zone **zone = NULL;

    struct oshfile *shf = osh_fopen(filename);
    if (!shf) {
        return 0;
    }

    /* test_format() also counts the number of bodies and zones */
    if (!_test_format(shf, &nbody, &nzone)) {
        osh_error("Unknown format of %s\n", filename);
        osh_fclose(shf);
        return 0;
    }

    g->filename = calloc(strlen(filename) + 1, sizeof(char));
    if (!g->filename) {
        osh_alloc_failed("g->filename");
        osh_fclose(shf);
        return 0;
    }
    snprintf(g->filename, strlen(filename) + 1, "%s", filename);

    /* allocate memory for pointers to lists */
    body = (struct body **) calloc(nbody, sizeof(struct body *));
    if (body == NULL) {
        osh_alloc_failed("*body");
        osh_fclose(shf);
        return 0;
    }
    zone = (struct zone **) calloc(nzone, sizeof(struct zone *));
    if (zone == NULL) {
        osh_alloc_failed("*zone");
        free(body);
        osh_fclose(shf);
        return 0;
    }

    /* allocate memory for every list item */
    for (i = 0; i < nbody; i++) {
        if (!osh_gemca_body_init(&body[i])) {
            osh_fclose(shf);
            return 0;
        }
    }
    for (i = 0; i < nzone; i++) {
        if (!osh_gemca_zone_init(&zone[i])) {
            osh_fclose(shf);
            return 0;
        }
    }

    /* initialize gemca workspace */
    g->bodies = body;
    g->zones = zone;
    g->nbodies = nbody;
    g->nzones = nzone;

    if (!osh_gemca_parse_bodies(shf, g)) {
        osh_fclose(shf);
        return 0;
    }
    if (!osh_gemca_parse_zones(shf, g)) {
        osh_fclose(shf);
        return 0;
    }
    if (!osh_gemca_parse_media(shf, g)) {
        osh_fclose(shf);
        return 0;
    }

    osh_fclose(shf);

    return 1;
}

/**
 * @brief Checks if file looks like a proper geo.dat file and count the number of bodies and zones.
 *
 * @param[in] fp - file pointer to file open for reading.
 * @param[out] nbody - number of bodies found
 * @param[out] nzone - number of zones found
 *
 * @returns 1 if format is OK, 0 otherwise
 *
 * @author Niels Bassler
 */
static int _test_format(struct oshfile *shf, size_t *nbody, size_t *nzone) {

    int ret = 0;

    if (!_rewind_oshfile(shf)) {
        return 0;
    }

    *nbody = osh_gemca_parse_count_bodies(shf);
    *nzone = osh_gemca_parse_count_zones(shf);

    if ((*nzone > 1) || (*nbody > 1)) {
        ret = 1;
    }
    // printf("GEMCA test_format() : Found %lu bodies and %lu zones. %i\n", *nbody, *nzone, ret);
    return ret;
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
        osh_error("Failed to rewind geometry file '%s'\n", shf->filename);
        return 0;
    }
    shf->lineno = 0;
    return 1;
}
