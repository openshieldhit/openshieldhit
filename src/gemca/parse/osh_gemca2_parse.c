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
int osh_gemca_parse(const char *filename, struct gemca_workspace *g) {

    size_t i;
    size_t nbody; /* number of bodies */
    size_t nzone; /* number of zones */

    struct body **body = NULL;
    struct zone **zone = NULL;

    struct oshfile *shf;

    shf = calloc(1, sizeof(struct oshfile));
    if (shf == NULL) {
        osh_error(EX_UNAVAILABLE, "could not allocate memory.\n");
    }

    shf = osh_fopen(filename);

    /* test_format() also counts the number of bodies and zones */
    if (!_test_format(shf, &nbody, &nzone)) {
        osh_error(EX_CONFIG, "Unknown format of %s\n", filename);
    }

    g->filename = calloc(strlen(filename) + 1, sizeof(char));
    snprintf(g->filename, strlen(filename) + 1, "%s", filename);

    /* allocate memory for pointers to lists */
    body = calloc(nbody, sizeof(struct body *));
    if (body == NULL) {
        fprintf(stderr, "*** Error in memory allocation, *body.\n");
        exit(EX_UNAVAILABLE);
    }
    zone = calloc(nzone, sizeof(struct zone *));
    if (zone == NULL) {
        fprintf(stderr, "*** Error in memory allocation, *zone.\n");
        exit(EX_UNAVAILABLE);
    }

    /* allocate memory for every list item */
    for (i = 0; i < nbody; i++) {
        osh_gemca_body_init(&body[i]);
    }
    for (i = 0; i < nzone; i++) {
        osh_gemca_zone_init(&zone[i]);
    }

    /* initialize gemca workspace */
    g->bodies = body;
    g->zones = zone;
    g->nbodies = nbody;
    g->nzones = nzone;

    osh_gemca_parse_bodies(shf, g);
    osh_gemca_parse_zones(shf, g);
    osh_gemca_parse_media(shf, g);

    osh_fclose(shf);
    free(shf);

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

    rewind(shf->fp);

    *nbody = osh_gemca_parse_count_bodies(shf);
    *nzone = osh_gemca_parse_count_zones(shf);

    if ((*nzone > 1) || (*nbody > 1)) {
        ret = 1;
    }
    // printf("GEMCA test_format() : Found %lu bodies and %lu zones. %i\n", *nbody, *nzone, ret);
    return ret;
}
