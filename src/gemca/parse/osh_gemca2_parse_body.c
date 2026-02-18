#include "gemca/parse/osh_gemca2_parse_body.h"

#include <ctype.h> /* isalpha() */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_logger.h"
#include "common/osh_readline.h"
#include "gemca/osh_gemca2.h"
#include "gemca/osh_gemca2_defines.h"
#include "gemca/parse/osh_gemca2_parse_keys.h"

static int _save_body(struct body *b, char *nstr, double *par, int npar, int btype);
static int _body_from_key(const char *key);

/**
 * @brief Initialize a body.
 *
 * @param[in,out] **body - body to be initialized, memory will be allocated and zeroed.
 *
 * @returns 1
 *
 * @author Niels Bassler
 */
int osh_gemca_body_init(struct body **body) {

    *body = calloc(1, sizeof(struct body));
    if (*body == NULL) {
        osh_error(EX_UNAVAILABLE, "osh_gemca_body_init() cannot allocate memory\n");
    }
    return 1;
}

/**
 * @brief  Return number of bodies found in the given geo.dat file.
 *
 * @details This function only parses the keys (the first word on each valid line)
 *          File pointer will be rewinded first by this funtion, but not when done.
 *
 * @param[in] shf - a struct oshfile with an open file pointer.
 *
 * @returns Number of bodies found.
 *
 * @author Niels Bassler
 */
size_t osh_gemca_parse_count_bodies(struct oshfile *shf) {
    int lineno;
    char *key = NULL;
    char *args = NULL;
    char *line = NULL;

    int nbody = 0;

    rewind(shf->fp);

    while (osh_readline_key(shf, &line, &key, &args, &lineno) > 0) {

        if (strcasecmp(OSH_GEMCA_KEY_END, key) != 0) {
            if (_body_from_key(key)) {
                ++nbody;
            }
        }
        free(line);
    }
    free(line);
    return nbody;
}

/**
 * @brief Reads from the geo.dat file and populates the g->bodies list.
 *
 * @details
 *
 * @param[in] *shf - pointer to oshfile which also keeps track of line numbers.
 * @param[out] *g - workspace struct, where g->bodies list will be made.
 *
 * @returns 1 // TODO could also be number of bodies found.
 *
 * @author Niels Bassler
 */
int osh_gemca_parse_bodies(struct oshfile *shf, struct gemca_workspace *g) {
    char *key = NULL;
    char *args = NULL;
    char *line = NULL;
    int lineno;
    int lineno_b; /* valid only when current_body != NULL */

    int nt;
    int btype = OSH_GEMCA_BODY_NONE;
    int btype_new;

    int nend = 0; /* keep if you plan to use it; otherwise remove */

    double par[OSH_GEMCA_NARGS_MAX];
    int npar = 0;
    int off = 0;

    char nstr[OSH_GEMCA_BODY_NAME_MAXLEN];

    struct body *current_body = NULL;

    size_t ibody = 0;
    size_t _ib;

    rewind(shf->fp);

    /* read line by line and parse the keys and arguments */
    while (osh_readline_key(shf, &line, &key, &args, &lineno) > 0) {

        /* END check early so we do not touch parsing state on END lines */
        if ((strcasecmp(key, OSH_GEMCA_KEY_END) == 0) && (nend == 0)) {
            if (current_body == NULL) {
                osh_error(EX_CONFIG,
                          "Error parsing geometry line %li - END encountered before any body definition\n",
                          (long int) lineno);
            }
            _save_body(current_body, nstr, par, npar, btype);
            current_body->lineno = lineno_b;
            break;
        }

        btype_new = _body_from_key(key);

        if (btype_new != OSH_GEMCA_BODY_NONE) {

            /* save previous body if any */
            if (current_body != NULL) {
                _save_body(current_body, nstr, par, npar, btype);
                current_body->lineno = lineno_b;
            }

            /* capacity guard */
            if (ibody >= g->nbodies) {
                osh_error(EX_CONFIG,
                          "Error parsing geometry line %li - too many bodies (max=%li)\n",
                          (long int) lineno,
                          (long int) g->nbodies);
            }

            /* start new body */
            current_body = g->bodies[ibody];
            btype = btype_new;
            lineno_b = lineno;

            if (args == NULL) {
                osh_error(
                    EX_CONFIG, "Error parsing geometry line %li - missing body name/parameters\n", (long int) lineno);
            }

            // TODO: sscanf("%s", ...) is unsafe because %s is unbounded and may overflow nstr.
            nt = sscanf(args, "%s %lf %lf %lf %lf %lf %lf", nstr, &par[0], &par[1], &par[2], &par[3], &par[4], &par[5]);

            /* check if body name is already used, and raise an error if that is the case */
            for (_ib = 0; _ib < ibody; _ib++) {
                if (strcmp(g->bodies[_ib]->name, nstr) == 0) {
                    osh_error(EX_CONFIG,
                              "Error parsing geometry line %li - body name '%s' already exists (defined at line %li)\n",
                              (long int) lineno,
                              nstr,
                              (long int) g->bodies[_ib]->lineno);
                }
            }

            npar = nt - 1;
            off = 6;

            ibody++;
        } else {
            /* continuation line */
            if (current_body == NULL) {
                osh_error(EX_CONFIG,
                          "Error parsing geometry line %li - parameters found before any body definition\n",
                          (long int) lineno);
            }

            /* We may write up to par[off+5] on a continuation line */
            if ((off + 5) >= OSH_GEMCA_NARGS_MAX) {
                osh_error(EX_CONFIG,
                          "Error parsing geometry line %li - too many arguments (need index %i, max index %i)\n",
                          (long int) lineno,
                          off + 5,
                          OSH_GEMCA_NARGS_MAX - 1);
            }

            nt = sscanf(key, "%lf", &par[off]);
            if (nt > 0) {
                npar += nt;
            }

            nt = 0;
            if (args != NULL) {
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

    free(line);
    line = NULL;
    return 1;
}

/**
 * @brief Checks if key is a valid body.
 *
 * @param[in] *key - character string holding the body 3-letter type
 *
 * @returns 0 (OSH_GEMCA_BODY_NONE) if no body was found, otherwise the type of the body. See osh_gemca2_defines.h
 *
 * @author Niels Bassler
 */
int _body_from_key(const char *key) {
    if (strcasecmp(key, OSH_GEMCA_KEY_SPH) == 0)
        return OSH_GEMCA_BODY_SPH;
    if (strcasecmp(key, OSH_GEMCA_KEY_WED) == 0)
        return OSH_GEMCA_BODY_WED;
    if (strcasecmp(key, OSH_GEMCA_KEY_ARB) == 0)
        return OSH_GEMCA_BODY_ARB;

    if (strcasecmp(key, OSH_GEMCA_KEY_BOX) == 0)
        return OSH_GEMCA_BODY_BOX;
    if (strcasecmp(key, OSH_GEMCA_KEY_VOX) == 0)
        return OSH_GEMCA_BODY_VOX;
    if (strcasecmp(key, OSH_GEMCA_KEY_RPP) == 0)
        return OSH_GEMCA_BODY_RPP;

    if (strcasecmp(key, OSH_GEMCA_KEY_RCC) == 0)
        return OSH_GEMCA_BODY_RCC;
    if (strcasecmp(key, OSH_GEMCA_KEY_REC) == 0)
        return OSH_GEMCA_BODY_REC;
    if (strcasecmp(key, OSH_GEMCA_KEY_TRC) == 0)
        return OSH_GEMCA_BODY_TRC;
    if (strcasecmp(key, OSH_GEMCA_KEY_ELL) == 0)
        return OSH_GEMCA_BODY_ELL;

    if (strcasecmp(key, OSH_GEMCA_KEY_YZP) == 0)
        return OSH_GEMCA_BODY_YZP;
    if (strcasecmp(key, OSH_GEMCA_KEY_XZP) == 0)
        return OSH_GEMCA_BODY_XZP;
    if (strcasecmp(key, OSH_GEMCA_KEY_XYP) == 0)
        return OSH_GEMCA_BODY_XYP;
    if (strcasecmp(key, OSH_GEMCA_KEY_PLA) == 0)
        return OSH_GEMCA_BODY_PLA;

    if (strcasecmp(key, OSH_GEMCA_KEY_ROT) == 0)
        return OSH_GEMCA_BODY_ROT;
    if (strcasecmp(key, OSH_GEMCA_KEY_CPY) == 0)
        return OSH_GEMCA_BODY_CPY;
    if (strcasecmp(key, OSH_GEMCA_KEY_MOV) == 0)
        return OSH_GEMCA_BODY_MOV;

    return OSH_GEMCA_BODY_NONE;
}

/**
 * @brief Store the body parameters into a body *b
 *
 * @details Memory will be allocated accordingly by this function.
 *
 * @param[in] *body - pointer to the body where the data will be stored
 * @param[in] *nstr - character string holding the user given body name
 * @param[in] *par - array with body arguments (double)
 * @param[in] *npar - length of *par array
 * @param[in] btype - type of body, see osh_gemca2_defines.h
 *
 * @returns always 1
 *
 * @author Niels Bassler
 */
static int _save_body(struct body *b, char *nstr, double *par, int npar, int btype) {

    int i = 0;

    b->type = btype;

    b->name = calloc(strlen(nstr) + 1, sizeof(char));
    snprintf(b->name, strlen(nstr) + 1, "%s", nstr);

    b->na = npar;
    b->a = calloc(npar, sizeof(double));

    memset(b->t, 0x00, 16 * sizeof(double)); // should not be necessary

    /* copy the arguments into the body array */
    for (i = 0; i < npar; i++) {
        b->a[i] = par[i];
    }
    return 1;
}
