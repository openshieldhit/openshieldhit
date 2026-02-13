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
    int lineno_b; /* save first line where this body was defined */

    int nt = 0;
    int btype = 0;
    int btype_new = 0;

    int nend = 0; /* number of end statements parsed so far */

    double par[OSH_GEMCA_NARGS_MAX]; /* temporary placeholder for body parsed arguments */
    int npar = 0;                    /* number of parameters */
    int off = 0;                     /* index offset for parameters while parsing */

    char nstr[OSH_GEMCA_BODY_NAME_MAXLEN]; /* string for body name */

    ssize_t ibody = -1; // TODO change to size_t
    size_t _ib;

    rewind(shf->fp);

    while (osh_readline_key(shf, &line, &key, &args, &lineno) > 0) {

        // printf("line %i  '%s':'%s'\n", lineno, key, args);
        btype_new = _body_from_key(key);
        if (OSH_GEMCA_BODY_NONE != btype_new) {
            /* we have found a new body */

            /* let us save the previous one, if there was one */
            if (ibody > -1) {
                _save_body(g->bodies[ibody], nstr, par, npar, btype);
                /* save the first line number where the body was defined */
                g->bodies[ibody]->lineno = lineno_b;
            }

            ibody++;
            // TODO: improve parser so all arguments may appear in a single line, if they want to
            // This corresponds to FLUKA-free format
            nt = sscanf(
                args, "%s %lf %lf %lf %lf %lf %lf\n", nstr, &par[0], &par[1], &par[2], &par[3], &par[4], &par[5]);
            /*
               // printf("parsed : %i parameters and incl. one name\n", nt);
               // printf("parsed : %lf %lf %lf - %lf %lf %lf\n", par[0], par[1], par[2], par[3], par[4], par[5]);
             */

            /* check if a body with name nstr already exist, and if yes, exit with an error */
            for (_ib = 0; _ib < (size_t) ibody; _ib++) {
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
            btype = btype_new;
        } else {
            if (btype_new == OSH_GEMCA_BODY_VOX) {
                /* this is the special case when we have a path to a filename */
                g->bodies[ibody]->filename_vox = calloc(strlen(key) + 1, sizeof(char));
                snprintf(g->bodies[ibody]->filename_vox, strlen(key), "%s", key);
                printf("Vox filename = %s\n", key);
            } else {
                nt = sscanf(key, "%lf", &par[off]);
            }

            /* check if there is a continuation line to the same body.
                This is the case if the next line is a number.  */
            if (nt > 0) {
                /* yes, there is, extend the number of parameters */
                npar += nt;
                // printf("parsed another : %i parameters\n", nt);
            }

            /* scan for remaining parameters in this line */
            nt = 0;
            if (args != NULL) {
                nt = sscanf(args,
                            "%lf %lf %lf %lf %lf\n",
                            &par[1 + off],
                            &par[2 + off],
                            &par[3 + off],
                            &par[4 + off],
                            &par[5 + off]);
            }
            if (nt > 0) {
                // printf("parsed another : %i parameters\n", nt);
                /* and update the number of arguments parsed */
                npar += nt;
            }
            off += 6;
        }
        if ((off + 5) > OSH_GEMCA_NARGS_MAX) {
            osh_error(EX_CONFIG,
                      "Error parsing geometry line %li - too many arguments  off+5 = %i\n",
                      (long int) lineno,
                      off + 5);
        }

        if ((strcasecmp(key, OSH_GEMCA_KEY_END) == 0) && (nend == 0)) {
            /*  save the last body we found */
            _save_body(g->bodies[ibody], nstr, par, npar, btype);
            break; /* break out of while loop */
        }
        free(line);
        /* save the line number before reading next body*/
        lineno_b = lineno;

    } /* end of while loop */
    free(line);
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
