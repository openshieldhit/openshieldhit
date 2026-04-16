#include "beam/osh_beam_spots.h"

#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_logger.h"
#include "common/osh_rc.h"
#include "common/osh_readline.h"

#define OSH_SPOTLIST_MAX_COLS 11

static int _spotlist_parse_values(char *line, double values[OSH_SPOTLIST_MAX_COLS], int *ncols_out);
static int _spotlist_detect_layout(FILE *fp, size_t *nspots_out, int *ncols_out);
static int _spotlist_fill_spot(struct beam_spot *spot, double const values[OSH_SPOTLIST_MAX_COLS], int ncols);

static int _spotlist_parse_values(char *line, double values[OSH_SPOTLIST_MAX_COLS], int *ncols_out) {
    char *p;
    char *endp;
    int ncols;

    if (!line || !values || !ncols_out) {
        return OSH_EINVAL;
    }

    p = line;
    ncols = 0;

    while (*p) {
        while (isspace((unsigned char) *p)) {
            p++;
        }
        if (*p == '\0' || strchr(OSH_READLINE_COMMENT, (unsigned char) *p)) {
            break;
        }
        if (ncols >= OSH_SPOTLIST_MAX_COLS) {
            return OSH_EINVAL;
        }

        values[ncols] = strtod(p, &endp);
        if (endp == p) {
            return OSH_EINVAL;
        }

        ncols++;
        p = endp;
    }

    *ncols_out = ncols;
    return OSH_OK;
}

static int _spotlist_detect_layout(FILE *fp, size_t *nspots_out, int *ncols_out) {
    char buff[OSH_MAX_LINE_LENGTH];
    double values[OSH_SPOTLIST_MAX_COLS];
    int line_ncols;
    int expected_ncols;
    size_t nspots;
    int rc;

    if (!fp || !nspots_out || !ncols_out) {
        return OSH_EINVAL;
    }

    rewind(fp);
    clearerr(fp);

    nspots = 0;
    expected_ncols = 0;

    while (fgets(buff, sizeof buff, fp)) {
        rc = _spotlist_parse_values(buff, values, &line_ncols);
        if (rc != OSH_OK) {
            return rc;
        }
        if (line_ncols == 0) {
            continue;
        }
        switch (line_ncols) {
        case 5:
        case 6:
        case 7:
        case 9:
        case 11:
            break;
        default:
            return OSH_EINVAL;
        }

        if (expected_ncols == 0) {
            expected_ncols = line_ncols;
        } else if (line_ncols != expected_ncols) {
            return OSH_EINVAL;
        }
        nspots++;
    }

    if (ferror(fp) || nspots == 0) {
        return OSH_EIO;
    }

    *nspots_out = nspots;
    *ncols_out = expected_ncols;
    return OSH_OK;
}

static int _spotlist_fill_spot(struct beam_spot *spot, double const values[OSH_SPOTLIST_MAX_COLS], int ncols) {
    double sx;
    double sy;
    double fwhm_to_sigma;

    if (!spot || !values) {
        return OSH_EINVAL;
    }

    fwhm_to_sigma = 1.0 / (2.0 * sqrt(2.0 * log(2.0)));

    spot->t0 = values[0] * 1000.0;
    spot->p0 = 0.0;
    spot->t0_per_nucleon = 1;

    switch (ncols) {
    case 5:
        sx = values[3] * fwhm_to_sigma;
        sy = values[3] * fwhm_to_sigma;
        spot->wt = values[4];
        break;

    case 6:
        sx = values[3] * fwhm_to_sigma;
        sy = values[4] * fwhm_to_sigma;
        spot->wt = values[5];
        break;

    case 7:
        spot->tsigma = values[1] * 1000.0;
        spot->psigma = 0.0;
        spot->tsigma_per_nucleon = 1;
        sx = values[4] * fwhm_to_sigma;
        sy = values[5] * fwhm_to_sigma;
        spot->wt = values[6];
        break;

    case 9:
        spot->tsigma = values[1] * 1000.0;
        spot->psigma = 0.0;
        spot->tsigma_per_nucleon = 1;
        sx = values[4] * fwhm_to_sigma;
        sy = values[5] * fwhm_to_sigma;
        spot->div[0] = values[6] * 0.001;
        spot->div[1] = values[7] * 0.001;
        spot->wt = values[8];
        break;

    case 11:
        spot->tsigma = values[1] * 1000.0;
        spot->psigma = 0.0;
        spot->tsigma_per_nucleon = 1;
        sx = values[4] * fwhm_to_sigma;
        sy = values[5] * fwhm_to_sigma;
        spot->div[0] = values[6] * 0.001;
        spot->div[1] = values[7] * 0.001;
        spot->cor[0] = values[8];
        spot->cor[1] = values[9];
        spot->wt = values[10];
        break;

    default:
        return OSH_EINVAL;
    }

    spot->p[0] = values[(ncols == 5 || ncols == 6) ? 1 : 2];
    spot->p[1] = values[(ncols == 5 || ncols == 6) ? 2 : 3];

    if (sx > 0.0 || sy > 0.0) {
        spot->shape = OSH_BEAM_SHAPE_GAUSSIAN;
        spot->size[0] = sx;
        spot->size[1] = sy;
    } else {
        spot->shape = OSH_BEAM_SHAPE_PENCIL;
        spot->size[0] = 0.0;
        spot->size[1] = 0.0;
    }

    return OSH_OK;
}

int osh_beam_spots_init(struct beam_spot **sl, size_t nspots) {
    if (!sl || nspots == 0) {
        return OSH_EINVAL;
    }

    *sl = calloc(nspots, sizeof(struct beam_spot));
    if (!*sl) {
        return OSH_ENOMEM;
    }

    return OSH_OK;
}

int osh_beam_spots_free(struct beam_spot *sl) {
    if (!sl) {
        return OSH_EINVAL;
    }
    free(sl);
    return OSH_OK;
}

int osh_beam_spotlist_load(struct beam_workspace *beam, char const *spotlist_path) {
    struct beam_spot template_spot;
    struct beam_spot *spots_new;
    FILE *fp;
    char buff[OSH_MAX_LINE_LENGTH];
    double values[OSH_SPOTLIST_MAX_COLS];
    size_t nspots;
    size_t i;
    int ncols;
    int line_ncols;
    int rc;

    if (!beam || !spotlist_path || !beam->spots || beam->nspots == 0) {
        return OSH_EINVAL;
    }

    fp = fopen(spotlist_path, "r");
    if (!fp) {
        return OSH_EIO;
    }

    osh_info("Loading external spotlist: %s", spotlist_path);

    rc = _spotlist_detect_layout(fp, &nspots, &ncols);
    if (rc != OSH_OK) {
        fclose(fp);
        return rc;
    }

    osh_info("External spotlist format: %d columns, %lu spots", ncols, (unsigned long) nspots);

    rc = osh_beam_spots_init(&spots_new, nspots);
    if (rc != OSH_OK) {
        fclose(fp);
        return rc;
    }

    template_spot = beam->spots[0];

    rewind(fp);
    clearerr(fp);

    i = 0;
    while (fgets(buff, sizeof buff, fp)) {
        rc = _spotlist_parse_values(buff, values, &line_ncols);
        if (rc != OSH_OK) {
            fclose(fp);
            free(spots_new);
            return rc;
        }
        if (line_ncols == 0) {
            continue;
        }
        if (line_ncols != ncols || i >= nspots) {
            fclose(fp);
            free(spots_new);
            return OSH_EINVAL;
        }

        spots_new[i] = template_spot;
        rc = _spotlist_fill_spot(&spots_new[i], values, ncols);
        if (rc != OSH_OK) {
            fclose(fp);
            free(spots_new);
            return rc;
        }
        spots_new[i].spot_id = (unsigned int) (i + 1u);
        i++;
    }

    fclose(fp);

    if (i != nspots) {
        free(spots_new);
        return OSH_EIO;
    }

    free(beam->spots);
    beam->spots = spots_new;
    beam->nspots = nspots;

    for (i = 0; i < beam->nspots; i++) {
        if (fabs(beam->spots[i].div[0]) > 0.0 || fabs(beam->spots[i].div[1]) > 0.0) {
            beam->shared.use_div = 1;
            break;
        }
    }

    return OSH_OK;
}

int osh_beam_shared_init(struct beam_shared *shared) {
    if (!shared) {
        return OSH_EINVAL;
    }
    memset(shared, 0, sizeof *shared);

    shared->sad[0] = 0.0;
    shared->sad[1] = 0.0;
    shared->focus = 0.0;
    shared->theta = 0.0; /* rad — along +Z by default */
    shared->phi = 0.0;   /* rad */
    shared->use_div = 0;
    shared->use_sad = 0;

    return OSH_OK;
}
