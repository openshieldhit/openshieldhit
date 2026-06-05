#include "apps/osh/osh_beam_spotlist.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_diag.h"
#include "common/osh_readline.h"
#include "openshieldhit/beam_defs.h"
#include "openshieldhit/const.h"

#define OSH_SPOTLIST_MAX_COLS 11

static enum osh_status spotlist_parse_values(char *line, double values[OSH_SPOTLIST_MAX_COLS], int *ncols_out);
static enum osh_status spotlist_detect_layout(FILE *fp, size_t *nspots_out, int *ncols_out);
static enum osh_status
spotlist_fill_spot(struct osh_beam_spot *spot, double const values[OSH_SPOTLIST_MAX_COLS], int ncols);

/**
 * @brief Parse one line of a spotlist file into an array of double values.
 *
 * @param[in,out] line      Line buffer; modified in place (NUL terminators may be inserted).
 * @param[out]    values    Receives up to OSH_SPOTLIST_MAX_COLS parsed values.
 * @param[out]    ncols_out Receives the number of values parsed from this line.
 *
 * @returns OSH_OK on success, OSH_EINVAL if too many columns or a non-numeric
 *          token is encountered.
 */
static enum osh_status spotlist_parse_values(char *line, double values[OSH_SPOTLIST_MAX_COLS], int *ncols_out) {
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

/**
 * @brief Scan the file to determine the column count and number of data rows.
 *
 * @details Rewinds @p fp before scanning. All non-blank, non-comment lines must
 * have the same column count; accepted counts are 5, 6, 7, 9, or 11.
 *
 * @param[in]  fp          Open file pointer (rewound on entry).
 * @param[out] nspots_out  Receives the number of data rows found.
 * @param[out] ncols_out   Receives the column count shared by all data rows.
 *
 * @returns OSH_OK on success, OSH_EINVAL if column counts are inconsistent or
 *          unsupported, OSH_EPARSE if no data rows are found, OSH_EIO on I/O error.
 */
static enum osh_status spotlist_detect_layout(FILE *fp, size_t *nspots_out, int *ncols_out) {
    char buff[OSH_MAX_LINE_LENGTH];
    double values[OSH_SPOTLIST_MAX_COLS];
    int line_ncols;
    int expected_ncols;
    size_t nspots;
    enum osh_status rc;

    if (!fp || !nspots_out || !ncols_out) {
        return OSH_EINVAL;
    }

    if (fseek(fp, 0L, SEEK_SET) != 0) {
        return OSH_EIO;
    }
    clearerr(fp);

    nspots = 0u;
    expected_ncols = 0;

    while (fgets(buff, sizeof buff, fp)) {
        rc = spotlist_parse_values(buff, values, &line_ncols);
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

    if (ferror(fp)) {
        return OSH_EIO;
    }
    if (nspots == 0u) {
        return OSH_EPARSE;
    }

    *nspots_out = nspots;
    *ncols_out = expected_ncols;
    return OSH_OK;
}

/**
 * @brief Populate a beam spot from one row of parsed spotlist values.
 *
 * @details Interprets the column layout determined by @p ncols (5, 6, 7, 9, or
 * 11). Energy is converted from GeV/nucleon to MeV/nucleon. Beam size is
 * converted from FWHM [cm] to sigma [cm]. Divergence is converted from mrad to
 * rad. Spot shape is set to Gaussian if either sigma is positive, otherwise pencil.
 *
 * @param[in,out] spot    Spot to fill; fields not covered by @p ncols are left
 *                        at their template values.
 * @param[in]     values  Array of at least @p ncols parsed doubles.
 * @param[in]     ncols   Column count; must be 5, 6, 7, 9, or 11.
 *
 * @returns OSH_OK on success, OSH_EINVAL for an unsupported column count or
 *          NULL pointer.
 */
static enum osh_status
spotlist_fill_spot(struct osh_beam_spot *spot, double const values[OSH_SPOTLIST_MAX_COLS], int ncols) {
    double sx;
    double sy;

    if (!spot || !values) {
        return OSH_EINVAL;
    }

    spot->t0 = values[0] * 1000.0;
    spot->p0 = 0.0;
    spot->t0_per_nucleon = 1;

    switch (ncols) {
    case 5:
        sx = values[3] * OSH_FWHM2SIGMA;
        sy = values[3] * OSH_FWHM2SIGMA;
        spot->wt = values[4];
        break;
    case 6:
        sx = values[3] * OSH_FWHM2SIGMA;
        sy = values[4] * OSH_FWHM2SIGMA;
        spot->wt = values[5];
        break;
    case 7:
        spot->tsigma = values[1] * 1000.0;
        spot->psigma = 0.0;
        spot->tsigma_per_nucleon = 1;
        sx = values[4] * OSH_FWHM2SIGMA;
        sy = values[5] * OSH_FWHM2SIGMA;
        spot->wt = values[6];
        break;
    case 9:
        spot->tsigma = values[1] * 1000.0;
        spot->psigma = 0.0;
        spot->tsigma_per_nucleon = 1;
        sx = values[4] * OSH_FWHM2SIGMA;
        sy = values[5] * OSH_FWHM2SIGMA;
        spot->div[0] = values[6] * 0.001;
        spot->div[1] = values[7] * 0.001;
        spot->wt = values[8];
        break;
    case 11:
        spot->tsigma = values[1] * 1000.0;
        spot->psigma = 0.0;
        spot->tsigma_per_nucleon = 1;
        sx = values[4] * OSH_FWHM2SIGMA;
        sy = values[5] * OSH_FWHM2SIGMA;
        spot->div[0] = values[6] * 0.001;
        spot->div[1] = values[7] * 0.001;
        spot->cor[0] = values[8];
        spot->cor[1] = values[9];
        spot->wt = values[10];
        break;
    default:
        return OSH_EINVAL;
    }

    /* x/y from the file are isocenter coordinates (SH12A / DICOM RT Plan
     * convention).  They are stored here as-is; the caller (osh_app_osh.c)
     * is responsible for back-projecting them to physical beam-start
     * coordinates using BEAMSAD before the spots enter the beam workspace. */
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

enum osh_status osh_beam_spotlist_import(char const *path,
                                         struct osh_diag_sink const *diag,
                                         struct osh_beam_spot const *template_spot,
                                         struct osh_beam_spot **spots_out,
                                         size_t *nspots_out) {
    struct osh_beam_spot *spots_new;
    FILE *fp;
    char buff[OSH_MAX_LINE_LENGTH];
    double values[OSH_SPOTLIST_MAX_COLS];
    size_t nspots;
    size_t i;
    int ncols;
    int line_ncols;
    enum osh_status rc;

    if (!path || !template_spot || !spots_out || !nspots_out) {
        return OSH_EINVAL;
    }
    *spots_out = NULL;
    *nspots_out = 0u;

    fp = fopen(path, "r");
    if (!fp) {
        return OSH_EIO;
    }

    OSH_DIAG_INFOF(diag, "Loading external spotlist: %s", path);

    rc = spotlist_detect_layout(fp, &nspots, &ncols);
    if (rc != OSH_OK) {
        fclose(fp);
        return rc;
    }

    OSH_DIAG_INFOF(diag, "External spotlist format: %d columns, %lu spots", ncols, (unsigned long) nspots);

    spots_new = (struct osh_beam_spot *) calloc(nspots, sizeof(*spots_new));
    if (!spots_new) {
        fclose(fp);
        return OSH_ENOMEM;
    }

    if (fseek(fp, 0L, SEEK_SET) != 0) {
        fclose(fp);
        free(spots_new);
        return OSH_EIO;
    }
    clearerr(fp);

    i = 0u;
    while (fgets(buff, sizeof buff, fp)) {
        rc = spotlist_parse_values(buff, values, &line_ncols);
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

        spots_new[i] = *template_spot;
        rc = spotlist_fill_spot(&spots_new[i], values, ncols);
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

    *spots_out = spots_new;
    *nspots_out = nspots;
    return OSH_OK;
}
