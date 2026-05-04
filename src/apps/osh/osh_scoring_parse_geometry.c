/**
 * @file osh_scoring_parse_geometry.c
 *
 * @brief Parse one tokenized line inside a scoring `Geometry` section.
 *
 * @details
 * Recognized keys:
 * - `Name <string>`
 * - axis lines: `x|y|z|r <lo> <hi> <nbins>`
 * - `Rotation <theta_deg> <phi_deg>` (and legacy alias `rot`)
 * - `Zones <start> <stop>`
 *
 * The section parser stores values in raw parsed form; consistency checks
 * against geometry type happen in later validation/finalization stages.
 */

#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_scoring_parse_internal.h"
#include "apps/osh/osh_scoring_parse_keys.h"
#include "common/osh_diag.h"

typedef enum osh_status (*geometry_handler_fn)(
    struct osh_scoring_geometry_def *, struct osh_diag_sink const *, char **, int, char const *, unsigned int);

struct geometry_entry {
    char const *key;
    geometry_handler_fn handler;
};

static enum osh_status
append_axis(struct osh_scoring_geometry_def *geo, char const *label, double lo, double hi, int nbins);
static enum osh_status geo_name(struct osh_scoring_geometry_def *geo,
                                struct osh_diag_sink const *diag,
                                char **words,
                                int nwords,
                                char const *path,
                                unsigned int lineno);
static enum osh_status geo_axis(struct osh_scoring_geometry_def *geo,
                                struct osh_diag_sink const *diag,
                                char **words,
                                int nwords,
                                char const *path,
                                unsigned int lineno);
static enum osh_status geo_rotation(struct osh_scoring_geometry_def *geo,
                                    struct osh_diag_sink const *diag,
                                    char **words,
                                    int nwords,
                                    char const *path,
                                    unsigned int lineno);
static enum osh_status geo_zones(struct osh_scoring_geometry_def *geo,
                                 struct osh_diag_sink const *diag,
                                 char **words,
                                 int nwords,
                                 char const *path,
                                 unsigned int lineno);
static enum osh_status geo_inputpath(struct osh_scoring_geometry_def *geo,
                                     struct osh_diag_sink const *diag,
                                     char **words,
                                     int nwords,
                                     char const *path,
                                     unsigned int lineno);
static enum osh_status geo_body(struct osh_scoring_geometry_def *geo,
                                struct osh_diag_sink const *diag,
                                char **words,
                                int nwords,
                                char const *path,
                                unsigned int lineno);

static struct geometry_entry geometry_table[] = {{OSH_SCORING_KEY_NAME, geo_name},
                                                 {OSH_SCORING_KEY_GEO_X, geo_axis},
                                                 {OSH_SCORING_KEY_GEO_Y, geo_axis},
                                                 {OSH_SCORING_KEY_GEO_Z, geo_axis},
                                                 {OSH_SCORING_KEY_GEO_R, geo_axis},
                                                 {OSH_SCORING_KEY_GEO_ROT, geo_rotation},
                                                 {"rot", geo_rotation},
                                                 {OSH_SCORING_KEY_GEO_ZONES, geo_zones},
                                                 {OSH_SCORING_KEY_GEO_INPUTPATH, geo_inputpath},
                                                 {OSH_SCORING_KEY_GEO_BODY, geo_body},
                                                 {NULL, NULL}};

/**
 * @brief Dispatch one tokenized line into the active geometry definition.
 */
enum osh_status osh_scoring_parse_geometry_line(struct osh_scoring_geometry_def *geo,
                                                struct osh_diag_sink const *diag,
                                                char **words,
                                                int nwords,
                                                char const *path,
                                                unsigned int lineno,
                                                int *found_out) {
    size_t i;

    if (found_out)
        *found_out = 0;
    for (i = 0; geometry_table[i].key != NULL; ++i) {
        if (strcmp(geometry_table[i].key, words[0]) == 0) {
            if (found_out)
                *found_out = 1;
            return geometry_table[i].handler(geo, diag, words, nwords, path, lineno);
        }
    }
    return OSH_OK;
}

/**
 * @brief Append one axis specification to the geometry.
 */
static enum osh_status
append_axis(struct osh_scoring_geometry_def *geo, char const *label, double lo, double hi, int nbins) {
    size_t len;
    struct osh_scoring_axis_def *tmp =
        (struct osh_scoring_axis_def *) realloc(geo->axes, (geo->naxes + 1u) * sizeof(*tmp));
    if (!tmp)
        return OSH_ENOMEM;
    geo->axes = tmp;
    memset(&geo->axes[geo->naxes], 0, sizeof(*tmp));
    len = strlen(label);
    if (len >= sizeof(geo->axes[0].label)) {
        len = sizeof(geo->axes[0].label) - 1u;
    }
    memcpy(geo->axes[geo->naxes].label, label, len);
    geo->axes[geo->naxes].label[len] = '\0';
    geo->axes[geo->naxes].lo = lo;
    geo->axes[geo->naxes].hi = hi;
    geo->axes[geo->naxes].nbins = nbins;
    geo->naxes++;
    return OSH_OK;
}

/**
 * @brief Parse `Name <value>` for a geometry section.
 */
static enum osh_status geo_name(struct osh_scoring_geometry_def *geo,
                                struct osh_diag_sink const *diag,
                                char **words,
                                int nwords,
                                char const *path,
                                unsigned int lineno) {
    if (nwords < 2) {
        OSH_DIAG_ERRORF(diag, "%s:%u: Geometry Name requires an argument", path, lineno);
        return OSH_EPARSE;
    }
    free(geo->name);
    geo->name = strdup(words[1]);
    return geo->name ? OSH_OK : OSH_ENOMEM;
}

/**
 * @brief Parse one axis line (`x|y|z|r lo hi nbins`).
 */
static enum osh_status geo_axis(struct osh_scoring_geometry_def *geo,
                                struct osh_diag_sink const *diag,
                                char **words,
                                int nwords,
                                char const *path,
                                unsigned int lineno) {
    char label[2];

    if (nwords < 4) {
        OSH_DIAG_ERRORF(diag, "%s:%u: Geometry axis '%s' requires lo hi nbins", path, lineno, words[0]);
        return OSH_EPARSE;
    }
    label[0] = (char) toupper((unsigned char) words[0][0]);
    label[1] = '\0';
    return append_axis(geo, label, atof(words[1]), atof(words[2]), atoi(words[3]));
}

/**
 * @brief Parse optional geometry rotation (`theta phi` in degrees).
 */
static enum osh_status geo_rotation(struct osh_scoring_geometry_def *geo,
                                    struct osh_diag_sink const *diag,
                                    char **words,
                                    int nwords,
                                    char const *path,
                                    unsigned int lineno) {
    if (nwords < 3) {
        OSH_DIAG_ERRORF(diag, "%s:%u: Geometry Rotation requires theta phi [deg]", path, lineno);
        return OSH_EPARSE;
    }
    geo->rot_theta_deg = atof(words[1]);
    geo->rot_phi_deg = atof(words[2]);
    geo->has_rotation = 1u;
    return OSH_OK;
}

/**
 * @brief Parse zone-range form used by `zone` scoring geometries.
 */
static enum osh_status geo_zones(struct osh_scoring_geometry_def *geo,
                                 struct osh_diag_sink const *diag,
                                 char **words,
                                 int nwords,
                                 char const *path,
                                 unsigned int lineno) {
    if (nwords < 3) {
        OSH_DIAG_ERRORF(diag, "%s:%u: Geometry Zones requires start stop", path, lineno);
        return OSH_EPARSE;
    }
    geo->zone_start = atoi(words[1]);
    geo->zone_stop = atoi(words[2]);
    return OSH_OK;
}

/**
 * @brief Parse `InputPath <path>` for DicomRTDOSE geometry sections.
 *
 * @details
 * The path is stored as-is on the cold geometry def.  The app layer reads the
 * file at orchestration time and populates vox_nbins before scoring compile.
 * The library never opens this path.
 */
static enum osh_status geo_inputpath(struct osh_scoring_geometry_def *geo,
                                     struct osh_diag_sink const *diag,
                                     char **words,
                                     int nwords,
                                     char const *path,
                                     unsigned int lineno) {
    if (nwords < 2) {
        OSH_DIAG_ERRORF(diag, "%s:%u: Geometry InputPath requires a file path", path, lineno);
        return OSH_EPARSE;
    }
    free(geo->vox_rtdose_path);
    geo->vox_rtdose_path = strdup(words[1]);
    return geo->vox_rtdose_path ? OSH_OK : OSH_ENOMEM;
}

/**
 * @brief Parse `Body <name>` for DicomCT / DicomRTDOSE geometry sections.
 *
 * @details
 * Stores the referenced CT body name for future multi-CT support.  Currently
 * the app auto-detects the single CT body in the geometry and ignores this
 * field, but it is validated and preserved so users can prepare their input
 * files for future use.
 */
static enum osh_status geo_body(struct osh_scoring_geometry_def *geo,
                                struct osh_diag_sink const *diag,
                                char **words,
                                int nwords,
                                char const *path,
                                unsigned int lineno) {
    if (nwords < 2) {
        OSH_DIAG_ERRORF(diag, "%s:%u: Geometry Body requires a body name", path, lineno);
        return OSH_EPARSE;
    }
    free(geo->vox_body_name);
    geo->vox_body_name = strdup(words[1]);
    return geo->vox_body_name ? OSH_OK : OSH_ENOMEM;
}
