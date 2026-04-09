#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_logger.h"
#include "scoring/parse/osh_scoring_parse_internal.h"
#include "scoring/parse/osh_scoring_parse_keys.h"

typedef enum osh_status (*geometry_handler_fn)(
    struct osh_scoring_geometry_def *, char **, int, char const *, unsigned int);

struct geometry_entry {
    char const *key;
    geometry_handler_fn handler;
};

static enum osh_status
append_axis(struct osh_scoring_geometry_def *geo, char const *label, double lo, double hi, int nbins);
static enum osh_status
geo_name(struct osh_scoring_geometry_def *geo, char **words, int nwords, char const *path, unsigned int lineno);
static enum osh_status
geo_axis(struct osh_scoring_geometry_def *geo, char **words, int nwords, char const *path, unsigned int lineno);
static enum osh_status
geo_rotation(struct osh_scoring_geometry_def *geo, char **words, int nwords, char const *path, unsigned int lineno);
static enum osh_status
geo_zones(struct osh_scoring_geometry_def *geo, char **words, int nwords, char const *path, unsigned int lineno);

static struct geometry_entry geometry_table[] = {{OSH_SCORING_KEY_NAME, geo_name},
                                                 {OSH_SCORING_KEY_GEO_X, geo_axis},
                                                 {OSH_SCORING_KEY_GEO_Y, geo_axis},
                                                 {OSH_SCORING_KEY_GEO_Z, geo_axis},
                                                 {OSH_SCORING_KEY_GEO_R, geo_axis},
                                                 {OSH_SCORING_KEY_GEO_ROT, geo_rotation},
                                                 {"rot", geo_rotation},
                                                 {OSH_SCORING_KEY_GEO_ZONES, geo_zones},
                                                 {NULL, NULL}};

enum osh_status osh_scoring_parse_geometry_line(struct osh_scoring_geometry_def *geo,
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
            return geometry_table[i].handler(geo, words, nwords, path, lineno);
        }
    }
    return OSH_OK;
}

static enum osh_status
append_axis(struct osh_scoring_geometry_def *geo, char const *label, double lo, double hi, int nbins) {
    struct osh_scoring_axis_def *tmp =
        (struct osh_scoring_axis_def *) realloc(geo->axes, (geo->naxes + 1u) * sizeof(*tmp));
    if (!tmp)
        return OSH_ENOMEM;
    geo->axes = tmp;
    memset(&geo->axes[geo->naxes], 0, sizeof(*tmp));
    strncpy(geo->axes[geo->naxes].label, label, sizeof(geo->axes[0].label) - 1u);
    geo->axes[geo->naxes].lo = lo;
    geo->axes[geo->naxes].hi = hi;
    geo->axes[geo->naxes].nbins = nbins;
    geo->naxes++;
    return OSH_OK;
}

static enum osh_status
geo_name(struct osh_scoring_geometry_def *geo, char **words, int nwords, char const *path, unsigned int lineno) {
    if (nwords < 2) {
        osh_error("%s:%u: Geometry Name requires an argument", path, lineno);
        return OSH_EPARSE;
    }
    free(geo->name);
    geo->name = strdup(words[1]);
    return geo->name ? OSH_OK : OSH_ENOMEM;
}

static enum osh_status
geo_axis(struct osh_scoring_geometry_def *geo, char **words, int nwords, char const *path, unsigned int lineno) {
    char label[2];

    if (nwords < 4) {
        osh_error("%s:%u: Geometry axis '%s' requires lo hi nbins", path, lineno, words[0]);
        return OSH_EPARSE;
    }
    label[0] = (char) toupper((unsigned char) words[0][0]);
    label[1] = '\0';
    return append_axis(geo, label, atof(words[1]), atof(words[2]), atoi(words[3]));
}

static enum osh_status
geo_rotation(struct osh_scoring_geometry_def *geo, char **words, int nwords, char const *path, unsigned int lineno) {
    if (nwords < 3) {
        osh_error("%s:%u: Geometry Rotation requires theta phi [deg]", path, lineno);
        return OSH_EPARSE;
    }
    geo->rot_theta_deg = atof(words[1]);
    geo->rot_phi_deg = atof(words[2]);
    geo->has_rotation = 1u;
    return OSH_OK;
}

static enum osh_status
geo_zones(struct osh_scoring_geometry_def *geo, char **words, int nwords, char const *path, unsigned int lineno) {
    if (nwords < 3) {
        osh_error("%s:%u: Geometry Zones requires start stop", path, lineno);
        return OSH_EPARSE;
    }
    geo->zone_start = atoi(words[1]);
    geo->zone_stop = atoi(words[2]);
    return OSH_OK;
}
