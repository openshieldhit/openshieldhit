/**
 * @file osh_scoring_parse_settings.c
 *
 * @brief Parse one tokenized line inside a scoring `Settings` section.
 *
 * @details
 * The settings parser accepts canonical keywords plus a few legacy aliases
 * (`sitediam`, `rho`, `npart`) and stores values in the active settings block.
 * Presence flags (`has_*`) are set per field so later stages can distinguish
 * explicit user input from defaults.
 */

#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_scoring_parse_internal.h"
#include "apps/osh/osh_scoring_parse_keys.h"
#include "openshieldhit/logger.h"

typedef enum osh_status (*settings_handler_fn)(
    struct osh_scoring_settings_def *, char **, int, char const *, unsigned int);

struct settings_entry {
    char const *key;
    settings_handler_fn handler;
};

static enum osh_status
settings_name(struct osh_scoring_settings_def *set, char **words, int nwords, char const *path, unsigned int lineno);
static enum osh_status
settings_rescale(struct osh_scoring_settings_def *set, char **words, int nwords, char const *path, unsigned int lineno);
static enum osh_status
settings_offset(struct osh_scoring_settings_def *set, char **words, int nwords, char const *path, unsigned int lineno);
static enum osh_status
settings_medium(struct osh_scoring_settings_def *set, char **words, int nwords, char const *path, unsigned int lineno);
static enum osh_status settings_nkmedium(
    struct osh_scoring_settings_def *set, char **words, int nwords, char const *path, unsigned int lineno);
static enum osh_status settings_sitediam(
    struct osh_scoring_settings_def *set, char **words, int nwords, char const *path, unsigned int lineno);
static enum osh_status
settings_density(struct osh_scoring_settings_def *set, char **words, int nwords, char const *path, unsigned int lineno);
static enum osh_status settings_maxcount(
    struct osh_scoring_settings_def *set, char **words, int nwords, char const *path, unsigned int lineno);

static struct settings_entry settings_table[] = {{OSH_SCORING_KEY_NAME, settings_name},
                                                 {OSH_SCORING_KEY_RESCALE, settings_rescale},
                                                 {OSH_SCORING_KEY_OFFSET, settings_offset},
                                                 {OSH_SCORING_KEY_MEDIUM, settings_medium},
                                                 {OSH_SCORING_KEY_NKMEDIUM, settings_nkmedium},
                                                 {OSH_SCORING_KEY_SITEDIAM, settings_sitediam},
                                                 {"sitediam", settings_sitediam},
                                                 {"rho", settings_density},
                                                 {OSH_SCORING_KEY_DENSITY, settings_density},
                                                 {OSH_SCORING_KEY_MAXCOUNT, settings_maxcount},
                                                 {"npart", settings_maxcount},
                                                 {NULL, NULL}};

/**
 * @brief Dispatch one tokenized line into the active settings definition.
 */
enum osh_status osh_scoring_parse_settings_line(struct osh_scoring_settings_def *set,
                                                char **words,
                                                int nwords,
                                                char const *path,
                                                unsigned int lineno,
                                                int *found_out) {
    size_t i;

    if (found_out)
        *found_out = 0;
    for (i = 0; settings_table[i].key != NULL; ++i) {
        if (strcmp(settings_table[i].key, words[0]) == 0) {
            if (found_out)
                *found_out = 1;
            return settings_table[i].handler(set, words, nwords, path, lineno);
        }
    }
    return OSH_OK;
}

/**
 * @brief Parse `Name <value>`.
 */
static enum osh_status
settings_name(struct osh_scoring_settings_def *set, char **words, int nwords, char const *path, unsigned int lineno) {
    if (nwords < 2) {
        osh_error("%s:%u: Settings Name requires an argument", path, lineno);
        return OSH_EPARSE;
    }
    free(set->name);
    set->name = strdup(words[1]);
    return set->name ? OSH_OK : OSH_ENOMEM;
}

/**
 * @brief Parse `Rescale <value>`.
 */
static enum osh_status settings_rescale(
    struct osh_scoring_settings_def *set, char **words, int nwords, char const *path, unsigned int lineno) {
    if (nwords < 2) {
        osh_error("%s:%u: Rescale requires a value", path, lineno);
        return OSH_EPARSE;
    }
    set->rescale = atof(words[1]);
    set->has_rescale = 1u;
    return OSH_OK;
}

/**
 * @brief Parse `Offset <value>`.
 */
static enum osh_status
settings_offset(struct osh_scoring_settings_def *set, char **words, int nwords, char const *path, unsigned int lineno) {
    if (nwords < 2) {
        osh_error("%s:%u: Offset requires a value", path, lineno);
        return OSH_EPARSE;
    }
    set->offset = atof(words[1]);
    set->has_offset = 1u;
    return OSH_OK;
}

/**
 * @brief Parse `Medium <id>`.
 */
static enum osh_status
settings_medium(struct osh_scoring_settings_def *set, char **words, int nwords, char const *path, unsigned int lineno) {
    if (nwords < 2) {
        osh_error("%s:%u: Medium requires a value", path, lineno);
        return OSH_EPARSE;
    }
    set->medium = atoi(words[1]);
    set->has_medium = 1u;
    return OSH_OK;
}

/**
 * @brief Parse `NKMedium <id>`.
 */
static enum osh_status settings_nkmedium(
    struct osh_scoring_settings_def *set, char **words, int nwords, char const *path, unsigned int lineno) {
    if (nwords < 2) {
        osh_error("%s:%u: NKMedium requires a value", path, lineno);
        return OSH_EPARSE;
    }
    set->nkmedium = atoi(words[1]);
    set->has_nkmedium = 1u;
    return OSH_OK;
}

/**
 * @brief Parse `SiteDiameter <value_um>` (and alias `sitediam`).
 */
static enum osh_status settings_sitediam(
    struct osh_scoring_settings_def *set, char **words, int nwords, char const *path, unsigned int lineno) {
    if (nwords < 2) {
        osh_error("%s:%u: SiteDiameter requires a value", path, lineno);
        return OSH_EPARSE;
    }
    set->site_diameter_um = atof(words[1]);
    set->has_site_diameter_um = 1u;
    return OSH_OK;
}

/**
 * @brief Parse `Density <value_g_cm3>` (and alias `rho`).
 */
static enum osh_status settings_density(
    struct osh_scoring_settings_def *set, char **words, int nwords, char const *path, unsigned int lineno) {
    if (nwords < 2) {
        osh_error("%s:%u: Density requires a value", path, lineno);
        return OSH_EPARSE;
    }
    set->density_g_cm3 = atof(words[1]);
    set->has_density_g_cm3 = 1u;
    return OSH_OK;
}

/**
 * @brief Parse `MaxCount <n>` (and alias `npart`).
 */
static enum osh_status settings_maxcount(
    struct osh_scoring_settings_def *set, char **words, int nwords, char const *path, unsigned int lineno) {
    if (nwords < 2) {
        osh_error("%s:%u: MaxCount requires a value", path, lineno);
        return OSH_EPARSE;
    }
    set->npart = (size_t) strtoull(words[1], NULL, 10);
    set->has_npart = 1u;
    return OSH_OK;
}
