/**
 * @file osh_scoring_parse_filter.c
 *
 * @brief Parse one tokenized line inside a scoring `Filter` section.
 *
 * @details
 * Recognized keys:
 * - `Name <string>`
 * - `z|a|e|gen|id <op> <value>`
 *
 * Rule keys are normalized to uppercase field names and appended to the
 * current filter as raw text/number rules. Semantic validation of operators
 * and field/value compatibility is done later during scoring finalize.
 */

#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_scoring_parse_internal.h"
#include "apps/osh/osh_scoring_parse_keys.h"
#include "common/osh_diag.h"

typedef enum osh_status (*filter_handler_fn)(
    struct osh_scoring_filter_def *, struct osh_diag_sink const *, char **, int, char const *, unsigned int);

struct filter_entry {
    char const *key;
    filter_handler_fn handler;
};

static enum osh_status
append_filter_rule(struct osh_scoring_filter_def *fil, char const *field, char const *op, double value);
static enum osh_status filter_name(struct osh_scoring_filter_def *fil,
                                   struct osh_diag_sink const *diag,
                                   char **words,
                                   int nwords,
                                   char const *path,
                                   unsigned int lineno);
static enum osh_status filter_rule(struct osh_scoring_filter_def *fil,
                                   struct osh_diag_sink const *diag,
                                   char **words,
                                   int nwords,
                                   char const *path,
                                   unsigned int lineno);

static struct filter_entry filter_table[] = {{OSH_SCORING_KEY_NAME, filter_name},
                                             {OSH_SCORING_KEY_FILTER_Z, filter_rule},
                                             {OSH_SCORING_KEY_FILTER_A, filter_rule},
                                             {OSH_SCORING_KEY_FILTER_E, filter_rule},
                                             {OSH_SCORING_KEY_FILTER_GEN, filter_rule},
                                             {OSH_SCORING_KEY_FILTER_ID, filter_rule},
                                             {NULL, NULL}};

/**
 * @brief Dispatch one tokenized line into the active filter definition.
 */
enum osh_status osh_scoring_parse_filter_line(struct osh_scoring_filter_def *fil,
                                              struct osh_diag_sink const *diag,
                                              char **words,
                                              int nwords,
                                              char const *path,
                                              unsigned int lineno,
                                              int *found_out) {
    size_t i;

    if (found_out)
        *found_out = 0;
    for (i = 0; filter_table[i].key != NULL; ++i) {
        if (strcmp(filter_table[i].key, words[0]) == 0) {
            if (found_out)
                *found_out = 1;
            return filter_table[i].handler(fil, diag, words, nwords, path, lineno);
        }
    }
    return OSH_OK;
}

/**
 * @brief Append one parsed rule to the filter's dynamic rule array.
 */
static enum osh_status
append_filter_rule(struct osh_scoring_filter_def *fil, char const *field, char const *op, double value) {
    size_t len;
    struct osh_scoring_filter_rule *tmp =
        (struct osh_scoring_filter_rule *) realloc(fil->rules, (fil->nrules + 1u) * sizeof(*tmp));
    if (!tmp)
        return OSH_ENOMEM;
    fil->rules = tmp;
    memset(&fil->rules[fil->nrules], 0, sizeof(*tmp));
    len = strlen(field);
    if (len >= sizeof(fil->rules[0].field)) {
        len = sizeof(fil->rules[0].field) - 1u;
    }
    memcpy(fil->rules[fil->nrules].field, field, len);
    fil->rules[fil->nrules].field[len] = '\0';
    len = strlen(op);
    if (len >= sizeof(fil->rules[0].op)) {
        len = sizeof(fil->rules[0].op) - 1u;
    }
    memcpy(fil->rules[fil->nrules].op, op, len);
    fil->rules[fil->nrules].op[len] = '\0';
    fil->rules[fil->nrules].value = value;
    fil->nrules++;
    return OSH_OK;
}

/**
 * @brief Parse `Name <value>` for a filter.
 */
static enum osh_status filter_name(struct osh_scoring_filter_def *fil,
                                   struct osh_diag_sink const *diag,
                                   char **words,
                                   int nwords,
                                   char const *path,
                                   unsigned int lineno) {
    if (nwords < 2) {
        OSH_DIAG_ERRORF(diag, "%s:%u: Filter Name requires an argument", path, lineno);
        return OSH_EPARSE;
    }
    free(fil->name);
    fil->name = strdup(words[1]);
    return fil->name ? OSH_OK : OSH_ENOMEM;
}

/**
 * @brief Parse `<field> <op> <value>` and append one filter rule.
 */
static enum osh_status filter_rule(struct osh_scoring_filter_def *fil,
                                   struct osh_diag_sink const *diag,
                                   char **words,
                                   int nwords,
                                   char const *path,
                                   unsigned int lineno) {
    char field[16];
    char *endptr;
    double value;
    size_t i;

    if (nwords < 3) {
        OSH_DIAG_ERRORF(diag, "%s:%u: filter rule requires: <field> <op> <value>", path, lineno);
        return OSH_EPARSE;
    }
    value = strtod(words[2], &endptr);
    if (endptr == words[2]) {
        OSH_DIAG_ERRORF(diag, "%s:%u: filter rule value '%s' is not a number", path, lineno, words[2]);
        return OSH_EPARSE;
    }
    for (i = 0; i < sizeof(field) - 1u && words[0][i]; ++i)
        field[i] = (char) toupper((unsigned char) words[0][i]);
    field[i] = '\0';
    return append_filter_rule(fil, field, words[1], value);
}
