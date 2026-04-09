#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_logger.h"
#include "scoring/parse/osh_scoring_parse_internal.h"
#include "scoring/parse/osh_scoring_parse_keys.h"

typedef enum osh_status (*filter_handler_fn)(struct osh_scoring_filter_def *, char **, int, char const *, unsigned int);

struct filter_entry {
    char const *key;
    filter_handler_fn handler;
};

static enum osh_status
append_filter_rule(struct osh_scoring_filter_def *fil, char const *field, char const *op, double value);
static enum osh_status
filter_name(struct osh_scoring_filter_def *fil, char **words, int nwords, char const *path, unsigned int lineno);
static enum osh_status
filter_rule(struct osh_scoring_filter_def *fil, char **words, int nwords, char const *path, unsigned int lineno);

static struct filter_entry filter_table[] = {{OSH_SCORING_KEY_NAME, filter_name},
                                             {OSH_SCORING_KEY_FILTER_Z, filter_rule},
                                             {OSH_SCORING_KEY_FILTER_A, filter_rule},
                                             {OSH_SCORING_KEY_FILTER_E, filter_rule},
                                             {OSH_SCORING_KEY_FILTER_GEN, filter_rule},
                                             {OSH_SCORING_KEY_FILTER_ID, filter_rule},
                                             {NULL, NULL}};

enum osh_status osh_scoring_parse_filter_line(struct osh_scoring_filter_def *fil,
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
            return filter_table[i].handler(fil, words, nwords, path, lineno);
        }
    }
    return OSH_OK;
}

static enum osh_status
append_filter_rule(struct osh_scoring_filter_def *fil, char const *field, char const *op, double value) {
    struct osh_scoring_filter_rule *tmp =
        (struct osh_scoring_filter_rule *) realloc(fil->rules, (fil->nrules + 1u) * sizeof(*tmp));
    if (!tmp)
        return OSH_ENOMEM;
    fil->rules = tmp;
    memset(&fil->rules[fil->nrules], 0, sizeof(*tmp));
    strncpy(fil->rules[fil->nrules].field, field, sizeof(fil->rules[0].field) - 1u);
    strncpy(fil->rules[fil->nrules].op, op, sizeof(fil->rules[0].op) - 1u);
    fil->rules[fil->nrules].value = value;
    fil->nrules++;
    return OSH_OK;
}

static enum osh_status
filter_name(struct osh_scoring_filter_def *fil, char **words, int nwords, char const *path, unsigned int lineno) {
    if (nwords < 2) {
        osh_error("%s:%u: Filter Name requires an argument", path, lineno);
        return OSH_EPARSE;
    }
    free(fil->name);
    fil->name = strdup(words[1]);
    return fil->name ? OSH_OK : OSH_ENOMEM;
}

static enum osh_status
filter_rule(struct osh_scoring_filter_def *fil, char **words, int nwords, char const *path, unsigned int lineno) {
    char field[16];
    char *endptr;
    double value;
    size_t i;

    if (nwords < 3) {
        osh_error("%s:%u: filter rule requires: <field> <op> <value>", path, lineno);
        return OSH_EPARSE;
    }
    value = strtod(words[2], &endptr);
    if (endptr == words[2]) {
        osh_error("%s:%u: filter rule value '%s' is not a number", path, lineno, words[2]);
        return OSH_EPARSE;
    }
    for (i = 0; i < sizeof(field) - 1u && words[0][i]; ++i)
        field[i] = (char) toupper((unsigned char) words[0][i]);
    field[i] = '\0';
    return append_filter_rule(fil, field, words[1], value);
}
