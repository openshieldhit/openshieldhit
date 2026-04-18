/**
 * @file osh_scoring_parse_output.c
 *
 * @brief Parse one tokenized line inside a scoring `Output` section.
 *
 * @details
 * Recognized keys:
 * - `Filename <path>`
 * - `Geo <name>`
 * - `Fileformat|Format <name>`
 * - `Quantity <name> [filter_name ...]`
 *
 * Each `Quantity` line creates one page entry and records optional filter-name
 * references for late resolution.
 */

#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_scoring_parse_internal.h"
#include "apps/osh/osh_scoring_parse_keys.h"
#include "common/osh_logger.h"
#include "common/osh_readline.h"

typedef enum osh_status (*output_handler_fn)(
    struct osh_scoring_output_def *, struct osh_diag_sink const *, char **, int, char const *, unsigned int);

struct output_entry {
    char const *key;
    output_handler_fn handler;
};

static enum osh_status append_page(struct osh_scoring_output_def *out);
static enum osh_status append_page_filter(struct osh_scoring_page_def *page, char const *name);
static enum osh_status output_filename(struct osh_scoring_output_def *out,
                                       struct osh_diag_sink const *diag,
                                       char **words,
                                       int nwords,
                                       char const *path,
                                       unsigned int lineno);
static enum osh_status output_geo(struct osh_scoring_output_def *out,
                                  struct osh_diag_sink const *diag,
                                  char **words,
                                  int nwords,
                                  char const *path,
                                  unsigned int lineno);
static enum osh_status output_fileformat(struct osh_scoring_output_def *out,
                                         struct osh_diag_sink const *diag,
                                         char **words,
                                         int nwords,
                                         char const *path,
                                         unsigned int lineno);
static enum osh_status output_quantity(struct osh_scoring_output_def *out,
                                       struct osh_diag_sink const *diag,
                                       char **words,
                                       int nwords,
                                       char const *path,
                                       unsigned int lineno);

static struct output_entry output_table[] = {{OSH_SCORING_KEY_FILENAME, output_filename},
                                             {OSH_SCORING_KEY_GEO_REF, output_geo},
                                             {OSH_SCORING_KEY_FILEFORMAT, output_fileformat},
                                             {"format", output_fileformat},
                                             {OSH_SCORING_KEY_QUANTITY, output_quantity},
                                             {NULL, NULL}};

/**
 * @brief Dispatch one tokenized line into the active output definition.
 */
enum osh_status osh_scoring_parse_output_line(struct osh_scoring_output_def *out,
                                              struct osh_diag_sink const *diag,
                                              char **words,
                                              int nwords,
                                              char const *path,
                                              unsigned int lineno,
                                              int *found_out) {
    size_t i;

    if (found_out)
        *found_out = 0;
    for (i = 0; output_table[i].key != NULL; ++i) {
        if (strcmp(output_table[i].key, words[0]) == 0) {
            if (found_out)
                *found_out = 1;
            return output_table[i].handler(out, diag, words, nwords, path, lineno);
        }
    }
    return OSH_OK;
}

/**
 * @brief Append one scoring page entry to an output.
 */
static enum osh_status append_page(struct osh_scoring_output_def *out) {
    struct osh_scoring_page_def *tmp =
        (struct osh_scoring_page_def *) realloc(out->pages, (out->npages + 1u) * sizeof(*tmp));
    if (!tmp)
        return OSH_ENOMEM;
    out->pages = tmp;
    memset(&out->pages[out->npages], 0, sizeof(*tmp));
    out->npages++;
    return OSH_OK;
}

/**
 * @brief Append one filter name reference to a page.
 */
static enum osh_status append_page_filter(struct osh_scoring_page_def *page, char const *name) {
    char **tmp = (char **) realloc(page->filter_names, (page->nfilter_names + 1u) * sizeof(*tmp));
    if (!tmp)
        return OSH_ENOMEM;
    page->filter_names = tmp;
    page->filter_names[page->nfilter_names] = strdup(name);
    if (!page->filter_names[page->nfilter_names])
        return OSH_ENOMEM;
    page->nfilter_names++;
    return OSH_OK;
}

/**
 * @brief Parse `Filename <value>`.
 */
static enum osh_status output_filename(struct osh_scoring_output_def *out,
                                       struct osh_diag_sink const *diag,
                                       char **words,
                                       int nwords,
                                       char const *path,
                                       unsigned int lineno) {
    if (nwords < 2) {
        OSH_DIAG_ERRORF(diag, "%s:%u: Output Filename requires an argument", path, lineno);
        return OSH_EPARSE;
    }
    free(out->filename);
    out->filename = strdup(words[1]);
    return out->filename ? OSH_OK : OSH_ENOMEM;
}

/**
 * @brief Parse `Geo <geometry_name>`.
 */
static enum osh_status output_geo(struct osh_scoring_output_def *out,
                                  struct osh_diag_sink const *diag,
                                  char **words,
                                  int nwords,
                                  char const *path,
                                  unsigned int lineno) {
    if (nwords < 2) {
        OSH_DIAG_ERRORF(diag, "%s:%u: Output Geo requires a geometry name", path, lineno);
        return OSH_EPARSE;
    }
    free(out->geometry_name);
    out->geometry_name = strdup(words[1]);
    return out->geometry_name ? OSH_OK : OSH_ENOMEM;
}

/**
 * @brief Parse `Fileformat <name>` (and alias `Format <name>`).
 */
static enum osh_status output_fileformat(struct osh_scoring_output_def *out,
                                         struct osh_diag_sink const *diag,
                                         char **words,
                                         int nwords,
                                         char const *path,
                                         unsigned int lineno) {
    if (nwords < 2) {
        OSH_DIAG_ERRORF(diag, "%s:%u: Output Fileformat requires a format name", path, lineno);
        return OSH_EPARSE;
    }
    free(out->fileformat);
    out->fileformat = strdup(words[1]);
    osh_lower_inplace(out->fileformat);
    return out->fileformat ? OSH_OK : OSH_ENOMEM;
}

/**
 * @brief Parse one `Quantity` line and append the corresponding page.
 */
static enum osh_status output_quantity(struct osh_scoring_output_def *out,
                                       struct osh_diag_sink const *diag,
                                       char **words,
                                       int nwords,
                                       char const *path,
                                       unsigned int lineno) {
    enum osh_status rc;
    size_t page_idx;
    int i;

    if (nwords < 2) {
        OSH_DIAG_ERRORF(diag, "%s:%u: Output Quantity requires a quantity name", path, lineno);
        return OSH_EPARSE;
    }
    rc = append_page(out);
    if (rc != OSH_OK)
        return rc;
    page_idx = out->npages - 1u;

    out->pages[page_idx].quantity = strdup(words[1]);
    if (!out->pages[page_idx].quantity)
        return OSH_ENOMEM;
    osh_lower_inplace(out->pages[page_idx].quantity);

    for (i = 2; i < nwords; ++i) {
        rc = append_page_filter(&out->pages[page_idx], words[i]);
        if (rc != OSH_OK)
            return rc;
    }
    return OSH_OK;
}
