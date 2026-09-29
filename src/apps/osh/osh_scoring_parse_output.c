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

#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_scoring_parse_internal.h"
#include "apps/osh/osh_scoring_parse_keys.h"
#include "common/osh_diag.h"
#include "common/osh_readline.h"

typedef enum osh_status (*output_handler_fn)(
    struct osh_scoring_output_def *, struct osh_diag_sink const *, char **, int, char const *, unsigned int);

struct output_entry {
    char const *key;
    output_handler_fn handler;
};

static enum osh_status append_page(struct osh_scoring_output_def *out);
static enum osh_status append_page_filter(struct osh_scoring_page_def *page, char const *name);
static enum osh_status append_fileformat(struct osh_scoring_output_def *out, char const *keyword);
static enum osh_status set_fileformat_filename(struct osh_scoring_output_def *out, size_t idx, char const *name);
static int ascii_case_equal(char const *a, char const *b);
static int is_known_fileformat_keyword(char const *word);
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
static enum osh_status output_diff1(struct osh_scoring_output_def *out,
                                    struct osh_diag_sink const *diag,
                                    char **words,
                                    int nwords,
                                    char const *path,
                                    unsigned int lineno);
static enum osh_status output_diff1type(struct osh_scoring_output_def *out,
                                        struct osh_diag_sink const *diag,
                                        char **words,
                                        int nwords,
                                        char const *path,
                                        unsigned int lineno);
static enum osh_status output_diff2(struct osh_scoring_output_def *out,
                                    struct osh_diag_sink const *diag,
                                    char **words,
                                    int nwords,
                                    char const *path,
                                    unsigned int lineno);
static enum osh_status output_diff2type(struct osh_scoring_output_def *out,
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
                                             {"diff1", output_diff1},
                                             {"diff1type", output_diff1type},
                                             {"diff2", output_diff2},
                                             {"diff2type", output_diff2type},
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
 * @brief Append one lowercased format keyword to an output's format list.
 */
static enum osh_status append_fileformat(struct osh_scoring_output_def *out, char const *keyword) {
    char **tmp = (char **) realloc((void *) out->fileformats, (out->nfileformats + 1u) * sizeof(*tmp));
    char **tmp_names;

    if (!tmp) {
        return OSH_ENOMEM;
    }
    out->fileformats = tmp;
    tmp_names = (char **) realloc((void *) out->fileformat_filenames, (out->nfileformats + 1u) * sizeof(*tmp_names));
    if (!tmp_names) {
        return OSH_ENOMEM;
    }
    out->fileformat_filenames = tmp_names;
    out->fileformats[out->nfileformats] = strdup(keyword);
    if (!out->fileformats[out->nfileformats]) {
        return OSH_ENOMEM;
    }
    osh_lower_inplace(out->fileformats[out->nfileformats]);
    out->fileformat_filenames[out->nfileformats] = NULL;
    out->nfileformats++;
    return OSH_OK;
}

static enum osh_status set_fileformat_filename(struct osh_scoring_output_def *out, size_t idx, char const *name) {
    char *copy;

    if (idx >= out->nfileformats) {
        return OSH_EINVAL;
    }
    copy = strdup(name);
    if (!copy) {
        return OSH_ENOMEM;
    }
    free(out->fileformat_filenames[idx]);
    out->fileformat_filenames[idx] = copy;
    return OSH_OK;
}

static int ascii_case_equal(char const *a, char const *b) {
    size_t i;
    size_t na;
    size_t nb;
    unsigned char ca;
    unsigned char cb;

    if (a == NULL || b == NULL) {
        return 0;
    }
    na = strlen(a);
    nb = strlen(b);
    if (na != nb) {
        return 0;
    }
    for (i = 0u; i < na; ++i) {
        ca = (unsigned char) a[i];
        cb = (unsigned char) b[i];
        if ((char) tolower(ca) != (char) tolower(cb)) {
            return 0;
        }
    }
    return 1;
}

static int is_known_fileformat_keyword(char const *word) {
    if (word == NULL) {
        return 0;
    }
    if (ascii_case_equal(word, "text") || ascii_case_equal(word, "txt") || ascii_case_equal(word, "ascii")
        || ascii_case_equal(word, "dat")) {
        return 1;
    }
    if (ascii_case_equal(word, "bdo") || ascii_case_equal(word, "bdo2019") || ascii_case_equal(word, "binary")
        || ascii_case_equal(word, "bin")) {
        return 1;
    }
    if (ascii_case_equal(word, "rtdose") || ascii_case_equal(word, "svg")) {
        return 1;
    }
    return 0;
}

/**
 * @brief Parse `Fileformat <name> [<name> ...]` (and alias `Format ...`).
 *
 * @details
 * Accepts one or more format keywords on a single line — e.g.
 * `FileFormat TEXT BDO` requests the same scored page-set be written as both a
 * TEXT and a BDO file.  Each keyword is appended to the output's format list;
 * repeated `FileFormat` lines accumulate rather than overwrite.  The page
 * accumulators are compiled once regardless of how many formats are requested.
 */
static enum osh_status output_fileformat(struct osh_scoring_output_def *out,
                                         struct osh_diag_sink const *diag,
                                         char **words,
                                         int nwords,
                                         char const *path,
                                         unsigned int lineno) {
    enum osh_status rc;
    int i;
    size_t last_format_idx;
    int have_pending_format;
    size_t before_count;
    size_t idx;

    if (nwords < 2) {
        OSH_DIAG_ERRORF(diag, "%s:%u: Output Fileformat requires a format name", path, lineno);
        return OSH_EPARSE;
    }
    before_count = out->nfileformats;
    last_format_idx = 0u;
    have_pending_format = 0;
    for (i = 1; i < nwords; ++i) {
        if (is_known_fileformat_keyword(words[i])) {
            rc = append_fileformat(out, words[i]);
            if (rc != OSH_OK) {
                return rc;
            }
            last_format_idx = out->nfileformats - 1u;
            have_pending_format = 1;
            continue;
        }

        if (!have_pending_format) {
            OSH_DIAG_ERRORF(
                diag, "%s:%u: Output Fileformat override '%s' must follow a format keyword", path, lineno, words[i]);
            return OSH_EPARSE;
        }

        idx = last_format_idx;
        if (out->fileformat_filenames[idx] != NULL) {
            OSH_DIAG_ERRORF(diag,
                            "%s:%u: Output Fileformat format '%s' already has an override filename",
                            path,
                            lineno,
                            out->fileformats[idx]);
            return OSH_EPARSE;
        }
        rc = set_fileformat_filename(out, idx, words[i]);
        if (rc != OSH_OK) {
            return rc;
        }
        have_pending_format = 0;
    }

    if (out->nfileformats == before_count) {
        OSH_DIAG_ERRORF(diag, "%s:%u: Output Fileformat requires at least one known format keyword", path, lineno);
        return OSH_EPARSE;
    }
    return OSH_OK;
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

/**
 * @brief Parse `Diff1 <lo> <hi> <nbins> [LOG]`.
 *
 * Applies to the most recently added Quantity page.  Activates differential
 * scoring along a single axis whose type is set by the subsequent Diff1Type line.
 */
static enum osh_status output_diff1(struct osh_scoring_output_def *out,
                                    struct osh_diag_sink const *diag,
                                    char **words,
                                    int nwords,
                                    char const *path,
                                    unsigned int lineno) {
    struct osh_scoring_page_def *page;
    double lo;
    double hi;
    double nbins_d;
    char log_buf[8];

    if (out->npages == 0u) {
        OSH_DIAG_ERRORF(diag, "%s:%u: Diff1 must follow a Quantity line", path, lineno);
        return OSH_EPARSE;
    }
    if (nwords < 4) {
        OSH_DIAG_ERRORF(diag, "%s:%u: Diff1 requires lo hi nbins [LOG]", path, lineno);
        return OSH_EPARSE;
    }
    lo = strtod(words[1], NULL);
    hi = strtod(words[2], NULL);
    nbins_d = strtod(words[3], NULL);
    if (!(hi > lo) || !(nbins_d >= 1.0)) {
        OSH_DIAG_ERRORF(diag, "%s:%u: Diff1: invalid lo/hi/nbins", path, lineno);
        return OSH_EPARSE;
    }

    page = &out->pages[out->npages - 1u];
    page->diff_lo = lo;
    page->diff_hi = hi;
    page->diff_nbins = (size_t) nbins_d;
    page->diff_log = 0;

    if (nwords >= 5) {
        strncpy(log_buf, words[4], sizeof(log_buf) - 1u);
        log_buf[sizeof(log_buf) - 1u] = '\0';
        osh_lower_inplace(log_buf);
        if (strcmp(log_buf, "log") == 0) {
            if (!(lo > 0.0)) {
                OSH_DIAG_ERRORF(diag, "%s:%u: Diff1 LOG requires lo > 0", path, lineno);
                return OSH_EPARSE;
            }
            page->diff_log = 1;
        }
    }
    return OSH_OK;
}

/**
 * @brief Parse `Diff1Type <type>`.
 *
 * Applies to the most recently added Quantity page.  Sets the physical quantity
 * used for differential binning (ekin, let, qeff, enuc, eamu, or synonyms e/dedx).
 */
static enum osh_status output_diff1type(struct osh_scoring_output_def *out,
                                        struct osh_diag_sink const *diag,
                                        char **words,
                                        int nwords,
                                        char const *path,
                                        unsigned int lineno) {
    struct osh_scoring_page_def *page;
    char *kind;

    if (out->npages == 0u) {
        OSH_DIAG_ERRORF(diag, "%s:%u: Diff1Type must follow a Quantity line", path, lineno);
        return OSH_EPARSE;
    }
    if (nwords < 2) {
        OSH_DIAG_ERRORF(diag, "%s:%u: Diff1Type requires a type keyword", path, lineno);
        return OSH_EPARSE;
    }

    page = &out->pages[out->npages - 1u];
    free(page->diff_kind_str);
    kind = strdup(words[1]);
    if (!kind) {
        return OSH_ENOMEM;
    }
    osh_lower_inplace(kind);
    page->diff_kind_str = kind;

    /* Optional third word is a Settings name that overrides the stopping-power
     * medium used for LET/QEFF axis binning — e.g. "Diff1Type DEDX in_Si".
     * Always clear any previous name first so that "Diff1Type EKIN" (no third
     * word) does not inherit the override from an earlier Diff1Type call on
     * the same page. */
    free(page->diff_kind_sset_name);
    page->diff_kind_sset_name = NULL;
    if (nwords >= 3) {
        page->diff_kind_sset_name = strdup(words[2]);
        if (!page->diff_kind_sset_name) {
            return OSH_ENOMEM;
        }
    }
    return OSH_OK;
}

/**
 * @brief Parse `Diff2 <lo> <hi> <nbins> [LOG]`.
 *
 * Applies to the most recently added Quantity page.  Activates a second differential
 * axis; Diff1 must have been set on the same page first.
 */
static enum osh_status output_diff2(struct osh_scoring_output_def *out,
                                    struct osh_diag_sink const *diag,
                                    char **words,
                                    int nwords,
                                    char const *path,
                                    unsigned int lineno) {
    struct osh_scoring_page_def *page;
    double lo;
    double hi;
    double nbins_d;
    char log_buf[8];

    if (out->npages == 0u) {
        OSH_DIAG_ERRORF(diag, "%s:%u: Diff2 must follow a Quantity line", path, lineno);
        return OSH_EPARSE;
    }
    if (nwords < 4) {
        OSH_DIAG_ERRORF(diag, "%s:%u: Diff2 requires lo hi nbins [LOG]", path, lineno);
        return OSH_EPARSE;
    }
    if (out->pages[out->npages - 1u].diff_nbins == 0u) {
        OSH_DIAG_ERRORF(diag, "%s:%u: Diff2 requires Diff1 to be set first on this Quantity", path, lineno);
        return OSH_EPARSE;
    }
    lo = strtod(words[1], NULL);
    hi = strtod(words[2], NULL);
    nbins_d = strtod(words[3], NULL);
    if (!(hi > lo) || !(nbins_d >= 1.0)) {
        OSH_DIAG_ERRORF(diag, "%s:%u: Diff2: invalid lo/hi/nbins", path, lineno);
        return OSH_EPARSE;
    }

    page = &out->pages[out->npages - 1u];
    page->diff2_lo = lo;
    page->diff2_hi = hi;
    page->diff2_nbins = (size_t) nbins_d;
    page->diff2_log = 0;

    if (nwords >= 5) {
        strncpy(log_buf, words[4], sizeof(log_buf) - 1u);
        log_buf[sizeof(log_buf) - 1u] = '\0';
        osh_lower_inplace(log_buf);
        if (strcmp(log_buf, "log") == 0) {
            if (!(lo > 0.0)) {
                OSH_DIAG_ERRORF(diag, "%s:%u: Diff2 LOG requires lo > 0", path, lineno);
                return OSH_EPARSE;
            }
            page->diff2_log = 1;
        }
    }
    return OSH_OK;
}

/**
 * @brief Parse `Diff2Type <type>`.
 *
 * Applies to the most recently added Quantity page.
 */
static enum osh_status output_diff2type(struct osh_scoring_output_def *out,
                                        struct osh_diag_sink const *diag,
                                        char **words,
                                        int nwords,
                                        char const *path,
                                        unsigned int lineno) {
    struct osh_scoring_page_def *page;
    char *kind;

    if (out->npages == 0u) {
        OSH_DIAG_ERRORF(diag, "%s:%u: Diff2Type must follow a Quantity line", path, lineno);
        return OSH_EPARSE;
    }
    if (nwords < 2) {
        OSH_DIAG_ERRORF(diag, "%s:%u: Diff2Type requires a type keyword", path, lineno);
        return OSH_EPARSE;
    }

    page = &out->pages[out->npages - 1u];
    free(page->diff2_kind_str);
    kind = strdup(words[1]);
    if (!kind) {
        return OSH_ENOMEM;
    }
    osh_lower_inplace(kind);
    page->diff2_kind_str = kind;

    /* Optional third word is a Settings name for the diff2 axis SP override
     * — e.g. "Diff2Type LET in_Water".  Always clear first so that a bare
     * "Diff2Type EKIN" does not inherit a previous override on the same page. */
    free(page->diff2_kind_sset_name);
    page->diff2_kind_sset_name = NULL;
    if (nwords >= 3) {
        page->diff2_kind_sset_name = strdup(words[2]);
        if (!page->diff2_kind_sset_name) {
            return OSH_ENOMEM;
        }
    }
    return OSH_OK;
}
