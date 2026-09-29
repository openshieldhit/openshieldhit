/**
 * @file osh_scoring_multiformat.c
 *
 * @brief Compile-time expansion of multi-format scoring `Output` blocks (issue #308).
 *
 * @details
 * Extracted from osh_scoring_compile.c to keep that translation unit focused on
 * building the page/geometry runtime.  The single entry point,
 * @ref osh_scoring_expand_multiformat_outputs, is the former "Phase 7": it turns
 * one cold block requesting K formats into K runtime outputs that share the
 * block's page indices, derives per-format filenames, resolves the RTDOSE
 * target's single page, and rejects blocks whose targets would collide or misuse
 * RTDOSE.
 */

#include "scoring/runtime/osh_scoring_multiformat.h"

#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_diag.h"
#include "openshieldhit/scoring.h"
#include "scoring/runtime/osh_scoring_defs.h"
#include "scoring/runtime/osh_scoring_runtime.h"

static char const *format_canonical_ext(char const *fmt);
static int format_is_rtdose(char const *fmt);
static int score_kind_is_rtdose_compatible(enum osh_scoring_score_kind kind);
static int suffix_case_equal(char const *text, size_t textlen, char const *suffix, size_t suffixlen);
static char *derive_format_filename(char const *stem, char const *fmt);
static char const *output_display_name(struct osh_scoring_output_def const *out);
static enum osh_status
count_extra_outputs(struct osh_scoring_workspace const *ws, struct osh_diag_sink const *diag, size_t *n_extra_out);
static enum osh_status resolve_rtdose_page(struct osh_scoring_workspace const *ws,
                                           struct osh_scoring_runtime const *rt,
                                           struct osh_diag_sink const *diag,
                                           size_t block_idx,
                                           size_t nformats,
                                           size_t const *pages,
                                           size_t npages,
                                           size_t *dose_page_out);
static enum osh_status expand_one_block(struct osh_scoring_workspace const *ws,
                                        struct osh_diag_sink const *diag,
                                        struct osh_scoring_runtime *rt,
                                        size_t block_idx,
                                        size_t *next);
static enum osh_status check_block_collisions(struct osh_scoring_workspace const *ws,
                                              struct osh_scoring_runtime const *rt,
                                              struct osh_diag_sink const *diag,
                                              size_t block_idx,
                                              size_t block_first_extra,
                                              size_t nformats,
                                              char const *primary_name);

/* Canonical filename extension for a (lowercased) format keyword, or NULL when
 * the keyword is not a recognised writer.  Mirrors the save-layer dispatch in
 * osh_scoring_save.c so a derived multi-format filename matches the writer that
 * will consume it. */
static char const *format_canonical_ext(char const *fmt) {
    if (!fmt) {
        return NULL;
    }
    if (strcmp(fmt, "text") == 0 || strcmp(fmt, "txt") == 0 || strcmp(fmt, "ascii") == 0 || strcmp(fmt, "dat") == 0) {
        return ".dat";
    }
    if (strcmp(fmt, "bdo") == 0 || strcmp(fmt, "bdo2019") == 0 || strcmp(fmt, "binary") == 0
        || strcmp(fmt, "bin") == 0) {
        return ".bdo";
    }
    if (strcmp(fmt, "rtdose") == 0) {
        return ".dcm";
    }
    if (strcmp(fmt, "svg") == 0) {
        return ".svg";
    }
    return NULL;
}

static int format_is_rtdose(char const *fmt) {
    return fmt && strcmp(fmt, "rtdose") == 0;
}

static int score_kind_is_rtdose_compatible(enum osh_scoring_score_kind kind) {
    return kind == OSH_SCORING_SCORE_DOSE || kind == OSH_SCORING_SCORE_DOSEGY;
}

/* Return 1 when @p text ends with @p suffix (case-insensitive), else 0.
 * Comparison length is provided by @p textlen / @p suffixlen to avoid repeated
 * strlen() calls in derive_format_filename(). */
static int suffix_case_equal(char const *text, size_t textlen, char const *suffix, size_t suffixlen) {
    size_t i;
    unsigned char a;
    unsigned char b;

    if (textlen < suffixlen) {
        return 0;
    }
    for (i = 0u; i < suffixlen; ++i) {
        a = (unsigned char) text[textlen - suffixlen + i];
        b = (unsigned char) suffix[i];
        if (tolower(a) != tolower(b)) {
            return 0;
        }
    }
    return 1;
}

/* Derive a per-format filename from a stem: strip one recognised trailing
 * extension, then append the canonical extension for @p fmt.  Caller owns the
 * returned string; returns NULL on allocation failure (the format is validated
 * to be known before this is called). */
static char *derive_format_filename(char const *stem, char const *fmt) {
    static char const *const known_ext[] = {".dat", ".txt", ".bdo", ".bdz", ".bin", ".dcm", ".svg"};
    char const *ext = format_canonical_ext(fmt);
    size_t stemlen;
    size_t base;
    size_t extlen;
    size_t e;
    size_t xl;
    char *result;

    if (!stem || !ext) {
        return NULL;
    }
    stemlen = strlen(stem);
    base = stemlen;
    for (e = 0u; e < sizeof(known_ext) / sizeof(known_ext[0]); ++e) {
        xl = strlen(known_ext[e]);
        if (suffix_case_equal(stem, stemlen, known_ext[e], xl)) {
            base = stemlen - xl;
            break;
        }
    }
    extlen = strlen(ext);
    result = (char *) malloc(base + extlen + 1u);
    if (!result) {
        return NULL;
    }
    memcpy(result, stem, base);
    memcpy(result + base, ext, extlen + 1u);
    return result;
}

/* Label used only inside diagnostic message text to identify a block — its
 * Filename, or the placeholder "(unnamed)" when no Filename was given.  The
 * parenthesised placeholder is never used as (nor written to) a real output path;
 * the brackets make clear in a message that no name exists.  Matches the
 * convention already used across osh_scoring_compile.c. */
static char const *output_display_name(struct osh_scoring_output_def const *out) {
    if (out->filename) {
        return out->filename;
    }
    return "(unnamed)";
}

/* First pass: validate every requested format keyword and count how many extra
 * runtime outputs the expansion will add (one per additional format per block). */
static enum osh_status
count_extra_outputs(struct osh_scoring_workspace const *ws, struct osh_diag_sink const *diag, size_t *n_extra_out) {
    size_t i;
    size_t k;
    size_t n_extra;

    n_extra = 0u;
    for (i = 0; i < ws->noutputs; ++i) {
        size_t nf_i = ws->outputs[i].nfileformats;

        if (nf_i == 0u) {
            continue;
        }
        for (k = 0; k < nf_i; ++k) {
            char const *fmt = ws->outputs[i].fileformats[k];
            char const *fmt_name = "(null)";

            if (fmt != NULL) {
                fmt_name = fmt;
            }
            if (format_canonical_ext(fmt) == NULL) {
                OSH_DIAG_ERRORF(diag,
                                "Scoring output '%s' requests unknown FileFormat '%s'",
                                output_display_name(&ws->outputs[i]),
                                fmt_name);
                return OSH_ENOTSUP;
            }
        }
        if (nf_i > 1u) {
            n_extra += nf_i - 1u;
        }
    }
    *n_extra_out = n_extra;
    return OSH_OK;
}

/* Pick the single page an RTDOSE target in this block will encode into the .dcm.
 *
 * Rules (see docs/user/detect.dat.md):
 *   - a single-format RTDOSE block must score exactly one page (of any quantity);
 *     RTDOSE writes one grid, so more than one page is rejected;
 *   - a block with exactly one page uses that page (any quantity) — this is how
 *     `FileFormat RTDOSE` + `Quantity Energy` writes energy into the .dcm;
 *   - a block with several pages (necessarily a mixed block) uses the first
 *     Dose/DoseGy page; the other pages feed the block's other formats.  A mixed
 *     block with no Dose/DoseGy page is rejected.
 */
static enum osh_status resolve_rtdose_page(struct osh_scoring_workspace const *ws,
                                           struct osh_scoring_runtime const *rt,
                                           struct osh_diag_sink const *diag,
                                           size_t block_idx,
                                           size_t nformats,
                                           size_t const *pages,
                                           size_t npages,
                                           size_t *dose_page_out) {
    char const *name = output_display_name(&ws->outputs[block_idx]);
    size_t a;

    if (nformats == 1u && npages > 1u) {
        OSH_DIAG_ERRORF(diag,
                        "Scoring output '%s' requests RTDOSE alone but scores %zu pages; RTDOSE writes exactly "
                        "one page — request a single Quantity, or combine RTDOSE with another format",
                        name,
                        npages);
        return OSH_ENOTSUP;
    }
    if (npages == 1u) {
        *dose_page_out = pages[0];
        return OSH_OK;
    }
    for (a = 0u; a < npages; ++a) {
        if (score_kind_is_rtdose_compatible(rt->pages[pages[a]].score_kind)) {
            *dose_page_out = pages[a];
            return OSH_OK;
        }
    }
    OSH_DIAG_ERRORF(
        diag, "Scoring output '%s' combines RTDOSE with other formats but has no Dose/DoseGy page to write", name);
    return OSH_ENOTSUP;
}

/* Reject two targets in a block that resolve to the same output path (e.g.
 * "FileFormat TEXT DAT" — both canonicalise to .dat, or two overrides that
 * name the same file). */
static enum osh_status check_block_collisions(struct osh_scoring_workspace const *ws,
                                              struct osh_scoring_runtime const *rt,
                                              struct osh_diag_sink const *diag,
                                              size_t block_idx,
                                              size_t block_first_extra,
                                              size_t nformats,
                                              char const *primary_name) {
    size_t a;
    size_t b;

    for (a = 0u; a < nformats; ++a) {
        char const *na;

        if (a == 0u) {
            na = primary_name;
        } else {
            na = rt->outputs[block_first_extra + a - 1u].filename;
        }
        for (b = a + 1u; b < nformats; ++b) {
            char const *nb;

            if (b == 0u) {
                nb = primary_name;
            } else {
                nb = rt->outputs[block_first_extra + b - 1u].filename;
            }
            if (na != NULL && nb != NULL && strcmp(na, nb) == 0) {
                OSH_DIAG_ERRORF(diag,
                                "Scoring output '%s' requests formats that resolve to the same file '%s'",
                                output_display_name(&ws->outputs[block_idx]),
                                na);
                return OSH_EINVAL;
            }
        }
    }
    return OSH_OK;
}

/* Expand one cold block (index @p block_idx) into its extra runtime outputs,
 * starting at rt->outputs[*next].  Advances *next past the extras it appends. */
static enum osh_status expand_one_block(struct osh_scoring_workspace const *ws,
                                        struct osh_diag_sink const *diag,
                                        struct osh_scoring_runtime *rt,
                                        size_t block_idx,
                                        size_t *next) {
    struct osh_scoring_output_runtime *primary = &rt->outputs[block_idx];
    struct osh_scoring_output_def const *cold = &ws->outputs[block_idx];
    char const *stem = cold->filename;
    size_t nf_i = cold->nfileformats;
    size_t block_first_extra = *next;
    size_t primary_npages_full;
    size_t *primary_pages_full;
    size_t dose_page_index = 0u;
    int primary_is_rtdose;
    int has_primary_override;
    int block_has_rtdose = 0;
    char *resolved_name;
    enum osh_status rc;
    size_t k;

    /* A block with no FileFormat line keeps its default (BDO) runtime output from
     * phases 1-6 as is — nothing to expand, and fileformats[] is empty. */
    if (nf_i == 0u) {
        return OSH_OK;
    }

    primary_npages_full = primary->npages;
    primary_pages_full = primary->page_indices;
    primary_is_rtdose = format_is_rtdose(cold->fileformats[0]);
    has_primary_override = cold->fileformat_filenames != NULL && cold->fileformat_filenames[0] != NULL;

    for (k = 0u; k < nf_i; ++k) {
        if (format_is_rtdose(cold->fileformats[k])) {
            block_has_rtdose = 1;
            break;
        }
    }
    if (block_has_rtdose) {
        rc = resolve_rtdose_page(
            ws, rt, diag, block_idx, nf_i, primary_pages_full, primary_npages_full, &dose_page_index);
        if (rc != OSH_OK) {
            return rc;
        }
    }

    resolved_name = NULL;
    if (has_primary_override) {
        resolved_name = strdup(cold->fileformat_filenames[0]);
    } else if (nf_i > 1u) {
        resolved_name = derive_format_filename(stem, cold->fileformats[0]);
    }
    if (resolved_name != NULL) {
        free(primary->filename);
        primary->filename = resolved_name;
    } else if (has_primary_override || nf_i > 1u) {
        return OSH_ENOMEM;
    }

    for (k = 1u; k < nf_i; ++k) {
        struct osh_scoring_output_runtime *extra = &rt->outputs[*next];
        char const *fmt = cold->fileformats[k];
        int extra_is_rtdose;

        if (!fmt) {
            OSH_DIAG_ERRORF(diag, "Scoring output '%s' has a null FileFormat entry", output_display_name(cold));
            return OSH_EINVAL;
        }
        extra_is_rtdose = format_is_rtdose(fmt);
        extra->fileformat = strdup(fmt);
        if (!extra->fileformat) {
            return OSH_ENOMEM;
        }
        if (cold->fileformat_filenames != NULL && cold->fileformat_filenames[k] != NULL) {
            extra->filename = strdup(cold->fileformat_filenames[k]);
        } else {
            extra->filename = derive_format_filename(stem, fmt);
        }
        if (!extra->filename) {
            return OSH_ENOMEM;
        }
        extra->geometry_idx = primary->geometry_idx;
        extra->npages = primary_npages_full;
        if (extra_is_rtdose) {
            extra->npages = 1u;
        }
        if (extra->npages > 0u) {
            extra->page_indices = (size_t *) malloc(extra->npages * sizeof(*extra->page_indices));
            if (!extra->page_indices) {
                return OSH_ENOMEM;
            }
            if (extra_is_rtdose) {
                extra->page_indices[0] = dose_page_index;
            } else {
                memcpy(extra->page_indices, primary_pages_full, extra->npages * sizeof(*extra->page_indices));
            }
        }
        (*next)++;
    }

    if (primary_is_rtdose && primary_npages_full > 0u) {
        size_t *single_page = (size_t *) malloc(sizeof(*single_page));

        if (!single_page) {
            return OSH_ENOMEM;
        }
        single_page[0] = dose_page_index;
        free(primary->page_indices);
        primary->page_indices = single_page;
        primary->npages = 1u;
    }

    return check_block_collisions(ws, rt, diag, block_idx, block_first_extra, nf_i, primary->filename);
}

enum osh_status osh_scoring_expand_multiformat_outputs(struct osh_scoring_workspace const *ws,
                                                       struct osh_diag_sink const *diag,
                                                       struct osh_scoring_runtime *rt) {
    enum osh_status rc;
    size_t n_extra;
    size_t next;
    size_t i;

    if (!ws || !rt) {
        return OSH_EINVAL;
    }

    rc = count_extra_outputs(ws, diag, &n_extra);
    if (rc != OSH_OK) {
        return rc;
    }

    if (n_extra > 0u) {
        struct osh_scoring_output_runtime *grown =
            (struct osh_scoring_output_runtime *) realloc(rt->outputs, (rt->noutputs + n_extra) * sizeof(*rt->outputs));
        if (!grown) {
            return OSH_ENOMEM;
        }
        rt->outputs = grown;
        /* Zero the new tail so a mid-expansion failure leaves every slot either
         * fully built or cleanly zeroed; bump noutputs now so
         * osh_scoring_runtime_free() (run by the caller on error) covers them. */
        memset(&rt->outputs[rt->noutputs], 0, n_extra * sizeof(*rt->outputs));
        rt->noutputs += n_extra;
    }

    next = rt->noutputs - n_extra;
    for (i = 0; i < ws->noutputs; ++i) {
        rc = expand_one_block(ws, diag, rt, i, &next);
        if (rc != OSH_OK) {
            return rc;
        }
    }
    return OSH_OK;
}
