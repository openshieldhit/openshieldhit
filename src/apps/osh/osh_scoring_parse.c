/**
 * @file osh_scoring_parse.c
 *
 * @brief Parser for the grouped detect.dat scoring configuration format.
 *
 * @details
 * The file format uses four section types, introduced by unindented keywords:
 *
 *   Filter
 *   Settings
 *   Geometry <type>
 *   Output
 *
 * Section content follows on indented lines.  Comments begin with '#' or '!'
 * and extend to the end of the line.  Blank lines are ignored everywhere.
 * Section order within the file is free; dangling name references are caught
 * by the post-parse validation step.
 *
 * Each section has its own dispatch table mapping lower-cased keyword strings
 * to handler functions.  Adding a new keyword means adding a key constant to
 * osh_scoring_parse_keys.h and a handler + table entry here — nothing else
 * changes in the main parse loop.
 *
 * Handler signature:
 *
 *   enum osh_status handler(section_ctx *ctx,
 *                           char **words, int nwords,
 *                           char const *path, unsigned int lineno)
 *
 * @p words[0] is already the (lower-cased) key; @p words[1..] are the
 * arguments.  Handlers receive tokenised words rather than a raw args string
 * because most scoring keywords take several space-separated values (axis
 * bounds, filter rules, etc.) that would need re-tokenising from a string.
 */

#include "apps/osh/osh_scoring_parse.h"

#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_scoring_parse_internal.h"
#include "common/osh_readline.h"
#include "openshieldhit/logger.h"
#include "openshieldhit/scoring.h"

/* ---- Internal section enum ----------------------------------------------- */

enum scoring_section { SECTION_NONE = 0, SECTION_FILTER, SECTION_SETTINGS, SECTION_GEOMETRY, SECTION_OUTPUT };

/* ---- Internal helper declarations ---------------------------------------- */

static enum osh_status append_filter(struct osh_scoring_workspace *ws);
static enum osh_status append_settings(struct osh_scoring_workspace *ws);
static enum osh_status append_geometry(struct osh_scoring_workspace *ws);
static enum osh_status append_output(struct osh_scoring_workspace *ws);
static enum osh_status validate(struct osh_scoring_workspace const *ws);

/* ---- Public entry point -------------------------------------------------- */

enum osh_status osh_scoring_parse(struct oshfile *oshf, struct osh_scoring_workspace *ws) {
    char const *path;
    enum scoring_section section;
    char *line;
    int lineno;
    enum osh_status rc;

    if (!oshf || !ws) {
        return OSH_EINVAL;
    }

    path = oshf->filename;
    section = SECTION_NONE;
    line = NULL;

    while (osh_readline(oshf, &line, &lineno) >= 0) {
        char *words[32];
        int nwords;
        int found;

        nwords = osh_tokenise(line, words, (int) (sizeof(words) / sizeof(words[0])));
        if (nwords <= 0)
            continue;

        osh_lower_inplace(words[0]);

        /* ---- Section headers ---------------------------------------------- */
        if (strcmp(words[0], "filter") == 0) {
            rc = append_filter(ws);
            if (rc != OSH_OK)
                goto fail;
            section = SECTION_FILTER;
            continue;
        }
        if (strcmp(words[0], "settings") == 0) {
            rc = append_settings(ws);
            if (rc != OSH_OK)
                goto fail;
            section = SECTION_SETTINGS;
            continue;
        }
        if (strcmp(words[0], "geometry") == 0) {
            rc = append_geometry(ws);
            if (rc != OSH_OK)
                goto fail;
            section = SECTION_GEOMETRY;
            if (nwords >= 2) {
                ws->geometries[ws->ngeometries - 1u].kind = strdup(words[1]);
                if (!ws->geometries[ws->ngeometries - 1u].kind) {
                    rc = OSH_ENOMEM;
                    goto fail;
                }
                osh_lower_inplace(ws->geometries[ws->ngeometries - 1u].kind);
            }
            continue;
        }
        if (strcmp(words[0], "output") == 0) {
            rc = append_output(ws);
            if (rc != OSH_OK)
                goto fail;
            section = SECTION_OUTPUT;
            continue;
        }

        /* ---- Section body: dispatch to the appropriate table -------------- */
        if (section == SECTION_NONE) {
            osh_error("%s:%u: keyword '%s' before any section header", path, lineno, words[0]);
            rc = OSH_EPARSE;
            goto fail;
        }

        found = 0;
        rc = OSH_OK;

        switch (section) {
        case SECTION_FILTER:
            rc = osh_scoring_parse_filter_line(
                &ws->filters[ws->nfilters - 1u], words, nwords, path, (unsigned int) lineno, &found);
            break;
        case SECTION_SETTINGS:
            rc = osh_scoring_parse_settings_line(
                &ws->settings[ws->nsettings - 1u], words, nwords, path, (unsigned int) lineno, &found);
            break;
        case SECTION_GEOMETRY:
            rc = osh_scoring_parse_geometry_line(
                &ws->geometries[ws->ngeometries - 1u], words, nwords, path, (unsigned int) lineno, &found);
            break;
        case SECTION_OUTPUT:
            rc = osh_scoring_parse_output_line(
                &ws->outputs[ws->noutputs - 1u], words, nwords, path, (unsigned int) lineno, &found);
            break;
        default:
            break;
        }

        if (rc != OSH_OK)
            goto fail;

        if (!found) {
            osh_warn("%s:%u: unknown key '%s' in section — ignored", path, lineno, words[0]);
        }
    }

    free(line);

    return validate(ws);

fail:
    free(line);
    return rc;
}

/* ---- Append helpers ------------------------------------------------------ */

/**
 * @brief Append a zero-initialized filter slot to the workspace.
 *
 * @param[in,out] ws  Scoring workspace; filters array is grown by one.
 *
 * @returns OSH_OK on success, OSH_ENOMEM on allocation failure.
 */
static enum osh_status append_filter(struct osh_scoring_workspace *ws) {
    struct osh_scoring_filter_def *tmp =
        (struct osh_scoring_filter_def *) realloc(ws->filters, (ws->nfilters + 1u) * sizeof(*tmp));
    if (!tmp)
        return OSH_ENOMEM;
    ws->filters = tmp;
    memset(&ws->filters[ws->nfilters], 0, sizeof(*tmp));
    ws->nfilters++;
    return OSH_OK;
}

/**
 * @brief Append a default-initialized settings slot to the workspace.
 *
 * @details Sets medium and nkmedium sentinel values to -1.
 *
 * @param[in,out] ws  Scoring workspace; settings array is grown by one.
 *
 * @returns OSH_OK on success, OSH_ENOMEM on allocation failure.
 */
static enum osh_status append_settings(struct osh_scoring_workspace *ws) {
    struct osh_scoring_settings_def *tmp =
        (struct osh_scoring_settings_def *) realloc(ws->settings, (ws->nsettings + 1u) * sizeof(*tmp));
    if (!tmp)
        return OSH_ENOMEM;
    ws->settings = tmp;
    memset(&ws->settings[ws->nsettings], 0, sizeof(*tmp));
    ws->settings[ws->nsettings].medium = -1;
    ws->settings[ws->nsettings].nkmedium = -1;
    ws->nsettings++;
    return OSH_OK;
}

/**
 * @brief Append a zero-initialized geometry slot to the workspace.
 *
 * @param[in,out] ws  Scoring workspace; geometries array is grown by one.
 *
 * @returns OSH_OK on success, OSH_ENOMEM on allocation failure.
 */
static enum osh_status append_geometry(struct osh_scoring_workspace *ws) {
    struct osh_scoring_geometry_def *tmp =
        (struct osh_scoring_geometry_def *) realloc(ws->geometries, (ws->ngeometries + 1u) * sizeof(*tmp));
    if (!tmp)
        return OSH_ENOMEM;
    ws->geometries = tmp;
    memset(&ws->geometries[ws->ngeometries], 0, sizeof(*tmp));
    ws->ngeometries++;
    return OSH_OK;
}

/**
 * @brief Append a zero-initialized output slot to the workspace.
 *
 * @param[in,out] ws  Scoring workspace; outputs array is grown by one.
 *
 * @returns OSH_OK on success, OSH_ENOMEM on allocation failure.
 */
static enum osh_status append_output(struct osh_scoring_workspace *ws) {
    struct osh_scoring_output_def *tmp =
        (struct osh_scoring_output_def *) realloc(ws->outputs, (ws->noutputs + 1u) * sizeof(*tmp));
    if (!tmp)
        return OSH_ENOMEM;
    ws->outputs = tmp;
    memset(&ws->outputs[ws->noutputs], 0, sizeof(*tmp));
    ws->noutputs++;
    return OSH_OK;
}

/* ---- Post-parse validation ----------------------------------------------- */

/**
 * @brief Validate the fully parsed scoring workspace for required fields.
 *
 * @details Checks that every filter has a Name, every settings block has a
 * Name, every geometry has both a type keyword and a Name, and every output
 * has a Filename, a Geo reference, and at least one Quantity page.
 *
 * @param[in] ws  Completed scoring workspace to validate.
 *
 * @returns OSH_OK if all required fields are present, OSH_EPARSE otherwise.
 */
static enum osh_status validate(struct osh_scoring_workspace const *ws) {
    size_t i, j;

    for (i = 0; i < ws->nfilters; ++i) {
        if (!ws->filters[i].name || ws->filters[i].name[0] == '\0') {
            osh_error("scoring: filter %zu is missing a Name", i);
            return OSH_EPARSE;
        }
    }

    for (i = 0; i < ws->nsettings; ++i) {
        if (!ws->settings[i].name || ws->settings[i].name[0] == '\0') {
            osh_error("scoring: settings %zu is missing a Name", i);
            return OSH_EPARSE;
        }
    }

    for (i = 0; i < ws->ngeometries; ++i) {
        if (!ws->geometries[i].kind || ws->geometries[i].kind[0] == '\0') {
            osh_error("scoring: geometry %zu is missing a type keyword", i);
            return OSH_EPARSE;
        }
        if (!ws->geometries[i].name || ws->geometries[i].name[0] == '\0') {
            osh_error("scoring: geometry %zu (%s) is missing a Name", i, ws->geometries[i].kind);
            return OSH_EPARSE;
        }
        if (ws->geometries[i].naxes == 0u && strcmp(ws->geometries[i].kind, "zone") != 0) {
            osh_warn("scoring: geometry '%s' has no axis definitions", ws->geometries[i].name);
        }
    }

    for (i = 0; i < ws->noutputs; ++i) {
        if (!ws->outputs[i].filename || ws->outputs[i].filename[0] == '\0') {
            osh_error("scoring: output %zu is missing a Filename", i);
            return OSH_EPARSE;
        }
        if (!ws->outputs[i].geometry_name || ws->outputs[i].geometry_name[0] == '\0') {
            osh_error("scoring: output '%s' is missing a Geo reference", ws->outputs[i].filename);
            return OSH_EPARSE;
        }
        if (ws->outputs[i].npages == 0u) {
            osh_error("scoring: output '%s' has no Quantity lines", ws->outputs[i].filename);
            return OSH_EPARSE;
        }
        for (j = 0; j < ws->outputs[i].npages; ++j) {
            if (!ws->outputs[i].pages[j].quantity || ws->outputs[i].pages[j].quantity[0] == '\0') {
                osh_error("scoring: output '%s' page %zu has an empty Quantity", ws->outputs[i].filename, j);
                return OSH_EPARSE;
            }
        }
    }

    return OSH_OK;
}
