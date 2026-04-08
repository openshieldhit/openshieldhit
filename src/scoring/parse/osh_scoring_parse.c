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
 * Section order within the file is free (single-pass; dangling name references
 * are caught by the post-parse validation step).
 *
 * The parser produces raw osh_scoring_workspace with all fields filled from
 * the file.  No name resolution, volume computation, or runtime allocation
 * happens here.
 */

#include "scoring/parse/osh_scoring_parse.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_logger.h"
#include "common/osh_readline.h"
#include "scoring/osh_scoring.h"

/* ---- Internal section enum ----------------------------------------------- */

enum scoring_section {
    SECTION_NONE = 0,
    SECTION_FILTER,
    SECTION_SETTINGS,
    SECTION_GEOMETRY,
    SECTION_OUTPUT
};

/* ---- String helpers ------------------------------------------------------- */

static char *dupstr(char const *src) {
    size_t n;
    char *dst;
    if (!src) return NULL;
    n   = strlen(src);
    dst = (char *) malloc(n + 1u);
    if (!dst) return NULL;
    memcpy(dst, src, n + 1u);
    return dst;
}

/** Strip inline '#' or '!' comments in-place by replacing the marker with NUL. */
static void strip_comment(char *s) {
    for (; s && *s; ++s) {
        if (*s == '#' || *s == '!') {
            *s = '\0';
            return;
        }
    }
}

/** Advance past leading whitespace; trim trailing whitespace in-place.
 *  Returns pointer into the original buffer (never allocates). */
static char *trim(char *s) {
    char *end;
    if (!s) return NULL;
    while (*s && isspace((unsigned char)*s)) ++s;
    if (*s == '\0') return s;
    end = s + strlen(s);
    while (end > s && isspace((unsigned char) end[-1])) --end;
    *end = '\0';
    return s;
}

static int is_indented(char const *raw_line) {
    return raw_line && (*raw_line == ' ' || *raw_line == '\t');
}


/* ---- Append helpers ------------------------------------------------------- */

static enum osh_status append_filter(struct osh_scoring_workspace *ws) {
    struct osh_scoring_filter_def *tmp;
    tmp = (struct osh_scoring_filter_def *) realloc(ws->filters,
                                                     (ws->nfilters + 1u) * sizeof(*tmp));
    if (!tmp) return OSH_ENOMEM;
    ws->filters = tmp;
    memset(&ws->filters[ws->nfilters], 0, sizeof(*tmp));
    ws->nfilters++;
    return OSH_OK;
}

static enum osh_status append_settings(struct osh_scoring_workspace *ws) {
    struct osh_scoring_settings_def *tmp;
    tmp = (struct osh_scoring_settings_def *) realloc(ws->settings,
                                                       (ws->nsettings + 1u) * sizeof(*tmp));
    if (!tmp) return OSH_ENOMEM;
    ws->settings = tmp;
    memset(&ws->settings[ws->nsettings], 0, sizeof(*tmp));
    ws->settings[ws->nsettings].medium   = -1;
    ws->settings[ws->nsettings].nkmedium = -1;
    ws->nsettings++;
    return OSH_OK;
}

static enum osh_status append_geometry(struct osh_scoring_workspace *ws) {
    struct osh_scoring_geometry_def *tmp;
    tmp = (struct osh_scoring_geometry_def *) realloc(ws->geometries,
                                                       (ws->ngeometries + 1u) * sizeof(*tmp));
    if (!tmp) return OSH_ENOMEM;
    ws->geometries = tmp;
    memset(&ws->geometries[ws->ngeometries], 0, sizeof(*tmp));
    ws->ngeometries++;
    return OSH_OK;
}

static enum osh_status append_output(struct osh_scoring_workspace *ws) {
    struct osh_scoring_output_def *tmp;
    tmp = (struct osh_scoring_output_def *) realloc(ws->outputs,
                                                     (ws->noutputs + 1u) * sizeof(*tmp));
    if (!tmp) return OSH_ENOMEM;
    ws->outputs = tmp;
    memset(&ws->outputs[ws->noutputs], 0, sizeof(*tmp));
    ws->noutputs++;
    return OSH_OK;
}

static enum osh_status append_page(struct osh_scoring_output_def *out) {
    struct osh_scoring_page_def *tmp;
    tmp = (struct osh_scoring_page_def *) realloc(out->pages,
                                                   (out->npages + 1u) * sizeof(*tmp));
    if (!tmp) return OSH_ENOMEM;
    out->pages = tmp;
    memset(&out->pages[out->npages], 0, sizeof(*tmp));
    out->npages++;
    return OSH_OK;
}

static enum osh_status append_page_filter(struct osh_scoring_page_def *page, char const *name) {
    char **tmp;
    tmp = (char **) realloc(page->filter_names,
                             (page->nfilter_names + 1u) * sizeof(*tmp));
    if (!tmp) return OSH_ENOMEM;
    page->filter_names = tmp;
    page->filter_names[page->nfilter_names] = dupstr(name);
    if (!page->filter_names[page->nfilter_names]) return OSH_ENOMEM;
    page->nfilter_names++;
    return OSH_OK;
}

static enum osh_status append_filter_rule(struct osh_scoring_filter_def *fil,
                                           char const *field, char const *op, double value) {
    struct osh_scoring_filter_rule *tmp;
    tmp = (struct osh_scoring_filter_rule *) realloc(fil->rules,
                                                      (fil->nrules + 1u) * sizeof(*tmp));
    if (!tmp) return OSH_ENOMEM;
    fil->rules = tmp;
    memset(&fil->rules[fil->nrules], 0, sizeof(*tmp));
    /* Truncate safely — field[16] and op[4] are sized for the known keywords. */
    strncpy(fil->rules[fil->nrules].field, field, sizeof(fil->rules[0].field) - 1u);
    strncpy(fil->rules[fil->nrules].op,    op,    sizeof(fil->rules[0].op)    - 1u);
    fil->rules[fil->nrules].value = value;
    fil->nrules++;
    return OSH_OK;
}

static enum osh_status append_axis(struct osh_scoring_geometry_def *geo,
                                    char const *label, double lo, double hi, int nbins) {
    struct osh_scoring_axis_def *tmp;
    tmp = (struct osh_scoring_axis_def *) realloc(geo->axes,
                                                   (geo->naxes + 1u) * sizeof(*tmp));
    if (!tmp) return OSH_ENOMEM;
    geo->axes = tmp;
    memset(&geo->axes[geo->naxes], 0, sizeof(*tmp));
    strncpy(geo->axes[geo->naxes].label, label, sizeof(geo->axes[0].label) - 1u);
    geo->axes[geo->naxes].lo    = lo;
    geo->axes[geo->naxes].hi    = hi;
    geo->axes[geo->naxes].nbins = nbins;
    geo->naxes++;
    return OSH_OK;
}

/* ---- Section-level line handlers ----------------------------------------- */

static enum osh_status handle_filter_line(struct osh_scoring_workspace *ws,
                                           char **words, int nwords,
                                           char const *path, unsigned int lineno) {
    struct osh_scoring_filter_def *fil = &ws->filters[ws->nfilters - 1u];

    if (nwords < 1) return OSH_OK;

    if (strcasecmp(words[0], "Name") == 0) {
        if (nwords < 2) {
            osh_error("%s:%u: Filter Name requires an argument", path, lineno);
            return OSH_EPARSE;
        }
        free(fil->name);
        fil->name = dupstr(words[1]);
        return fil->name ? OSH_OK : OSH_ENOMEM;
    }

    /*
     * Filter rule: <field> <op> <value>
     * e.g.  Z = 6   or   E > 0.1   or   A <= 12
     * The operator may be one of: =  !=  <  >  <=  >=
     */
    if (nwords >= 3) {
        char *endptr;
        double value = strtod(words[2], &endptr);
        if (endptr == words[2]) {
            osh_error("%s:%u: filter rule value '%s' is not a number", path, lineno, words[2]);
            return OSH_EPARSE;
        }
        return append_filter_rule(fil, words[0], words[1], value);
    }

    osh_warn("%s:%u: unrecognised Filter line ignored: '%s'", path, lineno, words[0]);
    return OSH_OK;
}

static enum osh_status handle_settings_line(struct osh_scoring_workspace *ws,
                                              char **words, int nwords,
                                              char const *path, unsigned int lineno) {
    struct osh_scoring_settings_def *set = &ws->settings[ws->nsettings - 1u];

    if (nwords < 2) return OSH_OK;

    if (strcasecmp(words[0], "Name") == 0) {
        free(set->name);
        set->name = dupstr(words[1]);
        return set->name ? OSH_OK : OSH_ENOMEM;
    }
    if (strcasecmp(words[0], "Rescale") == 0) {
        set->rescale     = atof(words[1]);
        set->has_rescale = 1u;
        return OSH_OK;
    }
    if (strcasecmp(words[0], "Offset") == 0) {
        set->offset     = atof(words[1]);
        set->has_offset = 1u;
        return OSH_OK;
    }
    if (strcasecmp(words[0], "Medium") == 0) {
        set->medium     = atoi(words[1]);
        set->has_medium = 1u;
        return OSH_OK;
    }
    if (strcasecmp(words[0], "NKMedium") == 0) {
        set->nkmedium     = atoi(words[1]);
        set->has_nkmedium = 1u;
        return OSH_OK;
    }
    if (strcasecmp(words[0], "SiteDiameter") == 0 || strcasecmp(words[0], "SiteDiam") == 0) {
        set->site_diameter_um     = atof(words[1]);
        set->has_site_diameter_um = 1u;
        return OSH_OK;
    }
    if (strcasecmp(words[0], "Density") == 0 || strcasecmp(words[0], "RHO") == 0) {
        set->density_g_cm3     = atof(words[1]);
        set->has_density_g_cm3 = 1u;
        return OSH_OK;
    }
    if (strcasecmp(words[0], "MaxCount") == 0 || strcasecmp(words[0], "Npart") == 0) {
        set->npart     = (size_t) strtoull(words[1], NULL, 10);
        set->has_npart = 1u;
        return OSH_OK;
    }

    osh_warn("%s:%u: unrecognised Settings key ignored: '%s'", path, lineno, words[0]);
    return OSH_OK;
}

static enum osh_status handle_geometry_line(struct osh_scoring_workspace *ws,
                                              char **words, int nwords,
                                              char const *path, unsigned int lineno) {
    struct osh_scoring_geometry_def *geo = &ws->geometries[ws->ngeometries - 1u];

    if (nwords < 1) return OSH_OK;

    if (strcasecmp(words[0], "Name") == 0) {
        if (nwords < 2) {
            osh_error("%s:%u: Geometry Name requires an argument", path, lineno);
            return OSH_EPARSE;
        }
        free(geo->name);
        geo->name = dupstr(words[1]);
        return geo->name ? OSH_OK : OSH_ENOMEM;
    }

    /* Axis lines: X lo hi nbins  |  Y lo hi nbins  |  Z lo hi nbins  |  R lo hi nbins */
    if (strcasecmp(words[0], "X") == 0 || strcasecmp(words[0], "Y") == 0 ||
        strcasecmp(words[0], "Z") == 0 || strcasecmp(words[0], "R") == 0) {
        double lo, hi;
        int nbins;
        if (nwords < 4) {
            osh_error("%s:%u: Geometry axis '%s' requires lo hi nbins", path, lineno, words[0]);
            return OSH_EPARSE;
        }
        lo    = atof(words[1]);
        hi    = atof(words[2]);
        nbins = atoi(words[3]);
        /* Uppercase the label for canonical storage. */
        char label[2] = { (char) toupper((unsigned char) words[0][0]), '\0' };
        return append_axis(geo, label, lo, hi, nbins);
    }

    if (strcasecmp(words[0], "Rotation") == 0 || strcasecmp(words[0], "Rot") == 0) {
        if (nwords < 3) {
            osh_error("%s:%u: Geometry Rotation requires theta phi [deg]", path, lineno);
            return OSH_EPARSE;
        }
        geo->rot_theta_deg = atof(words[1]);
        geo->rot_phi_deg   = atof(words[2]);
        geo->has_rotation  = 1u;
        return OSH_OK;
    }

    if (strcasecmp(words[0], "Zones") == 0) {
        if (nwords < 3) {
            osh_error("%s:%u: Geometry Zones requires start stop", path, lineno);
            return OSH_EPARSE;
        }
        geo->zone_start = atoi(words[1]);
        geo->zone_stop  = atoi(words[2]);
        return OSH_OK;
    }

    osh_warn("%s:%u: unrecognised Geometry key ignored: '%s'", path, lineno, words[0]);
    return OSH_OK;
}

static enum osh_status handle_output_line(struct osh_scoring_workspace *ws,
                                           char **words, int nwords,
                                           char const *path, unsigned int lineno) {
    struct osh_scoring_output_def *out = &ws->outputs[ws->noutputs - 1u];
    int i;

    if (nwords < 2) return OSH_OK;

    if (strcasecmp(words[0], "Filename") == 0) {
        free(out->filename);
        out->filename = dupstr(words[1]);
        return out->filename ? OSH_OK : OSH_ENOMEM;
    }
    if (strcasecmp(words[0], "Geo") == 0 || strcasecmp(words[0], "Geometry") == 0) {
        free(out->geometry_name);
        out->geometry_name = dupstr(words[1]);
        return out->geometry_name ? OSH_OK : OSH_ENOMEM;
    }
    if (strcasecmp(words[0], "Fileformat") == 0 || strcasecmp(words[0], "Format") == 0) {
        free(out->fileformat);
        out->fileformat = dupstr(words[1]);
        return out->fileformat ? OSH_OK : OSH_ENOMEM;
    }
    if (strcasecmp(words[0], "Quantity") == 0) {
        enum osh_status rc;
        size_t page_idx;

        rc = append_page(out);
        if (rc != OSH_OK) return rc;
        page_idx = out->npages - 1u;

        out->pages[page_idx].quantity = dupstr(words[1]);
        if (!out->pages[page_idx].quantity) return OSH_ENOMEM;

        /* Remaining words on the Quantity line are filter names. */
        for (i = 2; i < nwords; ++i) {
            rc = append_page_filter(&out->pages[page_idx], words[i]);
            if (rc != OSH_OK) return rc;
        }
        return OSH_OK;
    }

    osh_warn("%s:%u: unrecognised Output key ignored: '%s'", path, lineno, words[0]);
    return OSH_OK;
}

/* ---- Post-parse validation ------------------------------------------------ */

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
            osh_error("scoring: geometry %zu (%s) is missing a Name",
                      i, ws->geometries[i].kind);
            return OSH_EPARSE;
        }
        if (ws->geometries[i].naxes == 0u &&
            strcasecmp(ws->geometries[i].kind, "Zone") != 0) {
            osh_warn("scoring: geometry '%s' has no axis definitions",
                     ws->geometries[i].name);
        }
    }

    for (i = 0; i < ws->noutputs; ++i) {
        if (!ws->outputs[i].filename || ws->outputs[i].filename[0] == '\0') {
            osh_error("scoring: output %zu is missing a Filename", i);
            return OSH_EPARSE;
        }
        if (!ws->outputs[i].geometry_name || ws->outputs[i].geometry_name[0] == '\0') {
            osh_error("scoring: output '%s' is missing a Geo reference",
                      ws->outputs[i].filename);
            return OSH_EPARSE;
        }
        if (ws->outputs[i].npages == 0u) {
            osh_error("scoring: output '%s' has no Quantity lines",
                      ws->outputs[i].filename);
            return OSH_EPARSE;
        }
        for (j = 0; j < ws->outputs[i].npages; ++j) {
            if (!ws->outputs[i].pages[j].quantity ||
                ws->outputs[i].pages[j].quantity[0] == '\0') {
                osh_error("scoring: output '%s' page %zu has an empty Quantity",
                          ws->outputs[i].filename, j);
                return OSH_EPARSE;
            }
        }
    }

    return OSH_OK;
}

/* ---- Public entry point --------------------------------------------------- */

enum osh_status osh_scoring_parse_file(char const *path, struct osh_scoring_workspace **ws_out) {
    FILE *fp;
    struct osh_scoring_workspace *ws;
    enum scoring_section section;
    char line[4096];
    unsigned int lineno;
    enum osh_status rc;

    if (!path || !ws_out) return OSH_EINVAL;
    *ws_out = NULL;

    fp = fopen(path, "r");
    if (!fp) {
        osh_error("scoring: cannot open '%s'", path);
        return OSH_EIO;
    }

    ws = (struct osh_scoring_workspace *) calloc(1, sizeof(*ws));
    if (!ws) { fclose(fp); return OSH_ENOMEM; }

    ws->fname = dupstr(path);
    if (!ws->fname) { fclose(fp); osh_scoring_workspace_free(ws); return OSH_ENOMEM; }

    section = SECTION_NONE;
    lineno  = 0u;

    while (fgets(line, (int) sizeof(line), fp) != NULL) {
        char raw[4096];
        char *p;
        char *words[32];
        int nwords;
        int indented;

        ++lineno;
        /* Keep an unmodified copy for error messages only. */
        strncpy(raw, line, sizeof(raw) - 1u);
        raw[sizeof(raw) - 1u] = '\0';

        strip_comment(line);
        indented = is_indented(line);
        p        = trim(line);
        if (!p || p[0] == '\0') continue;

        nwords = osh_tokenise(p, words, (int)(sizeof(words) / sizeof(words[0])));
        if (nwords <= 0) continue;

        /* ---- Section header (unindented keyword) ------------------------- */
        if (!indented) {
            if (strcasecmp(words[0], "Filter") == 0) {
                rc = append_filter(ws);
                if (rc != OSH_OK) goto fail;
                section = SECTION_FILTER;
                continue;
            }
            if (strcasecmp(words[0], "Settings") == 0) {
                rc = append_settings(ws);
                if (rc != OSH_OK) goto fail;
                section = SECTION_SETTINGS;
                continue;
            }
            if (strcasecmp(words[0], "Geometry") == 0) {
                rc = append_geometry(ws);
                if (rc != OSH_OK) goto fail;
                section = SECTION_GEOMETRY;
                /* Geometry type is the word immediately after the keyword. */
                if (nwords >= 2) {
                    ws->geometries[ws->ngeometries - 1u].kind = dupstr(words[1]);
                    if (!ws->geometries[ws->ngeometries - 1u].kind) {
                        rc = OSH_ENOMEM;
                        goto fail;
                    }
                }
                continue;
            }
            if (strcasecmp(words[0], "Output") == 0) {
                rc = append_output(ws);
                if (rc != OSH_OK) goto fail;
                section = SECTION_OUTPUT;
                continue;
            }

            osh_error("%s:%u: unexpected unindented line outside a section: '%s'",
                      path, lineno, trim(raw));
            rc = OSH_EPARSE;
            goto fail;
        }

        /* ---- Section body (indented line) -------------------------------- */
        switch (section) {
        case SECTION_FILTER:
            rc = handle_filter_line(ws, words, nwords, path, lineno);
            break;
        case SECTION_SETTINGS:
            rc = handle_settings_line(ws, words, nwords, path, lineno);
            break;
        case SECTION_GEOMETRY:
            rc = handle_geometry_line(ws, words, nwords, path, lineno);
            break;
        case SECTION_OUTPUT:
            rc = handle_output_line(ws, words, nwords, path, lineno);
            break;
        default:
            osh_error("%s:%u: indented line with no active section", path, lineno);
            rc = OSH_EPARSE;
            break;
        }
        if (rc != OSH_OK) goto fail;
    }

    fclose(fp);

    rc = validate(ws);
    if (rc != OSH_OK) {
        osh_scoring_workspace_free(ws);
        return rc;
    }

    *ws_out = ws;
    return OSH_OK;

fail:
    fclose(fp);
    osh_scoring_workspace_free(ws);
    return rc;
}
