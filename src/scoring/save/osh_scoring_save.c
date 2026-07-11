#include "scoring/save/osh_scoring_save.h"

#include <string.h>

#include "common/osh_diag.h"
#include "scoring/save/osh_scoring_save_ascii.h"
#include "scoring/save/osh_scoring_save_bdo2019.h"
#include "scoring/save/osh_scoring_save_plot.h"
#include "scoring/save/osh_scoring_save_rtdose.h"

static enum osh_status save_dispatch(struct osh_scoring_workspace const *ws,
                                     struct osh_scoring_runtime const *rt,
                                     unsigned long long nstat,
                                     size_t const *want,
                                     size_t n_want,
                                     struct osh_diag_sink const *diag);
static enum osh_status save_one_output(struct osh_scoring_workspace const *ws,
                                       struct osh_scoring_runtime const *rt,
                                       unsigned long long nstat,
                                       size_t output_idx,
                                       struct osh_diag_sink const *diag);
static int fileformat_is_ascii(char const *fileformat);
static int fileformat_is_bdo2019(char const *fileformat);
static int fileformat_is_rtdose(char const *fileformat);
static int fileformat_is_plot(char const *fileformat);
static char const *fileformat_label(char const *fileformat);

enum osh_status osh_scoring_save(struct osh_scoring_workspace const *ws,
                                 struct osh_scoring_runtime const *rt,
                                 unsigned long long nstat) {
    return save_dispatch(ws, rt, nstat, NULL, 0u, NULL);
}

enum osh_status osh_scoring_save_diag(struct osh_scoring_workspace const *ws,
                                      struct osh_scoring_runtime const *rt,
                                      unsigned long long nstat,
                                      struct osh_diag_sink const *diag) {
    return save_dispatch(ws, rt, nstat, NULL, 0u, diag);
}

enum osh_status osh_scoring_save_outputs(struct osh_scoring_workspace const *ws,
                                         struct osh_scoring_runtime const *rt,
                                         unsigned long long nstat,
                                         size_t const *want,
                                         size_t n_want) {
    return save_dispatch(ws, rt, nstat, want, n_want, NULL);
}

/* Shared worker for every public entry point.  Argument/state validation is
 * fail-fast; the actual writing is best-effort (issue #308): every requested
 * target is attempted and the first error observed is returned, so one bad path
 * neither aborts the others nor hides itself from the exit status. */
static enum osh_status save_dispatch(struct osh_scoring_workspace const *ws,
                                     struct osh_scoring_runtime const *rt,
                                     unsigned long long nstat,
                                     size_t const *want,
                                     size_t n_want,
                                     struct osh_diag_sink const *diag) {
    enum osh_status rc;
    enum osh_status first_rc;
    size_t i;

    if (!ws || !rt) {
        return OSH_EINVAL;
    }
    if (nstat == 0ull) {
        return OSH_EINVAL;
    }
    /* Multi-format blocks fan out into extra runtime outputs (issue #308), so the
     * runtime holds at least as many outputs as the cold workspace — never fewer.
     * Fewer means a workspace/runtime mismatch. */
    if (rt->noutputs < ws->noutputs) {
        return OSH_ESTATE;
    }

    first_rc = OSH_OK;
    if (want == NULL) {
        for (i = 0; i < rt->noutputs; ++i) {
            rc = save_one_output(ws, rt, nstat, i, diag);
            if (rc != OSH_OK && first_rc == OSH_OK) {
                first_rc = rc;
            }
        }
        return first_rc;
    }

    for (i = 0; i < n_want; ++i) {
        if (want[i] >= rt->noutputs) {
            return OSH_EINVAL;
        }
        rc = save_one_output(ws, rt, nstat, want[i], diag);
        if (rc != OSH_OK && first_rc == OSH_OK) {
            first_rc = rc;
        }
    }
    return first_rc;
}

static enum osh_status save_one_output(struct osh_scoring_workspace const *ws,
                                       struct osh_scoring_runtime const *rt,
                                       unsigned long long nstat,
                                       size_t output_idx,
                                       struct osh_diag_sink const *diag) {
    char const *fileformat;
    char const *filename;
    enum osh_status rc;

    if (output_idx >= rt->noutputs) {
        return OSH_EINVAL;
    }

    /* Dispatch on the runtime output's own format, not the cold workspace's:
     * a fanned-out multi-format target (issue #308) has no cold counterpart. */
    fileformat = rt->outputs[output_idx].fileformat;
    if (fileformat_is_ascii(fileformat)) {
        rc = osh_scoring_save_ascii_output(ws, rt, nstat, output_idx);
    } else if (fileformat_is_bdo2019(fileformat)) {
        rc = osh_scoring_save_bdo2019_output(ws, rt, nstat, output_idx);
    } else if (fileformat_is_rtdose(fileformat)) {
        rc = osh_scoring_save_rtdose_output(ws, rt, nstat, output_idx);
    } else if (fileformat_is_plot(fileformat)) {
        rc = osh_scoring_save_plot_output(ws, rt, nstat, output_idx);
    } else {
        rc = OSH_ENOTSUP;
    }

    if (rc != OSH_OK) {
        filename = rt->outputs[output_idx].filename;
        if (!filename) {
            filename = "(unnamed)";
        }
        OSH_DIAG_WARNF(diag,
                       "scoring: could not write output '%s' as %s (status %d); continuing with remaining targets",
                       filename,
                       fileformat_label(fileformat),
                       (int) rc);
    }
    return rc;
}

static int fileformat_is_ascii(char const *fileformat) {
    if (!fileformat) {
        return 0;
    }
    return (strcmp(fileformat, "text") == 0) || (strcmp(fileformat, "txt") == 0) || (strcmp(fileformat, "ascii") == 0)
           || (strcmp(fileformat, "dat") == 0);
}

static int fileformat_is_bdo2019(char const *fileformat) {
    if (!fileformat) {
        return 1;
    }
    return (strcmp(fileformat, "bdo") == 0) || (strcmp(fileformat, "bdo2019") == 0)
           || (strcmp(fileformat, "binary") == 0) || (strcmp(fileformat, "bin") == 0);
}

static int fileformat_is_rtdose(char const *fileformat) {
    if (!fileformat) {
        return 0;
    }
    /* Format keywords are lowercased at parse time, so only "rtdose" can match —
     * keep this consistent with the ascii/bdo/plot predicates above. */
    return strcmp(fileformat, "rtdose") == 0;
}

/* Native "quick-look" plot formats (issue #238).  Only SVG is implemented; the
 * format keyword is stored lowercased by the parser. */
static int fileformat_is_plot(char const *fileformat) {
    if (!fileformat) {
        return 0;
    }
    return strcmp(fileformat, "svg") == 0;
}

/* Display label for diagnostics; a NULL runtime format means the default BDO. */
static char const *fileformat_label(char const *fileformat) {
    if (!fileformat) {
        return "bdo (default)";
    }
    return fileformat;
}
