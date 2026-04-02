#include "openshieldhit/openshieldhit.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_version.h"
#include "gemca/osh_gemca2.h"

/* ---- Internal context definition ---------------------------------------- */

struct openshieldhit_context {
    /* Input paths — owned by the context; freed in destroy(). */
    char *workdir;
    char *out_dir;
    char *geo_path;
    char *beam_path;
    char *mat_path;
    char *detect_path;

    /* Run options */
    unsigned long long nstat;
    int has_nstat;
    enum openshieldhit_run_mode run_mode;
    int log_level;

    /* Last error message, populated by openshieldhit_run() on failure.
     * Future: could be extended to a linked list of diagnostics. */
    char last_error[256];
};

/* ---- Default file names -------------------------------------------------- */

static char const *const OSH_DEFAULT_WORKDIR = ".";
static char const *const OSH_GEO_FILENAME = "geo.dat";
static char const *const OSH_BEAM_FILENAME = "beam.dat";
static char const *const OSH_MAT_FILENAME = "mat.dat";
static char const *const OSH_DETECT_FILENAME = "detect.dat";

/* ---- Internal helpers ---------------------------------------------------- */

static char *resolve_input_path(char const *workdir, char const *override_path, char const *filename);
static int file_exists(char const *path);

/* Sets ctx->last_error and optionally prints to err (NULL = silent). */
static void ctx_set_error(openshieldhit_context_t *ctx, FILE *err, char const *fmt, char const *detail) {
    if (ctx) {
        snprintf(ctx->last_error, sizeof(ctx->last_error), fmt, detail);
    }
    if (err) {
        fprintf(err, "Error: ");
        fprintf(err, fmt, detail);
        fprintf(err, "\n");
    }
}

/* Deep-copies src into *dst. On OOM the old value is preserved and 0 is
 * returned. Passing NULL/empty clears the field (always succeeds). */
static int ctx_set_string(char **dst, char const *src) {
    char *copy;
    if (!src || !src[0]) {
        free(*dst);
        *dst = NULL;
        return 1;
    }
    copy = strdup(src);
    if (!copy) {
        return 0; /* old value untouched */
    }
    free(*dst);
    *dst = copy;
    return 1;
}

/* ---- Lifecycle ----------------------------------------------------------- */

openshieldhit_context_t *openshieldhit_context_create(void) {
    return (openshieldhit_context_t *) calloc(1, sizeof(openshieldhit_context_t));
}

void openshieldhit_context_destroy(openshieldhit_context_t *ctx) {
    if (!ctx) {
        return;
    }
    free(ctx->workdir);
    free(ctx->out_dir);
    free(ctx->geo_path);
    free(ctx->beam_path);
    free(ctx->mat_path);
    free(ctx->detect_path);
    free(ctx);
}

/* ---- Configuration setters ----------------------------------------------- */

enum openshieldhit_status openshieldhit_context_set_workdir(openshieldhit_context_t *ctx, char const *path) {
    if (!ctx) {
        return OPENSHIELDHIT_STATUS_INVALID_ARGUMENT;
    }
    return ctx_set_string(&ctx->workdir, path) ? OPENSHIELDHIT_STATUS_OK : OPENSHIELDHIT_STATUS_NO_MEMORY;
}

enum openshieldhit_status openshieldhit_context_set_out_dir(openshieldhit_context_t *ctx, char const *path) {
    if (!ctx) {
        return OPENSHIELDHIT_STATUS_INVALID_ARGUMENT;
    }
    return ctx_set_string(&ctx->out_dir, path) ? OPENSHIELDHIT_STATUS_OK : OPENSHIELDHIT_STATUS_NO_MEMORY;
}

enum openshieldhit_status openshieldhit_context_set_run_mode(openshieldhit_context_t *ctx,
                                                             enum openshieldhit_run_mode mode) {
    if (!ctx) {
        return OPENSHIELDHIT_STATUS_INVALID_ARGUMENT;
    }
    ctx->run_mode = mode;
    return OPENSHIELDHIT_STATUS_OK;
}

enum openshieldhit_status openshieldhit_context_set_log_level(openshieldhit_context_t *ctx, int level) {
    if (!ctx) {
        return OPENSHIELDHIT_STATUS_INVALID_ARGUMENT;
    }
    ctx->log_level = level;
    return OPENSHIELDHIT_STATUS_OK;
}

enum openshieldhit_status openshieldhit_context_set_nstat(openshieldhit_context_t *ctx, unsigned long long nstat) {
    if (!ctx) {
        return OPENSHIELDHIT_STATUS_INVALID_ARGUMENT;
    }
    ctx->nstat = nstat;
    ctx->has_nstat = 1;
    return OPENSHIELDHIT_STATUS_OK;
}

enum openshieldhit_status openshieldhit_context_set_geo_path(openshieldhit_context_t *ctx, char const *path) {
    if (!ctx) {
        return OPENSHIELDHIT_STATUS_INVALID_ARGUMENT;
    }
    return ctx_set_string(&ctx->geo_path, path) ? OPENSHIELDHIT_STATUS_OK : OPENSHIELDHIT_STATUS_NO_MEMORY;
}

enum openshieldhit_status openshieldhit_context_set_beam_path(openshieldhit_context_t *ctx, char const *path) {
    if (!ctx) {
        return OPENSHIELDHIT_STATUS_INVALID_ARGUMENT;
    }
    return ctx_set_string(&ctx->beam_path, path) ? OPENSHIELDHIT_STATUS_OK : OPENSHIELDHIT_STATUS_NO_MEMORY;
}

enum openshieldhit_status openshieldhit_context_set_mat_path(openshieldhit_context_t *ctx, char const *path) {
    if (!ctx) {
        return OPENSHIELDHIT_STATUS_INVALID_ARGUMENT;
    }
    return ctx_set_string(&ctx->mat_path, path) ? OPENSHIELDHIT_STATUS_OK : OPENSHIELDHIT_STATUS_NO_MEMORY;
}

enum openshieldhit_status openshieldhit_context_set_detect_path(openshieldhit_context_t *ctx, char const *path) {
    if (!ctx) {
        return OPENSHIELDHIT_STATUS_INVALID_ARGUMENT;
    }
    return ctx_set_string(&ctx->detect_path, path) ? OPENSHIELDHIT_STATUS_OK : OPENSHIELDHIT_STATUS_NO_MEMORY;
}

/* ---- Error retrieval ----------------------------------------------------- */

char const *openshieldhit_last_error(openshieldhit_context_t const *ctx) {
    if (!ctx) {
        return "";
    }
    return ctx->last_error;
}

/* ---- Version ------------------------------------------------------------- */

char const *openshieldhit_version_string(void) {
    return OSH_VERSION;
}

int openshieldhit_version_major(void) {
    return OSH_VERSION_MAJOR;
}

int openshieldhit_version_minor(void) {
    return OSH_VERSION_MINOR;
}

int openshieldhit_version_patch(void) {
    return OSH_VERSION_PATCH;
}

/* ---- Run ----------------------------------------------------------------- */

enum openshieldhit_status openshieldhit_run(openshieldhit_context_t *ctx, FILE *out, FILE *err) {
    char const *workdir;
    char const *outdir;
    char *geo_path = NULL;
    char *beam_path = NULL;
    char *mat_path = NULL;
    char *detect_path = NULL;
    struct gemca_workspace *geom = NULL;
    enum openshieldhit_status rc = OPENSHIELDHIT_STATUS_OK;

    if (!ctx) {
        return OPENSHIELDHIT_STATUS_INVALID_ARGUMENT;
    }

    ctx->last_error[0] = '\0';

    workdir = (ctx->workdir && ctx->workdir[0]) ? ctx->workdir : OSH_DEFAULT_WORKDIR;
    outdir = (ctx->out_dir && ctx->out_dir[0]) ? ctx->out_dir : workdir;

    geo_path = resolve_input_path(workdir, ctx->geo_path, OSH_GEO_FILENAME);
    beam_path = resolve_input_path(workdir, ctx->beam_path, OSH_BEAM_FILENAME);
    mat_path = resolve_input_path(workdir, ctx->mat_path, OSH_MAT_FILENAME);
    detect_path = resolve_input_path(workdir, ctx->detect_path, OSH_DETECT_FILENAME);

    if (!geo_path || !beam_path || !mat_path || !detect_path) {
        ctx_set_error(ctx, err, "%s", "out of memory while resolving input paths");
        rc = OPENSHIELDHIT_STATUS_NO_MEMORY;
        goto cleanup;
    }

    if (ctx->run_mode != OPENSHIELDHIT_RUN_VALIDATE) {
        /* Full transport is not yet implemented. Return NOT_SUPPORTED so
         * callers can distinguish this from a successful run. */
        ctx_set_error(ctx, err, "%s", "run mode NORMAL is not yet implemented");
        rc = OPENSHIELDHIT_STATUS_NOT_SUPPORTED;
        goto cleanup;
    }

    if (out) {
        fprintf(out, "Validate configuration\n");
        fprintf(out, "  Log level        : %d\n", ctx->log_level);
        if (ctx->has_nstat) {
            fprintf(out, "  Requested nstat  : %llu\n", ctx->nstat);
        }
        fprintf(out, "  Working directory: %s\n", workdir);
        fprintf(out, "  Output directory : %s\n", outdir);
        fprintf(out, "  Geometry input   : %s\n", geo_path);
        fprintf(out, "  Beam input       : %s\n", beam_path);
        fprintf(out, "  Material input   : %s\n", mat_path);
        fprintf(out, "  Detect input     : %s\n", detect_path);
    }

    if (!file_exists(geo_path)) {
        ctx_set_error(ctx, err, "geometry file not found: %s", geo_path);
        rc = OPENSHIELDHIT_STATUS_IO_ERROR;
        goto cleanup;
    }

    if (!osh_gemca_workspace_init(&geom)) {
        ctx_set_error(ctx, err, "%s", "could not allocate geometry workspace");
        rc = OPENSHIELDHIT_STATUS_NO_MEMORY;
        goto cleanup;
    }

    if (!osh_gemca_load(geo_path, geom)) {
        ctx_set_error(ctx, err, "failed to load geometry: %s", geo_path);
        rc = OPENSHIELDHIT_STATUS_PARSE_ERROR;
        goto cleanup;
    }
    if (out) {
        fprintf(out, "Loaded geometry: %s\n", geo_path);
    }

    /* TODO: wire beam loader */
    if (out) {
        fprintf(
            out, "Beam file %s (loader wiring pending): %s\n", file_exists(beam_path) ? "found" : "missing", beam_path);
    }

    /* TODO: wire material loader */
    if (out) {
        fprintf(out,
                "Material file %s (loader wiring pending): %s\n",
                file_exists(mat_path) ? "found" : "missing",
                mat_path);
    }

    /* TODO: wire scoring/detect loader */
    if (out) {
        fprintf(out,
                "Detect file %s (loader wiring pending): %s\n",
                file_exists(detect_path) ? "found" : "missing",
                detect_path);
    }

    if (out) {
        fprintf(out, "Validation completed.\n");
    }

cleanup:
    if (geom) {
        osh_gemca_workspace_free(geom);
    }
    free(geo_path);
    free(beam_path);
    free(mat_path);
    free(detect_path);
    return rc;
}

/* ---- Internal helpers ---------------------------------------------------- */

static char *resolve_input_path(char const *workdir, char const *override_path, char const *filename) {
    char *path;
    size_t wlen;
    size_t flen;

    if (override_path && override_path[0]) {
        path = strdup(override_path);
        return path;
    }

    if (!workdir || !filename) {
        return NULL;
    }

    wlen = strlen(workdir);
    flen = strlen(filename);
    path = (char *) malloc(wlen + 1 + flen + 1);
    if (!path) {
        return NULL;
    }

    strcpy(path, workdir);
    if ((wlen > 0) && (workdir[wlen - 1] != '/')) {
        path[wlen] = '/';
        path[wlen + 1] = '\0';
    }
    strcat(path, filename);
    return path;
}

static int file_exists(char const *path) {
    FILE *fp;
    if (!path || !path[0]) {
        return 0;
    }
    fp = fopen(path, "r");
    if (!fp) {
        return 0;
    }
    fclose(fp);
    return 1;
}
