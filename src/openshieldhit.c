#include "openshieldhit/openshieldhit.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_logger.h"
#include "common/osh_version.h"
#include "gemca/osh_gemca2.h"

/* ---- Internal constants -------------------------------------------------- */

/* Size of the last_error buffer inside the context. Large enough for a path
 * plus a short diagnostic; adjust here if longer messages are ever needed. */
#define OSH_LAST_ERROR_SIZE 256

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
    char last_error[OSH_LAST_ERROR_SIZE];
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
static int duplicate_string(char const *src, char **dst_out);
static int config_has_field(openshieldhit_config_t const *cfg, size_t field_end);

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

/* ---- Configuration ------------------------------------------------------- */

enum openshieldhit_status openshieldhit_context_configure(openshieldhit_context_t *ctx,
                                                          openshieldhit_config_t const *cfg) {
    static openshieldhit_config_t const default_cfg = OPENSHIELDHIT_CONFIG_INIT;
    char *workdir = NULL;
    char *out_dir = NULL;
    char *geo_path = NULL;
    char *beam_path = NULL;
    char *mat_path = NULL;
    char *detect_path = NULL;
    unsigned long long nstat = 0;
    int has_nstat = 0;
    enum openshieldhit_run_mode run_mode = OPENSHIELDHIT_RUN_NORMAL;
    int log_level = 0;

    if (!ctx) {
        return OPENSHIELDHIT_STATUS_INVALID_ARGUMENT;
    }

    if (!cfg) {
        cfg = &default_cfg;
    }
    if (cfg->size < sizeof(cfg->size)) {
        return OPENSHIELDHIT_STATUS_INVALID_ARGUMENT;
    }

    if (config_has_field(cfg, offsetof(openshieldhit_config_t, workdir) + sizeof(cfg->workdir))
        && !duplicate_string(cfg->workdir, &workdir)) {
        goto oom;
    }
    if (config_has_field(cfg, offsetof(openshieldhit_config_t, out_dir) + sizeof(cfg->out_dir))
        && !duplicate_string(cfg->out_dir, &out_dir)) {
        goto oom;
    }
    if (config_has_field(cfg, offsetof(openshieldhit_config_t, geo_path) + sizeof(cfg->geo_path))
        && !duplicate_string(cfg->geo_path, &geo_path)) {
        goto oom;
    }
    if (config_has_field(cfg, offsetof(openshieldhit_config_t, beam_path) + sizeof(cfg->beam_path))
        && !duplicate_string(cfg->beam_path, &beam_path)) {
        goto oom;
    }
    if (config_has_field(cfg, offsetof(openshieldhit_config_t, mat_path) + sizeof(cfg->mat_path))
        && !duplicate_string(cfg->mat_path, &mat_path)) {
        goto oom;
    }
    if (config_has_field(cfg, offsetof(openshieldhit_config_t, detect_path) + sizeof(cfg->detect_path))
        && !duplicate_string(cfg->detect_path, &detect_path)) {
        goto oom;
    }
    if (config_has_field(cfg, offsetof(openshieldhit_config_t, run_mode) + sizeof(cfg->run_mode))) {
        run_mode = cfg->run_mode;
    }
    if (config_has_field(cfg, offsetof(openshieldhit_config_t, log_level) + sizeof(cfg->log_level))) {
        log_level = cfg->log_level;
    }
    if (config_has_field(cfg, offsetof(openshieldhit_config_t, nstat) + sizeof(cfg->nstat))) {
        nstat = cfg->nstat;
    }
    if (config_has_field(cfg, offsetof(openshieldhit_config_t, has_nstat) + sizeof(cfg->has_nstat))) {
        has_nstat = cfg->has_nstat;
    }

    free(ctx->workdir);
    free(ctx->out_dir);
    free(ctx->geo_path);
    free(ctx->beam_path);
    free(ctx->mat_path);
    free(ctx->detect_path);

    ctx->workdir = workdir;
    ctx->out_dir = out_dir;
    ctx->geo_path = geo_path;
    ctx->beam_path = beam_path;
    ctx->mat_path = mat_path;
    ctx->detect_path = detect_path;
    ctx->run_mode = run_mode;
    ctx->log_level = log_level;
    ctx->nstat = nstat;
    ctx->has_nstat = has_nstat;

    return OPENSHIELDHIT_STATUS_OK;

oom:
    free(workdir);
    free(out_dir);
    free(geo_path);
    free(beam_path);
    free(mat_path);
    free(detect_path);
    return OPENSHIELDHIT_STATUS_NO_MEMORY;
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

    /* Initialise the default logger from ctx->log_level.
     *   0  → WARN  (silent — only warnings and errors)
     *   1  → INFO  (normal informational output)
     *  ≥2  → DEBUG (verbose debug output)
     * Stdout is enabled so info/debug output reaches the terminal.
     * osh_log_init() is idempotent: safe to call on every run. */
    {
        int lvl = (ctx->log_level == 0) ? OSH_LOG_WARN : (ctx->log_level == 1) ? OSH_LOG_INFO : OSH_LOG_DEBUG;
        osh_log_init(lvl, OSH_LOG_F_NONE);
        osh_log_enable_stdout(1);
    }

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

static int duplicate_string(char const *src, char **dst_out) {
    char *copy = NULL;

    if (!dst_out) {
        return 0;
    }
    if (src && src[0]) {
        copy = strdup(src);
        if (!copy) {
            return 0;
        }
    }
    *dst_out = copy;
    return 1;
}

static int config_has_field(openshieldhit_config_t const *cfg, size_t field_end) {
    return cfg && (cfg->size >= field_end);
}
