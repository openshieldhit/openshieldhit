#include "openshieldhit/openshieldhit.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "beam/osh_beam.h"
#include "common/osh_logger.h"
#include "common/osh_rc.h"
#include "common/osh_version.h"
#include "gemca/osh_gemca2.h"
#include "material/osh_material.h"
#include "material/runtime/osh_material_prepare.h"
#include "scoring/osh_scoring.h"
#include "scoring/runtime/osh_scoring_prepare.h"

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

/* Resolve either an explicit override path or workdir/filename into an owned
 * path string for the context. */
static char *resolve_input_path(char const *workdir, char const *override_path, char const *filename);
/* Duplicate a C string into owned heap storage. */
static char *duplicate_cstring(char const *src);
/* Open a file in read-only mode using a Windows-safe wrapper around fopen. */
static FILE *fopen_readonly(char const *path);
/* Return 1 if a readable file exists at path, 0 otherwise. */
static int file_exists(char const *path);
/* Duplicate a config string so the context owns stable storage. */
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
    struct beam_workspace *beam = NULL;
    struct gemca_workspace *geom = NULL;
    struct material_workspace *mat = NULL;
    struct osh_scoring_workspace *scoring = NULL;
    struct osh_material_runtime transport_tables;
    struct osh_scoring_runtime scoring_runtime;
    enum openshieldhit_status rc = OPENSHIELDHIT_STATUS_OK;

    if (!ctx) {
        return OPENSHIELDHIT_STATUS_INVALID_ARGUMENT;
    }

    ctx->last_error[0] = '\0';
    memset(&transport_tables, 0, sizeof(transport_tables));
    memset(&scoring_runtime, 0, sizeof(scoring_runtime));

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

    if (osh_gemca_workspace_init(&geom) != OSH_OK) {
        ctx_set_error(ctx, err, "%s", "could not allocate geometry workspace");
        rc = OPENSHIELDHIT_STATUS_NO_MEMORY;
        goto cleanup;
    }

    if (osh_gemca_load(geo_path, geom) != OSH_OK) {
        ctx_set_error(ctx, err, "failed to load geometry: %s", geo_path);
        rc = OPENSHIELDHIT_STATUS_PARSE_ERROR;
        goto cleanup;
    }
    if (out) {
        fprintf(out, "Loaded geometry: %s\n", geo_path);
    }

    if (!file_exists(beam_path)) {
        ctx_set_error(ctx, err, "beam file not found: %s", beam_path);
        rc = OPENSHIELDHIT_STATUS_IO_ERROR;
        goto cleanup;
    }

    if (osh_beam_setup_from_path(beam_path, NULL, &beam) != OSH_OK) {
        ctx_set_error(ctx, err, "failed to load beam: %s", beam_path);
        rc = OPENSHIELDHIT_STATUS_PARSE_ERROR;
        goto cleanup;
    }
    if (ctx->has_nstat) {
        beam->nstat = (size_t) ctx->nstat;
    }
    if (out) {
        fprintf(out, "Loaded beam: %s\n", beam_path);
        if (ctx->has_nstat) {
            fprintf(out, "Applied nstat override: %llu\n", ctx->nstat);
        }
    }

    if (!file_exists(mat_path)) {
        ctx_set_error(ctx, err, "material file not found: %s", mat_path);
        rc = OPENSHIELDHIT_STATUS_IO_ERROR;
        goto cleanup;
    }

    if (osh_material_setup_from_path(mat_path, NULL, &mat) != OSH_OK) {
        ctx_set_error(ctx, err, "failed to load materials: %s", mat_path);
        rc = OPENSHIELDHIT_STATUS_PARSE_ERROR;
        goto cleanup;
    }
    if (out) {
        fprintf(out, "Loaded materials: %s\n", mat_path);
    }

    /* Resolve zone->material_name strings to dense material indices. */
    {
        size_t iz;
        for (iz = 0; iz < geom->nzones; iz++) {
            struct zone *z;
            struct material const *m;

            z = geom->zones[iz];
            if (!z->material_name) {
                ctx_set_error(ctx, err, "zone '%s' has no material assigned", z->name ? z->name : "(unnamed)");
                rc = OPENSHIELDHIT_STATUS_PARSE_ERROR;
                goto cleanup;
            }
            m = osh_material_by_name(mat, z->material_name);
            if (!m) {
                char msg[OSH_LAST_ERROR_SIZE];
                snprintf(msg,
                         sizeof(msg),
                         "zone '%s': unknown material '%s'",
                         z->name ? z->name : "(unnamed)",
                         z->material_name);
                ctx_set_error(ctx, err, "%s", msg);
                rc = OPENSHIELDHIT_STATUS_PARSE_ERROR;
                goto cleanup;
            }
            z->material_idx = m->index;
        }
    }
    if (out) {
        fprintf(out, "Material assembly complete: %llu zones resolved.\n", (unsigned long long) geom->nzones);
    }

    {
        unsigned int z_max;

        z_max = (beam->primary.z > 0u) ? (unsigned int) beam->primary.z : 1u;
        if (osh_material_prepare(mat, z_max, &transport_tables) != OSH_OK) {
            ctx_set_error(ctx, err, "%s", "failed to prepare runtime transport tables");
            rc = OPENSHIELDHIT_STATUS_PARSE_ERROR;
            goto cleanup;
        }
    }
    if (out) {
        fprintf(out,
                "Transport tables prepared: %llu materials x %llu projectiles x %llu energies.\n",
                (unsigned long long) transport_tables.nmaterials,
                (unsigned long long) transport_tables.nprojectiles,
                (unsigned long long) transport_tables.nenergy);
    }

    if (!file_exists(detect_path)) {
        ctx_set_error(ctx, err, "detect file not found: %s", detect_path);
        rc = OPENSHIELDHIT_STATUS_IO_ERROR;
        goto cleanup;
    }

    if (osh_scoring_setup_from_path(detect_path, &scoring) != OSH_OK) {
        ctx_set_error(ctx, err, "failed to load scoring/detect input: %s", detect_path);
        rc = OPENSHIELDHIT_STATUS_PARSE_ERROR;
        goto cleanup;
    }
    if (out) {
        fprintf(out, "Loaded scoring: %s\n", detect_path);
    }

    if (osh_scoring_prepare(scoring, &scoring_runtime) != OSH_OK) {
        ctx_set_error(ctx, err, "failed to prepare scoring runtime: %s", detect_path);
        rc = OPENSHIELDHIT_STATUS_PARSE_ERROR;
        goto cleanup;
    }
    if (out) {
        fprintf(out,
                "Scoring runtime prepared: %llu geometries, %llu outputs, %llu pages.\n",
                (unsigned long long) scoring_runtime.ngeometries,
                (unsigned long long) scoring_runtime.noutputs,
                (unsigned long long) scoring_runtime.npages);
    }

    if (out) {
        fprintf(out, "Validation completed.\n");
    }

cleanup:
    osh_scoring_runtime_free(&scoring_runtime);
    if (scoring) {
        osh_scoring_workspace_free(scoring);
    }
    osh_material_runtime_free(&transport_tables);
    if (mat) {
        osh_material_workspace_free(mat);
    }
    if (beam) {
        osh_beam_workspace_free(beam);
    }
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

/**
 * @brief Resolve an input path into owned storage for the context.
 *
 * If @p override_path is given, it is duplicated directly. Otherwise the
 * function joins @p workdir and @p filename into a newly allocated path.
 *
 * @param[in] workdir       Working directory used for default file lookup.
 * @param[in] override_path Optional explicit path from the caller.
 * @param[in] filename      Default file name to append to @p workdir.
 *
 * @return Newly allocated path string, or NULL on allocation/argument error.
 */
static char *resolve_input_path(char const *workdir, char const *override_path, char const *filename) {
    char *path;
    size_t wlen;
    size_t flen;

    if (override_path && override_path[0]) {
        path = duplicate_cstring(override_path);
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

    memcpy(path, workdir, wlen);
    path[wlen] = '\0';
    if ((wlen > 0) && (workdir[wlen - 1] != '/')) {
        path[wlen] = '/';
        path[wlen + 1] = '\0';
        wlen += 1u;
    }
    memcpy(path + wlen, filename, flen + 1u);
    return path;
}

/**
 * @brief Duplicate a C string into owned heap storage.
 *
 * This helper exists to avoid platform-specific strdup variants and to make
 * ownership explicit in code paths where the context needs stable copies of
 * caller-provided strings.
 *
 * @param[in] src Source string to duplicate.
 *
 * @return Newly allocated copy, or NULL on error.
 */
static char *duplicate_cstring(char const *src) {
    char *dst;
    size_t len;

    if (!src) {
        return NULL;
    }

    len = strlen(src) + 1u;
    dst = (char *) malloc(len);
    if (!dst) {
        return NULL;
    }
    memcpy(dst, src, len);
    return dst;
}

/**
 * @brief Open a file for read-only probing in a portable way.
 *
 * Uses fopen_s on MSVC to avoid C4996 warnings, and falls back to fopen on
 * other platforms.
 *
 * @param[in] path File path to open.
 *
 * @return Open FILE* on success, or NULL on failure.
 */
static FILE *fopen_readonly(char const *path) {
#if defined(_MSC_VER)
    FILE *fp = NULL;
    if (fopen_s(&fp, path, "r") != 0) {
        return NULL;
    }
    return fp;
#else
    return fopen(path, "r");
#endif
}

/**
 * @brief Check whether a readable file exists at @p path.
 *
 * @param[in] path File path to probe.
 *
 * @return 1 if the file can be opened for reading, otherwise 0.
 */
static int file_exists(char const *path) {
    FILE *fp;
    if (!path || !path[0]) {
        return 0;
    }
    fp = fopen_readonly(path);
    if (!fp) {
        return 0;
    }
    fclose(fp);
    return 1;
}

/**
 * @brief Duplicate an optional config string into context-owned storage.
 *
 * The configuration struct may point to stack memory or other short-lived
 * storage. The context therefore copies any non-empty string it accepts and
 * frees that memory in openshieldhit_context_destroy().
 *
 * @param[in]  src     Optional source string; NULL/empty means "unset".
 * @param[out] dst_out Receives the allocated copy, or NULL when unset.
 *
 * @return 1 on success, 0 on allocation/argument error.
 */
static int duplicate_string(char const *src, char **dst_out) {
    char *copy = NULL;

    if (!dst_out) {
        return 0;
    }
    if (src && src[0]) {
        copy = duplicate_cstring(src);
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
