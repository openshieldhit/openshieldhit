#include "openshieldhit/openshieldhit.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cli/osh_cli.h"
#include "common/osh_rc.h"
#include "common/osh_version.h"
#include "gemca/osh_gemca2.h"

static char const *const OPENSHIELDHIT_DEFAULT_WORKDIR = ".";
static char const *const OPENSHIELDHIT_GEO_FILENAME = "geo.dat";
static char const *const OPENSHIELDHIT_BEAM_FILENAME = "beam.dat";
static char const *const OPENSHIELDHIT_MAT_FILENAME = "mat.dat";
static char const *const OPENSHIELDHIT_DETECT_FILENAME = "detect.dat";

static void copy_cli_options(struct openshieldhit_cli_options *dst, struct osh_cli_options const *src);
static int map_osh_status(int rc);
static char *resolve_input_path(char const *workdir, char const *override_path, char const *filename);
static int file_exists(char const *path);

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

int openshieldhit_cli_parse(int argc, char *argv[], struct openshieldhit_cli_options *opt, char *err, size_t err_cap) {
    int rc;
    struct osh_cli_options internal_opt;

    if (!opt) {
        return OPENSHIELDHIT_STATUS_INVALID_ARGUMENT;
    }

    rc = osh_cli_parse(argc, argv, &internal_opt, err, err_cap);
    if (rc != OSH_OK) {
        return map_osh_status(rc);
    }

    copy_cli_options(opt, &internal_opt);
    return OPENSHIELDHIT_STATUS_OK;
}

void openshieldhit_cli_print_help(FILE *out, char const *prog) {
    osh_cli_print_help(out, prog);
}

int openshieldhit_run(struct openshieldhit_cli_options const *opt, FILE *out, FILE *err) {
    char const *workdir;
    char const *outdir;
    char *geo_path = NULL;
    char *beam_path = NULL;
    char *mat_path = NULL;
    char *detect_path = NULL;
    struct gemca_workspace *geom = NULL;

    if (!opt) {
        return OPENSHIELDHIT_STATUS_INVALID_ARGUMENT;
    }

    if (!out) {
        out = stdout;
    }
    if (!err) {
        err = stderr;
    }

    workdir = (opt->workdir && opt->workdir[0]) ? opt->workdir : OPENSHIELDHIT_DEFAULT_WORKDIR;
    outdir = (opt->out_dir && opt->out_dir[0]) ? opt->out_dir : workdir;

    geo_path = resolve_input_path(workdir, opt->geo_path, OPENSHIELDHIT_GEO_FILENAME);
    beam_path = resolve_input_path(workdir, opt->beam_path, OPENSHIELDHIT_BEAM_FILENAME);
    mat_path = resolve_input_path(workdir, opt->mat_path, OPENSHIELDHIT_MAT_FILENAME);
    detect_path = resolve_input_path(workdir, opt->detect_path, OPENSHIELDHIT_DETECT_FILENAME);
    if (!geo_path || !beam_path || !mat_path || !detect_path) {
        fprintf(err, "Error: out of memory while resolving input paths\n");
        free(geo_path);
        free(beam_path);
        free(mat_path);
        free(detect_path);
        return OPENSHIELDHIT_STATUS_NO_MEMORY;
    }

    if (!opt->dry_run) {
        fprintf(out, "OpenShieldHIT version %s\n", openshieldhit_version_string());
        fprintf(out, "No run mode selected. Use --dry-run to validate input loading.\n");
        fprintf(out, "Verbosity level  : %d\n", opt->verbose);
        if (opt->has_nstat) {
            fprintf(out, "Requested nstat  : %llu\n", opt->nstat);
        }
        fprintf(out, "Working directory: %s\n", workdir);
        fprintf(out, "Output directory : %s\n", outdir);
        fprintf(out, "Geometry input   : %s\n", geo_path);
        fprintf(out, "Beam input       : %s\n", beam_path);
        fprintf(out, "Material input   : %s\n", mat_path);
        fprintf(out, "Detect input     : %s\n", detect_path);
        free(geo_path);
        free(beam_path);
        free(mat_path);
        free(detect_path);
        return OPENSHIELDHIT_STATUS_OK;
    }

    fprintf(out, "Dry-run configuration\n");
    fprintf(out, "  Verbosity level  : %d\n", opt->verbose);
    if (opt->has_nstat) {
        fprintf(out, "  Requested nstat  : %llu\n", opt->nstat);
    }
    fprintf(out, "  Working directory: %s\n", workdir);
    fprintf(out, "  Output directory : %s\n", outdir);
    fprintf(out, "  Geometry input   : %s\n", geo_path);
    fprintf(out, "  Beam input       : %s\n", beam_path);
    fprintf(out, "  Material input   : %s\n", mat_path);
    fprintf(out, "  Detect input     : %s\n", detect_path);

    if (!file_exists(geo_path)) {
        fprintf(err, "Error: geometry file not found: %s\n", geo_path);
        free(geo_path);
        free(beam_path);
        free(mat_path);
        free(detect_path);
        return OPENSHIELDHIT_STATUS_IO_ERROR;
    }

    if (!osh_gemca_workspace_init(&geom)) {
        fprintf(err, "Error: could not allocate geometry workspace\n");
        free(geo_path);
        free(beam_path);
        free(mat_path);
        free(detect_path);
        return OPENSHIELDHIT_STATUS_NO_MEMORY;
    }

    if (!osh_gemca_load(geo_path, geom)) {
        fprintf(err, "Error: failed to load geometry '%s'\n", geo_path);
        osh_gemca_workspace_free(geom);
        free(geo_path);
        free(beam_path);
        free(mat_path);
        free(detect_path);
        return OPENSHIELDHIT_STATUS_PARSE_ERROR;
    }
    fprintf(out, "Loaded geometry: %s\n", geo_path);

    if (file_exists(beam_path)) {
        fprintf(out, "Beam file found (loader wiring pending): %s\n", beam_path);
    } else {
        fprintf(out, "Beam file missing (loader wiring pending): %s\n", beam_path);
    }
    if (file_exists(mat_path)) {
        fprintf(out, "Material file found (loader wiring pending): %s\n", mat_path);
    } else {
        fprintf(out, "Material file missing (loader wiring pending): %s\n", mat_path);
    }
    if (file_exists(detect_path)) {
        fprintf(out, "Detect file found (loader wiring pending): %s\n", detect_path);
    } else {
        fprintf(out, "Detect file missing (loader wiring pending): %s\n", detect_path);
    }

    if (geom) {
        osh_gemca_workspace_free(geom);
    }

    free(geo_path);
    free(beam_path);
    free(mat_path);
    free(detect_path);

    fprintf(out, "Dry-run completed.\n");
    return OPENSHIELDHIT_STATUS_OK;
}

static void copy_cli_options(struct openshieldhit_cli_options *dst, struct osh_cli_options const *src) {
    dst->action = (enum openshieldhit_cli_action) src->action;
    dst->dry_run = src->dry_run;
    dst->verbose = src->verbose;
    dst->workdir = src->workdir;
    dst->geo_path = src->geo_path;
    dst->beam_path = src->beam_path;
    dst->mat_path = src->mat_path;
    dst->detect_path = src->detect_path;
    dst->out_dir = src->out_dir;
    dst->nstat = src->nstat;
    dst->has_nstat = src->has_nstat;
}

static int map_osh_status(int rc) {
    switch (rc) {
    case OSH_OK:
        return OPENSHIELDHIT_STATUS_OK;
    case OSH_EINVAL:
        return OPENSHIELDHIT_STATUS_INVALID_ARGUMENT;
    case OSH_ENOMEM:
        return OPENSHIELDHIT_STATUS_NO_MEMORY;
    case OSH_EIO:
        return OPENSHIELDHIT_STATUS_IO_ERROR;
    case OSH_EPARSE:
        return OPENSHIELDHIT_STATUS_PARSE_ERROR;
    case OSH_EINCOMPLETE:
        return OPENSHIELDHIT_STATUS_INCOMPLETE;
    case OSH_ENOTSUP:
        return OPENSHIELDHIT_STATUS_NOT_SUPPORTED;
    case OSH_ESTATE:
        return OPENSHIELDHIT_STATUS_STATE_ERROR;
    default:
        return OPENSHIELDHIT_STATUS_STATE_ERROR;
    }
}

static char *resolve_input_path(char const *workdir, char const *override_path, char const *filename) {
    char *path = NULL;
    size_t wlen;
    size_t flen;

    if (override_path && override_path[0]) {
        path = (char *) malloc(strlen(override_path) + 1);
        if (!path) {
            return NULL;
        }
        strcpy(path, override_path);
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
