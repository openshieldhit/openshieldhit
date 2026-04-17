#include "apps/osh/osh_run.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_app_osh.h"
#include "openshieldhit/simulation.h"

/* ---- Default file names -------------------------------------------------- */

static char const *const OSH_DEFAULT_WORKDIR = ".";
static char const *const OSH_GEO_FILENAME = "geo.dat";
static char const *const OSH_BEAM_FILENAME = "beam.dat";
static char const *const OSH_MAT_FILENAME = "mat.dat";
static char const *const OSH_DETECT_FILENAME = "detect.dat";

/* ---- Internal helpers ---------------------------------------------------- */

static char *run_resolve_path(char const *workdir, char const *override_path, char const *filename);
static int run_file_exists(char const *path);

/* ---- Run ----------------------------------------------------------------- */

enum osh_status osh_run(struct osh_run_options const *opt, FILE *out, FILE *err) {
    char const *workdir;
    char const *outdir;
    char *geo_path = NULL;
    char *beam_path = NULL;
    char *mat_path = NULL;
    char *detect_path = NULL;
    struct osh_beam_workspace *beam = NULL;
    struct osh_geometry_workspace *geom = NULL;
    struct osh_material_workspace *mat = NULL;
    struct osh_scoring_workspace *scoring = NULL;
    struct osh_simulation *sim = NULL;
    enum osh_status rc = OSH_OK;

    if (!opt) {
        return OSH_EINVAL;
    }

    workdir = (opt->workdir && opt->workdir[0]) ? opt->workdir : OSH_DEFAULT_WORKDIR;
    outdir = (opt->out_dir && opt->out_dir[0]) ? opt->out_dir : workdir;

    geo_path = run_resolve_path(workdir, opt->geo_path, OSH_GEO_FILENAME);
    beam_path = run_resolve_path(workdir, opt->beam_path, OSH_BEAM_FILENAME);
    mat_path = run_resolve_path(workdir, opt->mat_path, OSH_MAT_FILENAME);
    detect_path = run_resolve_path(workdir, opt->detect_path, OSH_DETECT_FILENAME);

    if (!geo_path || !beam_path || !mat_path || !detect_path) {
        if (err) {
            fprintf(err, "Error: out of memory while resolving input paths\n");
        }
        rc = OSH_ENOMEM;
        goto cleanup;
    }

    if (out) {
        fprintf(out, "%s\n", opt->validate_only ? "Validate configuration" : "Run simulation");
        if (opt->has_nstat) {
            fprintf(out, "  Requested nstat  : %llu\n", opt->nstat);
        }
        if (opt->has_seed_offset) {
            fprintf(out, "  Requested seed offset: %llu\n", opt->seed_offset);
        }
        fprintf(out, "  Working directory: %s\n", workdir);
        fprintf(out, "  Output directory : %s\n", outdir);
        fprintf(out, "  Geometry input   : %s\n", geo_path);
        fprintf(out, "  Beam input       : %s\n", beam_path);
        fprintf(out, "  Material input   : %s\n", mat_path);
        fprintf(out, "  Detect input     : %s\n", detect_path);
    }

    if (osh_geometry_setup_from_path(geo_path, NULL, &geom) != OSH_OK) {
        if (err) {
            fprintf(err, "Error: failed to load geometry: %s\n", geo_path);
        }
        rc = OSH_EPARSE;
        goto cleanup;
    }
    if (out) {
        fprintf(out, "Loaded geometry: %s\n", geo_path);
    }

    if (!run_file_exists(beam_path)) {
        if (err) {
            fprintf(err, "Error: beam file not found: %s\n", beam_path);
        }
        rc = OSH_EIO;
        goto cleanup;
    }

    if (osh_beam_setup_from_path(beam_path, NULL, &beam) != OSH_OK) {
        if (err) {
            fprintf(err, "Error: failed to load beam: %s\n", beam_path);
        }
        rc = OSH_EPARSE;
        goto cleanup;
    }
    if (opt->has_nstat) {
        beam->nstat = (size_t) opt->nstat;
    }
    if (opt->has_seed_offset) {
        if (opt->seed_offset > 9999ull) {
            if (err) {
                fprintf(err, "Error: seed offset must be <= 9999 (got %llu)\n", opt->seed_offset);
            }
            rc = OSH_EINVAL;
            goto cleanup;
        }
        beam->rndoffset = (int) opt->seed_offset;
    }
    if (out) {
        fprintf(out, "Loaded beam: %s\n", beam_path);
        if (opt->has_nstat) {
            fprintf(out, "Applied nstat override: %llu\n", opt->nstat);
        }
        if (opt->has_seed_offset) {
            fprintf(out, "Applied seed offset override: %llu\n", opt->seed_offset);
        }
    }

    if (!run_file_exists(mat_path)) {
        if (err) {
            fprintf(err, "Error: material file not found: %s\n", mat_path);
        }
        rc = OSH_EIO;
        goto cleanup;
    }

    if (osh_material_setup_from_path(mat_path, NULL, &mat) != OSH_OK) {
        if (err) {
            fprintf(err, "Error: failed to load materials: %s\n", mat_path);
        }
        rc = OSH_EPARSE;
        goto cleanup;
    }
    if (out) {
        fprintf(out, "Loaded materials: %s\n", mat_path);
    }

    if (!run_file_exists(detect_path)) {
        if (err) {
            fprintf(err, "Error: detect file not found: %s\n", detect_path);
        }
        rc = OSH_EIO;
        goto cleanup;
    }

    if (osh_scoring_setup_from_path(detect_path, NULL, &scoring) != OSH_OK) {
        if (err) {
            fprintf(err, "Error: failed to load scoring/detect input: %s\n", detect_path);
        }
        rc = OSH_EPARSE;
        goto cleanup;
    }
    if (out) {
        fprintf(out, "Loaded scoring: %s\n", detect_path);
    }

    if (opt->validate_only) {
        if (out) {
            fprintf(out, "Validation completed.\n");
        }
        goto cleanup;
    }

    rc = osh_simulation_create(beam, geom, mat, scoring, &sim);
    if (rc != OSH_OK) {
        if (err) {
            fprintf(err, "Error: failed to compile simulation\n");
        }
        goto cleanup;
    }

    rc = osh_simulation_run(sim, outdir);
    if (rc != OSH_OK) {
        if (err) {
            fprintf(err, "Error: simulation run failed\n");
        }
        goto cleanup;
    }
    if (out) {
        fprintf(out, "Run completed. Outputs saved to %s\n", outdir);
    }

cleanup:
    osh_simulation_free(sim);
    if (scoring) {
        osh_scoring_workspace_free(scoring);
    }
    if (mat) {
        osh_material_workspace_free(mat);
    }
    if (beam) {
        osh_beam_workspace_free(beam);
    }
    if (geom) {
        osh_geometry_workspace_free(geom);
    }
    free(geo_path);
    free(beam_path);
    free(mat_path);
    free(detect_path);
    return rc;
}

/* ---- Internal helpers ---------------------------------------------------- */

static char *run_resolve_path(char const *workdir, char const *override_path, char const *filename) {
    char *path;
    size_t wlen;
    size_t flen;

    if (override_path && override_path[0]) {
        size_t len = strlen(override_path) + 1u;
        path = (char *) malloc(len);
        if (path) {
            memcpy(path, override_path, len);
        }
        return path;
    }

    if (!workdir || !filename) {
        return NULL;
    }

    wlen = strlen(workdir);
    flen = strlen(filename);
    path = (char *) malloc(wlen + 1u + flen + 1u);
    if (!path) {
        return NULL;
    }

    memcpy(path, workdir, wlen);
    path[wlen] = '\0';
    if (wlen > 0u && workdir[wlen - 1u] != '/') {
        path[wlen++] = '/';
        path[wlen] = '\0';
    }
    memcpy(path + wlen, filename, flen + 1u);
    return path;
}

static int run_file_exists(char const *path) {
    FILE *fp;
    if (!path || !path[0]) {
        return 0;
    }
#if defined(_MSC_VER)
    if (fopen_s(&fp, path, "r") != 0) {
        return 0;
    }
#else
    fp = fopen(path, "r");
#endif
    if (!fp) {
        return 0;
    }
    fclose(fp);
    return 1;
}
