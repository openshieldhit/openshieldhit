#include "apps/osh/osh_run.h"

#if defined(_WIN32)
#include <direct.h>
#define getcwd _getcwd
#else
#include <unistd.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_app_osh.h"
#include "common/osh_diag.h"
#include "common/osh_file.h"
#include "openshieldhit/geometry.h"
#include "openshieldhit/geometry_defs.h"
#include "openshieldhit/scoring.h"
#include "openshieldhit/simulation.h"

/* ---- Default file names -------------------------------------------------- */

static char const *const OSH_DEFAULT_WORKDIR = ".";
static char const *const OSH_GEO_FILENAME = "geo.dat";
static char const *const OSH_BEAM_FILENAME = "beam.dat";
static char const *const OSH_MAT_FILENAME = "mat.dat";
static char const *const OSH_DETECT_FILENAME = "detect.dat";

/* ---- Internal helpers ---------------------------------------------------- */

static char *run_resolve_path(char const *workdir, char const *override_path, char const *filename);
static char *run_resolve_absolute_path(char const *path);
static enum osh_status run_resolve_output_paths(struct osh_scoring_workspace *scoring, char const *out_dir);
static enum osh_status run_setup_voxel_scoring(struct osh_geometry_workspace const *geom,
                                               struct osh_scoring_workspace *scoring,
                                               struct osh_diag_sink const *diag);
static int run_file_exists(char const *path);

/* ---- Run ----------------------------------------------------------------- */

/**
 * @brief Top-level app orchestration for one CLI run or validation pass.
 *
 * @details
 * This is the policy boundary between the file-oriented app layer and the
 * pure simulation API. The function resolves default input names, loads the
 * four cold workspaces via `setup_from_path` helpers, rewrites all scoring
 * output names to fully resolved paths owned by the scoring workspace, then
 * calls the simulation API in explicit phases:
 *
 *   create -> run -> save -> free
 *
 * The library never decides where files are written; that policy stays here.
 * The library does decide which concrete save writer to use for each scoring
 * output block based on the parsed format keyword.
 *
 * @param[in] out  Optional human-readable progress stream.
 * @param[in] err  Optional human-readable error stream.
 */
enum osh_status osh_run(struct osh_run_options const *opt, FILE *out, FILE *err) {
    char const *workdir;
    char const *outdir;
    char *geo_path = NULL;
    char *beam_path = NULL;
    char *mat_path = NULL;
    char *detect_path = NULL;
    char *abs_outdir = NULL;
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
    abs_outdir = run_resolve_absolute_path(outdir);

    if (!geo_path || !beam_path || !mat_path || !detect_path) {
        if (err) {
            fprintf(err, "Error: out of memory resolving input paths\n");
        }
        rc = OSH_ENOMEM;
        goto cleanup;
    }
    if (!abs_outdir) {
        if (err) {
            fprintf(err, "Error: failed to resolve output directory (getcwd failed or invalid path)\n");
        }
        rc = OSH_EIO;
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
        fprintf(out, "  Output directory : %s\n", abs_outdir);
        fprintf(out, "  Geometry input   : %s\n", geo_path);
        fprintf(out, "  Beam input       : %s\n", beam_path);
        fprintf(out, "  Material input   : %s\n", mat_path);
        fprintf(out, "  Detect input     : %s\n", detect_path);
    }

    if (!run_file_exists(beam_path)) {
        if (err) {
            fprintf(err, "Error: beam file not found: %s\n", beam_path);
        }
        rc = OSH_EIO;
        goto cleanup;
    }

    if (osh_beam_setup_from_path(beam_path, opt->diag, &beam) != OSH_OK) {
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

    if (osh_geometry_setup_from_path(geo_path, opt->diag, &geom) != OSH_OK) {
        if (err) {
            fprintf(err, "Error: failed to load geometry: %s\n", geo_path);
        }
        rc = OSH_EPARSE;
        goto cleanup;
    }
    if (out) {
        fprintf(out, "Loaded geometry: %s\n", geo_path);
    }

    if (!run_file_exists(mat_path)) {
        if (err) {
            fprintf(err, "Error: material file not found: %s\n", mat_path);
        }
        rc = OSH_EIO;
        goto cleanup;
    }

    if (osh_material_setup_from_path(mat_path, opt->diag, &mat) != OSH_OK) {
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

    if (osh_scoring_setup_from_path(detect_path, opt->diag, &scoring) != OSH_OK) {
        if (err) {
            fprintf(err, "Error: failed to load scoring/detect input: %s\n", detect_path);
        }
        rc = OSH_EPARSE;
        goto cleanup;
    }
    if (out) {
        fprintf(out, "Loaded scoring: %s\n", detect_path);
    }

    rc = run_resolve_output_paths(scoring, abs_outdir);
    if (rc != OSH_OK) {
        if (err) {
            fprintf(err, "Error: failed to resolve scoring output paths\n");
        }
        goto cleanup;
    }

    rc = run_setup_voxel_scoring(geom, scoring, opt->diag);
    if (rc != OSH_OK) {
        if (err) {
            fprintf(err, "Error: voxel scoring setup failed\n");
        }
        goto cleanup;
    }

    if (opt->validate_only) {
        if (out) {
            fprintf(out, "Validation completed.\n");
        }
        goto cleanup;
    }

    rc = osh_simulation_create(beam, geom, mat, scoring, opt->diag, &sim);
    if (rc != OSH_OK) {
        if (err) {
            fprintf(err, "Error: failed to compile simulation\n");
        }
        goto cleanup;
    }

    rc = osh_simulation_run(sim);
    if (rc != OSH_OK) {
        if (err) {
            fprintf(err, "Error: simulation run failed\n");
        }
        goto cleanup;
    }
    rc = osh_simulation_save(sim);
    if (rc != OSH_OK) {
        if (err) {
            fprintf(err, "Error: failed to save scoring outputs\n");
        }
        goto cleanup;
    }
    if (out) {
        fprintf(out, "Run completed. Outputs saved under %s\n", abs_outdir);
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
    free(abs_outdir);
    return rc;
}

/* ---- Internal helpers ---------------------------------------------------- */

/**
 * @brief Resolve one optional input override against the working directory.
 *
 * @details
 * If @p override_path is present it is copied as-is. Otherwise a default file
 * name is appended to @p workdir. The returned string is heap-allocated and
 * owned by the caller.
 */
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

/**
 * @brief Resolve a path to a normalized absolute path string.
 *
 * @details
 * The CLI may pass a relative output directory. The app layer resolves that
 * once against the process working directory so later save code can treat all
 * scoring output filenames as already-final paths.
 *
 * @returns Newly allocated normalized path, or NULL on failure.
 */
static char *run_resolve_absolute_path(char const *path) {
    char cwd[4096];
    char *resolved = NULL;

    if (!path || !path[0]) {
        return NULL;
    }
    if (!getcwd(cwd, sizeof(cwd))) {
        return NULL;
    }
    if (osh_relative_path_to_file(&resolved, cwd, path) != 0) {
        return NULL;
    }
    osh_path_normalize(resolved);
    return resolved;
}

/**
 * @brief Resolve voxel scoring geometry against the parsed geometry workspace.
 *
 * @details
 * For each DicomCT or DicomRTDOSE scoring geometry found in @p scoring, this
 * function verifies that exactly one CT (VOX) body exists in @p geom, then
 * reads the grid dimensions from its raw argument array and writes
 * @c vox_nbins = nx*ny*nz onto the cold scoring geometry definition so that
 * the subsequent scoring compile step can allocate the correct data arrays.
 *
 * If no voxel scoring geometries are present this function is a no-op.
 * If more or fewer than one CT body is found when a voxel scoring geometry
 * exists, an error is returned.
 *
 * @note
 * The VOX body arg layout is: a[0..2]=origin, a[3..5]=spacing, a[6..8]=nx,ny,nz.
 * Grid population (transform matrix, raytrace grid) is deferred to the M6
 * scoring hot-path implementation.
 */
static enum osh_status run_setup_voxel_scoring(struct osh_geometry_workspace const *geom,
                                               struct osh_scoring_workspace *scoring,
                                               struct osh_diag_sink const *diag) {
    struct osh_geometry_body const *b;
    struct osh_scoring_geometry_def *g;
    char const *kind;
    size_t i;
    size_t nvox_geo = 0u;
    size_t nvox_body = 0u;
    size_t vox_body_idx = 0u;
    size_t nx, ny, nz;

    if (!geom || !scoring) {
        return OSH_EINVAL;
    }

    for (i = 0; i < scoring->ngeometries; ++i) {
        kind = scoring->geometries[i].kind;
        if (kind && (strcmp(kind, "dicomct") == 0 || strcmp(kind, "dicomrtdose") == 0)) {
            ++nvox_geo;
        }
    }

    if (nvox_geo == 0u) {
        return OSH_OK;
    }

    for (i = 0; i < geom->nbodies; ++i) {
        if (geom->bodies[i].type == OSH_GEOMETRY_BODY_VOX) {
            vox_body_idx = i;
            ++nvox_body;
        }
    }

    if (nvox_body != 1u) {
        OSH_DIAG_ERRORF(diag, "voxel scoring requires exactly one CT body in the geometry; found %zu", nvox_body);
        return OSH_EPARSE;
    }

    b = &geom->bodies[vox_body_idx];
    if (b->na < 9) {
        OSH_DIAG_ERRORF(diag, "CT body '%s' has too few arguments to read grid dimensions", b->name);
        return OSH_EPARSE;
    }
    nx = (size_t) b->a[6];
    ny = (size_t) b->a[7];
    nz = (size_t) b->a[8];

    for (i = 0; i < scoring->ngeometries; ++i) {
        kind = scoring->geometries[i].kind;
        if (!kind || (strcmp(kind, "dicomct") != 0 && strcmp(kind, "dicomrtdose") != 0)) {
            continue;
        }
        g = &scoring->geometries[i];
        g->vox_origin[0] = b->a[0];
        g->vox_origin[1] = b->a[1];
        g->vox_origin[2] = b->a[2];
        g->vox_spacing[0] = b->a[3];
        g->vox_spacing[1] = b->a[4];
        g->vox_spacing[2] = b->a[5];
        g->vox_nx = nx;
        g->vox_ny = ny;
        g->vox_nz = nz;
    }

    return OSH_OK;
}

/**
 * @brief Rewrite scoring output filenames to full paths under @p out_dir.
 *
 * @details
 * Each `detect.dat` output keeps its own file name, but the library save
 * layer now expects that name to already be a resolved destination path.
 * This helper performs that rewrite in-place on the scoring workspace by
 * replacing each owned `filename` string with a newly allocated full path.
 */
static enum osh_status run_resolve_output_paths(struct osh_scoring_workspace *scoring, char const *out_dir) {
    size_t i;

    if (!scoring || !out_dir) {
        return OSH_EINVAL;
    }

    for (i = 0; i < scoring->noutputs; ++i) {
        char *resolved = NULL;

        if (!scoring->outputs[i].filename || !scoring->outputs[i].filename[0]) {
            return OSH_EPARSE;
        }
        if (osh_relative_path_to_file(&resolved, out_dir, scoring->outputs[i].filename) != 0) {
            return OSH_ENOMEM;
        }
        osh_path_normalize(resolved);
        free(scoring->outputs[i].filename);
        scoring->outputs[i].filename = resolved;
    }

    return OSH_OK;
}

/**
 * @brief Lightweight existence check used before parser-specific setup.
 */
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
