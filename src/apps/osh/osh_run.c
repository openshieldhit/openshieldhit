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
#include "common/osh_vect.h"
#include "openshieldhit/const.h"
#include "openshieldhit/dicom.h"
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
                                               char const *detect_path,
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

    rc = run_setup_voxel_scoring(geom, scoring, detect_path, opt->diag);
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
 * @brief Append one axis to a scoring geometry definition.
 */
static enum osh_status
run_geo_append_axis(struct osh_scoring_geometry_def *geo, char const *label, double lo, double hi, int nbins) {
    size_t len;
    struct osh_scoring_axis_def *tmp =
        (struct osh_scoring_axis_def *) realloc(geo->axes, (geo->naxes + 1u) * sizeof(*tmp));
    if (!tmp) {
        return OSH_ENOMEM;
    }
    geo->axes = tmp;
    memset(&geo->axes[geo->naxes], 0, sizeof(*tmp));
    len = strlen(label);
    if (len >= sizeof(geo->axes[0].label)) {
        len = sizeof(geo->axes[0].label) - 1u;
    }
    memcpy(geo->axes[geo->naxes].label, label, len);
    geo->axes[geo->naxes].label[len] = '\0';
    geo->axes[geo->naxes].lo = lo;
    geo->axes[geo->naxes].hi = hi;
    geo->axes[geo->naxes].nbins = nbins;
    geo->naxes++;
    return OSH_OK;
}

/**
 * @brief Reconstruct the universe→local affine transform for a VOX body.
 *
 * @details
 * Mirrors the _setup_vox() calculation using the public body arguments:
 *   a[9]  = gantry angle [deg]
 *   a[10] = couch angle [deg]
 *   a[11..13] = universe position of local voxel-corner (x0,y0,z0) [cm]
 *
 * @param[in]  b  VOX body with na >= 14.
 * @param[out] t  Receives the 4×4 row-major transform (same layout as body->t).
 */
static void _vox_body_build_transform(struct osh_geometry_body const *b, double t[16]) {
    double tb[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    double gantry_rad = b->a[9] * OSH_M_PI_180;
    double couch_rad = b->a[10] * OSH_M_PI_180;
    double tx = b->a[11];
    double ty = b->a[12];
    double tz = b->a[13];
    int i;
    int j;

    for (i = 0; i < 3; i++) {
        osh_vect_rot_y(couch_rad, tb[i]);
        osh_vect_rot_z(gantry_rad, tb[i]);
    }
    for (j = 0; j < 3; j++) {
        for (i = 0; i < 3; i++) {
            t[j * 4 + i] = tb[j][i];
        }
        t[j * 4 + 3] = tx * tb[j][0] + ty * tb[j][1] + tz * tb[j][2];
    }
}

/**
 * @brief Resolve voxel scoring geometry against the parsed geometry workspace.
 *
 * @details
 * DicomCT and DicomRTDOSE scoring geometries are converted to Mesh geometries
 * whose axes are expressed in the CT body's local (BZALIGN) frame, and the
 * CT body's universe→local affine transform is stored on the scoring geometry.
 * The scoring step then transforms each particle position to local frame before
 * the raytrace bin lookup, so the result is correct for any gantry/couch angle.
 *
 * Plain Mesh geometries (has_rotation == 0) are unaffected; their axes and
 * particle coordinates are both in universe frame.
 *
 * If no DicomCT or DicomRTDOSE geometries are present this function is a no-op.
 */
static enum osh_status run_setup_voxel_scoring(struct osh_geometry_workspace const *geom,
                                               struct osh_scoring_workspace *scoring,
                                               char const *detect_path,
                                               struct osh_diag_sink const *diag) {
    struct osh_geometry_body const *b;  /* VOX/DCM body from geo.dat */
    struct osh_scoring_geometry_def *g; /* current scoring geometry being mutated */
    struct osh_dicom_rtdose rd;         /* RTDOSE metadata read from InputPath */
    char const *kind;
    char *detect_dir = NULL; /* directory of detect.dat; base for relative InputPath */
    char *resolved = NULL;   /* absolute InputPath after joining with detect_dir */
    char *new_kind = NULL;   /* "mesh" string allocated before overwriting g->kind */
    size_t i;
    size_t nct_geo = 0u;                 /* number of DicomCT scoring geometries found */
    size_t nvox_body = 0u;               /* number of VOX bodies in the transport geometry */
    size_t vox_body_idx = 0u;            /* index of the single VOX body (Phase 2) */
    size_t nx, ny, nz;                   /* RTDOSE grid dimensions: cols, rows, frames */
    double dx, dy, dz;                   /* voxel spacing [cm]: col, row, frame */
    double lo_x, lo_y, lo_z;             /* corner of the first voxel [cm] */
    double offset_x, offset_y, offset_z; /* patient→world offset from the CT body [cm] */
    double ct_t[16];                     /* universe→local transform reconstructed from CT body */
    int ct_has_rotation = 0;             /* set when gantry or couch angle is non-zero */
    int k;
    enum osh_status rc;

    if (!geom || !scoring) {
        return OSH_EINVAL;
    }

    /* Find the CT VOX body: grab patient→world offset (a[14..16]) and
     * reconstruct the universe→local affine transform from gantry/couch/tx.
     * Both are used by the DicomRTDOSE and DicomCT scoring paths below.
     * ct_has_rotation is set only when gantry or couch angle is non-zero;
     * for pure translations the scoring axes stay in universe frame so that
     * BDO output (which cannot represent a rotated frame) keeps working. */
    offset_x = offset_y = offset_z = 0.0;
    for (i = 0; i < geom->nbodies; ++i) {
        if (geom->bodies[i].type == OSH_GEOMETRY_BODY_VOX && geom->bodies[i].na >= 17) {
            offset_x = geom->bodies[i].a[14];
            offset_y = geom->bodies[i].a[15];
            offset_z = geom->bodies[i].a[16];
            _vox_body_build_transform(&geom->bodies[i], ct_t);
            if (geom->bodies[i].a[9] != 0.0 || geom->bodies[i].a[10] != 0.0) {
                ct_has_rotation = 1;
            }
            break;
        }
    }

    /* --- Phase 1: DicomRTDOSE → Mesh conversion --------------------------- */

    for (i = 0; i < scoring->ngeometries; ++i) {
        kind = scoring->geometries[i].kind;
        if (!kind || strcmp(kind, "dicomrtdose") != 0) {
            continue;
        }
        g = &scoring->geometries[i];

        if (!g->vox_rtdose_path || !g->vox_rtdose_path[0]) {
            OSH_DIAG_ERRORF(
                diag, "scoring geometry '%s': DicomRTDOSE requires InputPath", g->name ? g->name : "(unnamed)");
            return OSH_EPARSE;
        }

        /* Resolve InputPath relative to detect.dat directory. */
        detect_dir = detect_path ? osh_path_dirname(detect_path) : NULL;
        if (!detect_dir) {
            detect_dir = strdup(".");
            if (!detect_dir) {
                return OSH_ENOMEM;
            }
        }

        if (osh_relative_path_to_file(&resolved, detect_dir, g->vox_rtdose_path) != 0) {
            free(detect_dir);
            return OSH_ENOMEM;
        }
        osh_path_normalize(resolved);
        free(detect_dir);
        detect_dir = NULL;

        /* Update vox_rtdose_path to the resolved absolute path so the save
         * step can open it regardless of the process working directory. */
        free(g->vox_rtdose_path);
        g->vox_rtdose_path = resolved;
        resolved = NULL;

        rc = osh_dicom_rtdose_read(g->vox_rtdose_path, &rd, diag);
        if (rc != OSH_OK) {
            OSH_DIAG_ERRORF(
                diag, "scoring geometry '%s': failed to read RTDOSE DICOM", g->name ? g->name : "(unnamed)");
            return rc;
        }

        if (rd.rows < 1 || rd.cols < 1 || rd.n_frames < 1 || !rd.frame_offsets) {
            OSH_DIAG_ERRORF(
                diag, "scoring geometry '%s': RTDOSE has invalid dimensions", g->name ? g->name : "(unnamed)");
            osh_dicom_rtdose_free(&rd);
            return OSH_EPARSE;
        }

        nx = (size_t) rd.cols;
        ny = (size_t) rd.rows;
        nz = (size_t) rd.n_frames;
        dx = rd.pixel_spacing[1] / 10.0; /* col spacing mm→cm */
        dy = rd.pixel_spacing[0] / 10.0; /* row spacing mm→cm */
        dz = (nz > 1u) ? (rd.frame_offsets[1] - rd.frame_offsets[0]) / 10.0 : dx;

        if (!(dx > 0.0) || !(dy > 0.0) || !(dz > 0.0)) {
            OSH_DIAG_ERRORF(
                diag, "scoring geometry '%s': RTDOSE has non-positive voxel spacing", g->name ? g->name : "(unnamed)");
            osh_dicom_rtdose_free(&rd);
            return OSH_EPARSE;
        }

        /* Validate uniform z-spacing. */
        if (nz > 2u) {
            size_t iz;
            double expected_dz = rd.frame_offsets[1] - rd.frame_offsets[0];
            for (iz = 2u; iz < nz; ++iz) {
                double actual = rd.frame_offsets[iz] - rd.frame_offsets[iz - 1u];
                if (actual < expected_dz * 0.999 || actual > expected_dz * 1.001) {
                    OSH_DIAG_WARNF(diag,
                                   "scoring geometry '%s': RTDOSE z-spacing is non-uniform; using first interval",
                                   g->name ? g->name : "(unnamed)");
                    break;
                }
            }
        }

        /* DICOM origin is the center of the first voxel; convert to corner [cm].
         * Apply patient→world offset so the RTDOSE grid is in simulation coords.
         * If the CT body is rotated, also transform the corner to CT local frame. */
        lo_x = rd.origin[0] / 10.0 + offset_x - 0.5 * dx;
        lo_y = rd.origin[1] / 10.0 + offset_y - 0.5 * dy;
        lo_z = rd.origin[2] / 10.0 + rd.frame_offsets[0] / 10.0 + offset_z - 0.5 * dz;

        if (ct_has_rotation) {
            double lo_uni[3];
            double lo_loc[3];
            lo_uni[0] = lo_x;
            lo_uni[1] = lo_y;
            lo_uni[2] = lo_z;
            for (k = 0; k < 3; k++) {
                int row = k * 4;
                lo_loc[k] =
                    lo_uni[0] * ct_t[row] + lo_uni[1] * ct_t[row + 1] + lo_uni[2] * ct_t[row + 2] - ct_t[row + 3];
            }
            lo_x = lo_loc[0];
            lo_y = lo_loc[1];
            lo_z = lo_loc[2];
        }

        osh_dicom_rtdose_free(&rd);

        rc = run_geo_append_axis(g, "X", lo_x, lo_x + (double) nx * dx, (int) nx);
        if (rc == OSH_OK) {
            rc = run_geo_append_axis(g, "Y", lo_y, lo_y + (double) ny * dy, (int) ny);
        }
        if (rc == OSH_OK) {
            rc = run_geo_append_axis(g, "Z", lo_z, lo_z + (double) nz * dz, (int) nz);
        }
        if (rc != OSH_OK) {
            return rc;
        }

        if (ct_has_rotation) {
            memcpy(g->t, ct_t, sizeof(g->t));
            g->has_rotation = 1;
        }

        new_kind = strdup("mesh");
        if (!new_kind) {
            return OSH_ENOMEM;
        }
        free(g->kind);
        g->kind = new_kind;
        new_kind = NULL;
        /* vox_rtdose_path is preserved for the save step's RTDOSE writer. */
    }

    /* --- Phase 2: DicomCT → Mesh conversion -------------------------------- */

    for (i = 0; i < scoring->ngeometries; ++i) {
        kind = scoring->geometries[i].kind;
        if (kind && strcmp(kind, "dicomct") == 0) {
            ++nct_geo;
        }
    }

    if (nct_geo == 0u) {
        return OSH_OK;
    }

    for (i = 0; i < geom->nbodies; ++i) {
        if (geom->bodies[i].type == OSH_GEOMETRY_BODY_VOX) {
            vox_body_idx = i;
            ++nvox_body;
        }
    }

    if (nvox_body != 1u) {
        OSH_DIAG_ERRORF(diag, "DicomCT scoring requires exactly one CT body in the geometry; found %zu", nvox_body);
        return OSH_EPARSE;
    }

    b = &geom->bodies[vox_body_idx];
    if (b->na < 14) {
        OSH_DIAG_ERRORF(diag, "CT body '%s' has too few arguments to read grid dimensions", b->name);
        return OSH_EPARSE;
    }

    /* Axes and transform strategy depends on whether the CT body is rotated:
     *  - Rotation present: axes in CT local frame [b->a[0..2]], has_rotation=1.
     *    Scoring step transforms particle to local frame before raytrace.
     *  - No rotation (pure translation, gantry=couch=0): axes stay in universe
     *    frame [b->a[11..13]], has_rotation=0.  Scoring uses universe coords
     *    directly, and BDO output (which cannot encode a rotated frame) works.
     * b->a[3..5]  = dx, dy, dz [cm];  b->a[6..8]  = nx, ny, nz */
    for (i = 0; i < scoring->ngeometries; ++i) {
        kind = scoring->geometries[i].kind;
        if (!kind || strcmp(kind, "dicomct") != 0) {
            continue;
        }
        g = &scoring->geometries[i];

        if (ct_has_rotation) {
            rc = run_geo_append_axis(g, "X", b->a[0], b->a[0] + b->a[6] * b->a[3], (int) b->a[6]);
            if (rc == OSH_OK) {
                rc = run_geo_append_axis(g, "Y", b->a[1], b->a[1] + b->a[7] * b->a[4], (int) b->a[7]);
            }
            if (rc == OSH_OK) {
                rc = run_geo_append_axis(g, "Z", b->a[2], b->a[2] + b->a[8] * b->a[5], (int) b->a[8]);
            }
        } else {
            rc = run_geo_append_axis(g, "X", b->a[11], b->a[11] + b->a[6] * b->a[3], (int) b->a[6]);
            if (rc == OSH_OK) {
                rc = run_geo_append_axis(g, "Y", b->a[12], b->a[12] + b->a[7] * b->a[4], (int) b->a[7]);
            }
            if (rc == OSH_OK) {
                rc = run_geo_append_axis(g, "Z", b->a[13], b->a[13] + b->a[8] * b->a[5], (int) b->a[8]);
            }
        }
        if (rc != OSH_OK) {
            return rc;
        }

        if (ct_has_rotation) {
            memcpy(g->t, ct_t, sizeof(g->t));
            g->has_rotation = 1;
        }

        new_kind = strdup("mesh");
        if (!new_kind) {
            return OSH_ENOMEM;
        }
        free(g->kind);
        g->kind = new_kind;
        new_kind = NULL;
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
