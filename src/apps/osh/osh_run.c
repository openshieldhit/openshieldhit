#include "apps/osh/osh_run.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_app_osh.h"
#include "beam/osh_beamdef.h"
#include "beam/runtime/osh_beam_runtime.h"
#include "gemca/osh_gemca2.h"
#include "gemca/runtime/osh_gemca_runtime.h"
#include "material/runtime/osh_material_prepare.h"
#include "openshieldhit/geometry.h"
#include "openshieldhit/material.h"
#include "openshieldhit/scoring.h"
#include "scoring/runtime/osh_scoring_postprocess.h"
#include "scoring/runtime/osh_scoring_prepare.h"
#include "scoring/save/osh_scoring_save.h"
#include "transport/osh_transport.h"

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
    struct osh_beam_runtime *beam_rt = NULL;
    struct osh_transport_context transport_ctx;
    struct osh_geometry_workspace *geom = NULL;
    struct osh_gemca_runtime geom_rt;
    struct osh_material_workspace *mat = NULL;
    struct osh_scoring_workspace *scoring = NULL;
    struct osh_material_runtime transport_tables;
    struct osh_scoring_runtime scoring_runtime;
    enum osh_status rc = OSH_OK;

    if (!opt) {
        return OSH_EINVAL;
    }

    memset(&transport_tables, 0, sizeof(transport_tables));
    memset(&scoring_runtime, 0, sizeof(scoring_runtime));
    memset(&geom_rt, 0, sizeof(geom_rt));
    memset(&transport_ctx, 0, sizeof(transport_ctx));

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
    if (out) {
        fprintf(out, "Loaded beam: %s\n", beam_path);
        if (opt->has_nstat) {
            fprintf(out, "Applied nstat override: %llu\n", opt->nstat);
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

    /* Resolve zone → material_name strings to dense material indices.
     * The names were stored in the internal gemca workspace during prepare(). */
    {
        struct osh_gemca_prepared *gemca = geom->prepared;
        size_t iz;
        for (iz = 0u; iz < gemca->nzones; ++iz) {
            struct zone *z = gemca->zones[iz];
            struct osh_material const *m;

            if (!z->material_name) {
                if (err) {
                    fprintf(err, "Error: zone '%s' has no material assigned\n", z->name ? z->name : "(unnamed)");
                }
                rc = OSH_EPARSE;
                goto cleanup;
            }
            m = osh_material_by_name(mat, z->material_name);
            if (!m) {
                if (err) {
                    fprintf(err,
                            "Error: zone '%s': unknown material '%s'\n",
                            z->name ? z->name : "(unnamed)",
                            z->material_name);
                }
                rc = OSH_EPARSE;
                goto cleanup;
            }
            z->material_idx = m->index;
        }
    }
    if (out) {
        fprintf(out, "Material assembly complete: %llu zones resolved.\n", (unsigned long long) geom->prepared->nzones);
    }

    if (osh_gemca_runtime_setup(geom->prepared, &geom_rt) != OSH_OK) {
        if (err) {
            fprintf(err, "Error: failed to compile geometry runtime\n");
        }
        rc = OSH_ENOMEM;
        goto cleanup;
    }
    if (out) {
        fprintf(out,
                "Geometry runtime compiled: %llu bodies, %llu zones.\n",
                (unsigned long long) geom_rt.nbodies,
                (unsigned long long) geom_rt.nzones);
    }

    {
        unsigned int z_max = (beam->primary.z > 0u) ? (unsigned int) beam->primary.z : 1u;
        if (osh_material_prepare(mat, z_max, &transport_tables) != OSH_OK) {
            if (err) {
                fprintf(err, "Error: failed to prepare runtime transport tables\n");
            }
            rc = OSH_EPARSE;
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

    if (osh_scoring_prepare(scoring, &scoring_runtime) != OSH_OK) {
        if (err) {
            fprintf(err, "Error: failed to prepare scoring runtime\n");
        }
        rc = OSH_EPARSE;
        goto cleanup;
    }
    if (out) {
        fprintf(out,
                "Scoring runtime prepared: %llu geometries, %llu outputs, %llu pages.\n",
                (unsigned long long) scoring_runtime.ngeometries,
                (unsigned long long) scoring_runtime.noutputs,
                (unsigned long long) scoring_runtime.npages);
    }

    transport_ctx.params.nstat = beam->nstat;
    transport_ctx.params.deltae = beam->deltae;
    transport_ctx.params.demin = beam->demin;
    transport_ctx.params.tcut = beam->tcut;
    transport_ctx.params.rndseed = beam->rndseed;
    transport_ctx.params.rndoffset = beam->rndoffset;
    switch (beam->scatter) {
    case OSH_BEAM_MSCAT_GAUSS:
        transport_ctx.params.mcs_mode = OSH_TRANSPORT_MCS_GAUSSIAN;
        break;
    case OSH_BEAM_MSCAT_MOLIERE:
        transport_ctx.params.mcs_mode = OSH_TRANSPORT_MCS_MOLIERE;
        break;
    default:
        transport_ctx.params.mcs_mode = OSH_TRANSPORT_MCS_OFF;
        break;
    }
    switch (beam->straggl) {
    case OSH_BEAM_STRAGG_GAUSS:
        transport_ctx.params.straggling_mode = OSH_TRANSPORT_STRAGGLING_GAUSSIAN;
        break;
    case OSH_BEAM_STRAGG_VAVILOV:
        transport_ctx.params.straggling_mode = OSH_TRANSPORT_STRAGGLING_VAVILOV;
        break;
    default:
        transport_ctx.params.straggling_mode = OSH_TRANSPORT_STRAGGLING_OFF;
        break;
    }

    if (osh_beam_runtime_setup(beam, &beam_rt) != OSH_OK) {
        if (err) {
            fprintf(err, "Error: failed to initialise beam runtime\n");
        }
        rc = OSH_ESTATE;
        goto cleanup;
    }

    if (osh_transport_run_minimal(&transport_ctx, beam_rt, &geom_rt, &transport_tables, &scoring_runtime) != OSH_OK) {
        if (err) {
            fprintf(err, "Error: transport failed\n");
        }
        rc = OSH_ESTATE;
        goto cleanup;
    }
    if (out) {
        fprintf(out, "Transport completed: %llu primaries.\n", (unsigned long long) transport_ctx.params.nstat);
    }

    if (osh_scoring_postprocess(&scoring_runtime) != OSH_OK) {
        if (err) {
            fprintf(err, "Error: scoring postprocess failed\n");
        }
        rc = OSH_ESTATE;
        goto cleanup;
    }
    if (out) {
        fprintf(out, "Scoring postprocess completed.\n");
    }

    {
        struct osh_scoring_save_request save_req;

        memset(&save_req, 0, sizeof(save_req));
        save_req.out_dir = outdir;
        save_req.ws = scoring;
        save_req.rt = &scoring_runtime;
        save_req.nstat = beam->nstat;
        save_req.has_nstat = 1;

        if (osh_scoring_save(&save_req) != OSH_OK) {
            if (err) {
                fprintf(err, "Error: scoring save failed\n");
            }
            rc = OSH_ESTATE;
            goto cleanup;
        }
    }
    if (out) {
        fprintf(out, "Scoring outputs saved to %s\n", outdir);
        fprintf(out, "Run completed.\n");
    }

cleanup:
    osh_scoring_runtime_free(&scoring_runtime);
    if (scoring) {
        osh_scoring_workspace_free(scoring);
    }
    osh_gemca_runtime_free(&geom_rt);
    osh_material_runtime_free(&transport_tables);
    if (mat) {
        osh_material_workspace_free(mat);
    }
    osh_beam_runtime_free(beam_rt);
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
