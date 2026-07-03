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
#include "apps/osh/osh_membudget.h"
#include "common/osh_diag.h"
#include "common/osh_file.h"
#include "common/osh_patient_position.h"
#include "common/osh_sysinfo.h"
#include "common/osh_time.h"
#include "common/osh_vect.h"
#include "openshieldhit/const.h"
#include "openshieldhit/dicom.h"
#include "openshieldhit/geometry.h"
#include "openshieldhit/geometry_defs.h"
#include "openshieldhit/scoring.h"
#include "openshieldhit/simulation.h"
#include "openshieldhit/version.h"

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
static enum osh_status run_check_memory(struct osh_scoring_workspace const *scoring,
                                        char const *mem_budget_override,
                                        int reserve_shadow,
                                        FILE *out,
                                        FILE *err);
static enum osh_status run_write_profile_json(char const *path,
                                              struct osh_simulation_profile const *prof,
                                              size_t nstat,
                                              int rndseed,
                                              int rndoffset,
                                              double parse_s,
                                              double compile_s,
                                              double run_s,
                                              double save_s,
                                              struct osh_diag_sink const *diag);

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
    double t_mark;
    double parse_s = 0.0;
    double compile_s = 0.0;
    double run_s = 0.0;
    double save_s = 0.0;
    double eff_dump_every_s = 0.0;                      /* CLI --dump-every over beam.dat DUMPEVERY */
    unsigned long long eff_dump_every_primaries = 0ull; /* CLI --dump-every-primaries over NSTAT nsave */
    int scheduled_dump = 0;                             /* a periodic dump cadence is active → reserve shadow */
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
    rc = osh_path_ensure_dir(abs_outdir);
    if (rc != OSH_OK) {
        if (err) {
            fprintf(err, "Error: failed to create output directory: %s\n", abs_outdir);
        }
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

    t_mark = osh_monotonic_seconds();

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

    parse_s = osh_monotonic_seconds() - t_mark;

    if (opt->validate_only) {
        if (opt->profile_path && opt->profile_path[0]) {
            OSH_DIAG_WARNF(opt->diag, "%s", "profile requested with --dry-run; no profile written");
        }
        if (out) {
            fprintf(out, "Validation completed.\n");
        }
        goto cleanup;
    }

    /* Resolve the effective dump cadences so the memory check can reserve the
     * snapshot shadow when a periodic dump *will* happen, and
     * osh_simulation_set_dump_control() below can apply the same values.
     *
     * Precedence: a command-line flag always overrides the matching beam.dat card.
     * The time cadence is CLI --dump-every over the DUMPEVERY card; the count
     * cadence is CLI --dump-every-primaries over the NSTAT save step (nsave). */
    if (opt->has_dump_every) {
        eff_dump_every_s = opt->dump_every_s; /* CLI --dump-every wins */
    } else {
        eff_dump_every_s = beam->dump_every_s; /* else the beam.dat DUMPEVERY card */
    }
    if (opt->has_dump_every_primaries) {
        eff_dump_every_primaries = opt->dump_every_primaries; /* CLI --dump-every-primaries wins */
    } else {
        eff_dump_every_primaries = (unsigned long long) beam->nsave; /* else the NSTAT save step */
    }
    /* A dump is "scheduled" when either cadence is active (as opposed to only the
     * on-demand SIGUSR1 trigger); scheduled dumps are the ones that reserve memory. */
    scheduled_dump = 0;
    if (eff_dump_every_s > 0.0 || eff_dump_every_primaries > 0ull) {
        scheduled_dump = 1;
    }

    /* Detect host resources, report the scoring memory footprint, and refuse
     * the run up front if it would exceed the memory budget — before
     * osh_simulation_create() allocates any scoring buffers.  When a periodic dump
     * is scheduled, the snapshot shadow is reserved here too so a scheduled dump
     * can never OOM mid-run (issue #193 budget-reservation rule). */
    rc = run_check_memory(scoring, opt->mem_budget, scheduled_dump, out, err);
    if (rc != OSH_OK) {
        goto cleanup;
    }

    t_mark = osh_monotonic_seconds();
    rc = osh_simulation_create(beam, geom, mat, scoring, opt->diag, &sim);
    compile_s = osh_monotonic_seconds() - t_mark;
    if (rc != OSH_OK) {
        if (err) {
            fprintf(err, "Error: failed to compile simulation\n");
        }
        goto cleanup;
    }

    if (opt->profile_path && opt->profile_path[0]) {
        osh_simulation_set_profiling(sim, 1);
    }

    if (opt->has_pool_capacity) {
        osh_simulation_set_pool_capacity(sim, (size_t) opt->pool_capacity);
        if (out) {
            fprintf(out, "Applied pool capacity override: %llu\n", opt->pool_capacity);
        }
    }

    /* Sequential score-replica diagnostic (issue #230): splits the run into N
     * private-accumulator sub-ranges merged into the master.  Rejected here if it
     * exceeds nstat (the earliest point both are known), before the run starts. */
    if (opt->has_score_replicas) {
        rc = osh_simulation_set_score_replicas(sim, (size_t) opt->score_replicas);
        if (rc != OSH_OK) {
            if (err) {
                fprintf(
                    err, "Error: --score-replicas %llu is invalid (must be >= 1 and <= nstat)\n", opt->score_replicas);
            }
            goto cleanup;
        }
        if (out) {
            fprintf(out,
                    "Score replicas   : %llu sequential private-accumulator sub-range(s) "
                    "(diagnostic; merged before save)\n",
                    opt->score_replicas);
        }
    }

    /* Clean-stop / wall-budget policy.  The CLI --max-time overrides the
     * beam.dat MAXTIME card; the graceful-stop callback (e.g. from Ctrl-C) is
     * wired in unconditionally so an interactive interrupt always stops cleanly. */
    {
        double wall_budget_s = opt->has_max_time ? opt->max_time_s : beam->wall_budget_s;
        if (wall_budget_s > 0.0 || opt->should_stop) {
            osh_simulation_set_run_control(sim, wall_budget_s, opt->should_stop, opt->should_stop_user);
        }
        if (out && wall_budget_s > 0.0) {
            fprintf(out,
                    "Wall-time budget : %g s%s (stop cleanly, save partial result)\n",
                    wall_budget_s,
                    opt->has_max_time ? " (--max-time)" : " (beam.dat MAXTIME)");
        }
    }

    /* Periodic / on-demand partial-result dumps (issue #193).  Activated only when
     * a cadence is set: the run then checkpoints at that cadence and overwrites the
     * output files with the exact partial result each time.  The on-demand SIGUSR1
     * callback is passed through so it can trigger an extra dump at the next
     * checkpoint within a cadenced run. */
    if (scheduled_dump) {
        rc = osh_simulation_set_dump_control(
            sim, eff_dump_every_s, eff_dump_every_primaries, opt->should_dump, opt->should_dump_user);
        if (rc != OSH_OK) {
            if (err) {
                fprintf(err, "Error: failed to configure partial-result dumps\n");
            }
            goto cleanup;
        }
        if (out) {
            if (eff_dump_every_primaries > 0ull) {
                fprintf(out,
                        "Periodic dumps   : every %llu primaries%s (exact partial result)\n",
                        eff_dump_every_primaries,
                        opt->has_dump_every_primaries ? " (--dump-every-primaries)" : " (beam.dat NSTAT step)");
            }
            if (eff_dump_every_s > 0.0) {
                fprintf(out,
                        "Periodic dumps   : every %g s%s (exact partial result)\n",
                        eff_dump_every_s,
                        opt->has_dump_every ? " (--dump-every)" : " (beam.dat DUMPEVERY)");
            }
        }
    }

    t_mark = osh_monotonic_seconds();
    rc = osh_simulation_run(sim);
    run_s = osh_monotonic_seconds() - t_mark;
    if (rc != OSH_OK) {
        if (err) {
            fprintf(err, "Error: simulation run failed\n");
        }
        goto cleanup;
    }

    t_mark = osh_monotonic_seconds();
    rc = osh_simulation_save(sim);
    save_s = osh_monotonic_seconds() - t_mark;
    if (rc != OSH_OK) {
        if (err) {
            fprintf(err, "Error: failed to save scoring outputs\n");
        }
        goto cleanup;
    }
    if (out) {
        struct osh_results const *results = NULL;
        if (osh_simulation_get_results(sim, &results) == OSH_OK && results) {
            unsigned long long const requested = osh_results_requested_nstat(results);
            unsigned long long const completed = osh_results_completed_nstat(results);
            if (completed < requested) {
                fprintf(out,
                        "Stopped early: %llu of %llu primaries completed; "
                        "partial result saved (normalised by %llu).\n",
                        completed,
                        requested,
                        completed);
            }
        }
        fprintf(out, "Run completed. Outputs saved under %s\n", abs_outdir);
    }

    if (opt->profile_path && opt->profile_path[0]) {
        struct osh_simulation_profile prof;
        rc = osh_simulation_get_profile(sim, &prof);
        if (rc == OSH_OK) {
            rc = run_write_profile_json(opt->profile_path,
                                        &prof,
                                        beam->nstat,
                                        beam->rndseed,
                                        beam->rndoffset,
                                        parse_s,
                                        compile_s,
                                        run_s,
                                        save_s,
                                        opt->diag);
        }
        if (rc != OSH_OK) {
            if (err) {
                fprintf(err, "Error: failed to write profile JSON: %s\n", opt->profile_path);
            }
            goto cleanup;
        }
        if (out) {
            fprintf(out, "Profile written to %s\n", opt->profile_path);
        }
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
 * @brief Report host resources and the scoring memory footprint, and gate the
 *        run against a memory budget.
 *
 * @details
 * Detects CPU/RAM via osh_sysinfo, resolves a budget (the @p mem_budget_override
 * string when given, otherwise the default 80%-of-available policy), estimates
 * the scoring accumulator memory from the parsed configuration, and prints a
 * short report to @p out.  If a budget is known and the estimate exceeds it,
 * prints an actionable message to @p err and returns OSH_ENOMEM so the caller
 * aborts *before* any scoring buffers are allocated — turning a would-be
 * out-of-memory crash into a clear, recoverable refusal.
 *
 * When no memory figure can be detected and no override is given, the budget is
 * 0 and the run proceeds unguarded (best effort, never a false refusal).
 *
 * @param[in] scoring             Parsed scoring workspace (fully populated).
 * @param[in] mem_budget_override Budget string (e.g. "8GB", "80%") or NULL.
 * @param[in] reserve_shadow      Non-zero when a periodic dump is scheduled: add
 *                                the snapshot shadow (est.shadow_bytes) to the
 *                                footprint so a scheduled dump can never OOM
 *                                mid-run (issue #193).  On-demand-only dumps do
 *                                not reserve — they allocate lazily and fail-soft.
 * @param[in] out                 Human-readable report stream, or NULL.
 * @param[in] err                 Error stream for the refusal message, or NULL.
 *
 * @returns OSH_OK to proceed, OSH_ENOMEM if over budget, OSH_EINVAL on a
 *          malformed override string.
 */
static enum osh_status run_check_memory(struct osh_scoring_workspace const *scoring,
                                        char const *mem_budget_override,
                                        int reserve_shadow,
                                        FILE *out,
                                        FILE *err) {
    struct osh_sysinfo info;
    struct osh_scoring_mem_estimate est;
    uint64_t base;
    uint64_t budget = 0u;
    uint64_t footprint; /* accumulators + (reserved) snapshot shadow */
    char b_total[32];
    char b_avail[32];
    char b_scoring[32];
    char b_budget[32];
    char b_largest[32];
    char const *dump_drop_hint; /* extra "drop periodic dumps" advice in the refusal message */

    osh_sysinfo_query(&info);

    /* Use available RAM as the reference for "%" budget strings so that
     * "80%" means 80 % of what is actually free, not of what is installed. */
    if (info.ram_available_bytes > 0u) {
        base = info.ram_available_bytes;
    } else {
        base = info.ram_total_bytes; /* fall back when available is unknown */
    }

    if (mem_budget_override && mem_budget_override[0]) {
        if (!osh_membudget_parse(mem_budget_override, base, &budget)) {
            if (err) {
                fprintf(err,
                        "Error: invalid --mem-budget value '%s' (use e.g. 8GB, 512MiB, or 80%%)\n",
                        mem_budget_override);
            }
            return OSH_EINVAL;
        }
    } else {
        budget = osh_membudget_default(&info); /* 0 ⇒ unknown ⇒ no enforcement */
    }

    if (osh_scoring_estimate_memory(scoring, &est) != OSH_OK) {
        return OSH_OK; /* estimate unavailable: do not block the run */
    }

    /* Footprint gated against the budget: the accumulators, plus the snapshot
     * shadow when a periodic dump is scheduled (it will be allocated on the first
     * dump, so reserving it up front means a scheduled dump never OOMs mid-run).
     * Saturate rather than wrap on the (practically impossible) overflow. */
    footprint = est.accum_bytes;
    if (reserve_shadow && est.shadow_bytes > 0u) {
        footprint = (footprint <= UINT64_MAX - est.shadow_bytes) ? (footprint + est.shadow_bytes) : UINT64_MAX;
    }

    osh_sysinfo_format_bytes(info.ram_total_bytes, b_total, sizeof(b_total));
    osh_sysinfo_format_bytes(info.ram_available_bytes, b_avail, sizeof(b_avail));
    osh_sysinfo_format_bytes(footprint, b_scoring, sizeof(b_scoring));

    if (out) {
        if (info.logical_cores > 0u) {
            fprintf(out,
                    "Resources: %u core(s), RAM %s total / %s available\n",
                    info.logical_cores,
                    info.ram_total_bytes > 0u ? b_total : "unknown",
                    info.ram_available_bytes > 0u ? b_avail : "unknown");
        } else {
            fprintf(out,
                    "Resources: cores unknown, RAM %s total / %s available\n",
                    info.ram_total_bytes > 0u ? b_total : "unknown",
                    info.ram_available_bytes > 0u ? b_avail : "unknown");
        }
        if (reserve_shadow && est.shadow_bytes > 0u) {
            char b_shadow[32];
            osh_sysinfo_format_bytes(est.shadow_bytes, b_shadow, sizeof(b_shadow));
            fprintf(out,
                    "Scoring memory: %s across %zu page(s) (includes %s reserved for periodic-dump snapshot)\n",
                    b_scoring,
                    est.npages,
                    b_shadow);
        } else {
            fprintf(out, "Scoring memory: %s across %zu page(s)\n", b_scoring, est.npages);
        }
        if (budget > 0u) {
            osh_sysinfo_format_bytes(budget, b_budget, sizeof(b_budget));
            fprintf(
                out, "Memory budget: %s%s\n", b_budget, mem_budget_override ? " (--mem-budget)" : " (default policy)");
        }
    }

    if (budget > 0u && footprint > budget) {
        if (err) {
            osh_sysinfo_format_bytes(budget, b_budget, sizeof(b_budget));
            osh_sysinfo_format_bytes(est.largest_page_bytes, b_largest, sizeof(b_largest));
            /* Only suggest dropping periodic dumps when they are the reason the
             * footprint grew (a scheduled dump actually reserved a shadow). */
            if (reserve_shadow && est.shadow_bytes > 0u) {
                dump_drop_hint = ", drop periodic dumps";
            } else {
                dump_drop_hint = "";
            }
            fprintf(err,
                    "Error: scoring would allocate %s, exceeding the memory budget of %s.\n"
                    "  Largest scorer: geometry '%s' (%s).\n"
                    "  Detected RAM: %s total, %s available.\n"
                    "  To proceed, reduce the scoring mesh/bin counts%s, or raise the limit with\n"
                    "  --mem-budget (e.g. --mem-budget %s). Aborting before allocating memory.\n",
                    b_scoring,
                    b_budget,
                    est.largest_geometry,
                    b_largest,
                    info.ram_total_bytes > 0u ? b_total : "unknown",
                    info.ram_available_bytes > 0u ? b_avail : "unknown",
                    dump_drop_hint,
                    b_scoring);
        }
        return OSH_ENOMEM;
    }

    return OSH_OK;
}

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
 * Mirrors the _setup_vox() calculation exactly, using the public body arguments:
 *   a[9]  = gantry angle [deg]
 *   a[10] = couch angle [deg]
 *   a[11..13] = universe position of local voxel-corner (x0,y0,z0) [cm]
 *   a[17] = patient position code (enum osh_patient_position; present when na >= 18)
 *
 * The transform chain is identical to _setup_vox():
 *   1. Base rotation from patient position (IEC 61217 DICOM→universe mapping).
 *   2. Couch rotation around universe Z (vertical axis per IEC 61217).
 *   3. Gantry rotation around universe Y (sagittal-plane axis per IEC 61217).
 *
 * This function is used by run_setup_voxel_scoring() to reconstruct the body
 * transform for RTDOSE/DicomCT scoring geometries after parsing is complete.
 *
 * @param[in]  b  VOX body with na >= 14; a[17] read when na >= 18.
 * @param[out] t  Receives the 4x4 row-major transform (same layout as body->t).
 */
static void _vox_body_build_transform(struct osh_geometry_body const *b, double t[16]) {
    double tb[3][3];
    double gantry_rad;
    double couch_rad;
    double tx;
    double ty;
    double tz;
    enum osh_patient_position pp;
    int i;
    int j;

    gantry_rad = b->a[9] * OSH_M_PI_180;
    couch_rad = b->a[10] * OSH_M_PI_180;
    tx = b->a[11];
    ty = b->a[12];
    tz = b->a[13];

    /* Read patient position from a[17] (present when na >= 18); fall back to
     * HFS for old geo.dat files that produced only 17 arguments. */
    pp = (b->na >= 18) ? (enum osh_patient_position)(int) b->a[17] : OSH_PP_HFS;
    osh_patient_position_base_rotation(pp, tb);

    /* Apply couch then gantry rotations (IEC 61217 axis convention).
     * Must match _setup_vox() exactly to keep this reconstruction correct. */
    for (i = 0; i < 3; i++) {
        osh_vect_rot_z(-couch_rad, tb[i]); /* IEC couch: CCW from above; rot_z is CW, so negate */
        osh_vect_rot_y(gantry_rad, tb[i]); /* IEC gantry: around universe Y (cranial-caudal) */
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
 * RTDOSE grid placement (two-stage):
 *  Stage 1 here: a[14..16] = -ct.origin[j]/10 [cm] are added to the RTDOSE
 *   DICOM origin so the dose grid is expressed in the CT body's local
 *   DICOM-aligned frame.  No universe-frame transform is applied in this stage.
 *  Stage 2 at score time: the scoring geometry's has_rotation flag and ct_t
 *   matrix (copied to g->t) cause the scoring step to transform each particle
 *   from universe to CT-local before the voxel bin lookup.
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
    int ct_has_rotation = 0;             /* set when patient position code is present (na >= 18) */
    enum osh_status rc;

    if (!geom || !scoring) {
        return OSH_EINVAL;
    }

    /* Find the CT VOX body: read the CT-origin offset (a[14..16]) and
     * reconstruct the universe->local affine transform (used at score time).
     *
     * WHY ct_has_rotation is always set when na >= 18:
     * The patient-position base rotation is ALWAYS non-trivial (except HFS at
     * zero gantry/couch, which still has a different axis mapping than identity).
     * Previously ct_has_rotation was only set when gantry or couch != 0, but
     * that missed the base rotation.  Now: any DCM body with 18 arguments has
     * a known patient position, so the transform is always needed for correct
     * RTDOSE/DicomCT scoring.
     *
     * TWO-STAGE RTDOSE PLACEMENT:
     * Stage 1 (here): a[14..16] = -ct.origin[j]/10 [cm] are added to the RTDOSE
     *   DICOM absolute origin to give CT-local coordinates (DICOM-relative to CT
     *   first voxel).  lo_x/y/z come out in the CT-local (DICOM) frame.
     * Stage 2 (score time): the scoring geometry has has_rotation=1 and ct_t
     *   copied into g->t.  The scoring step transforms each particle position
     *   from universe to CT-local before the voxel bin lookup.
     * No patient-position rotation is needed in Stage 1 because RTDOSE and CT
     * are both in DICOM physical (LPS) space; the universe<->DICOM rotation is
     * handled entirely by ct_t at scoring time. */
    offset_x = offset_y = offset_z = 0.0;
    for (i = 0; i < geom->nbodies; ++i) {
        if (geom->bodies[i].type == OSH_GEOMETRY_BODY_VOX && geom->bodies[i].na >= 18) {
            offset_x = geom->bodies[i].a[14];
            offset_y = geom->bodies[i].a[15];
            offset_z = geom->bodies[i].a[16];
            _vox_body_build_transform(&geom->bodies[i], ct_t);
            /* Patient position base rotation is always non-trivial: any DCM
             * body with a patient position code always introduces a rotation. */
            ct_has_rotation = 1;
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

        /* Convert DICOM voxel-center origin to voxel-corner [cm] in CT-local frame.
         *
         * offset_{x,y,z} = a[14..16] = -ct.origin[j]/10 [cm].
         * Adding these to the RTDOSE DICOM origin gives the RTDOSE position
         * relative to the CT's first voxel, i.e. CT-local coordinates:
         *   lo_j = (rd.origin[j] - ct.origin[j]) / 10 - half_pixel
         *
         * No universe-frame rotation is applied here because both the RTDOSE
         * and the CT voxel grid are in DICOM physical (LPS) space.  The
         * has_rotation flag and ct_t transform handle the universe->DICOM
         * mapping at scoring time. */
        lo_x = rd.origin[0] / 10.0 + offset_x - 0.5 * dx;
        lo_y = rd.origin[1] / 10.0 + offset_y - 0.5 * dy;
        lo_z = rd.origin[2] / 10.0 + rd.frame_offsets[0] / 10.0 + offset_z - 0.5 * dz;

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
 * @brief Identify the compiler that produced this binary.
 *
 * @details
 * Profiling numbers are meaningless without their toolchain; this string is
 * embedded in the profile JSON so result files remain self-describing.
 */
static char const *run_compiler_string(void) {
#if defined(__clang__)
    return "clang " __clang_version__;
#elif defined(__GNUC__)
    return "gcc " __VERSION__;
#elif defined(_MSC_VER)
    return "msvc";
#else
    return "unknown";
#endif
}

/**
 * @brief Write a one-line JSON profile record for a completed run.
 *
 * @details
 * The record is self-contained: binary metadata (version, compiler), run
 * parameters (nstat, seed), the app-level phase wall times (parse, compile,
 * run, save), the transport-loop phase decomposition, and event counters.
 * The benchmark harness in tools/bench/ wraps this record with case and
 * machine metadata; the binary emits only what it can know itself.
 *
 * @returns OSH_OK on success, OSH_EIO when the file cannot be written.
 */
static enum osh_status run_write_profile_json(char const *path,
                                              struct osh_simulation_profile const *prof,
                                              size_t nstat,
                                              int rndseed,
                                              int rndoffset,
                                              double parse_s,
                                              double compile_s,
                                              double run_s,
                                              double save_s,
                                              struct osh_diag_sink const *diag) {
    FILE *fp;
    double const phase_sum_s = prof->phase_fill_s + prof->phase_zone_ref_s + prof->phase_distance_s + prof->phase_step_s
                               + prof->phase_compact_s;
    double const prim_per_s = (prof->transport_s > 0.0) ? ((double) nstat / prof->transport_s) : 0.0;
    double const steps_per_primary = (nstat > 0u) ? ((double) prof->steps / (double) nstat) : 0.0;

#if defined(_MSC_VER)
    if (fopen_s(&fp, path, "w") != 0) {
        fp = NULL;
    }
#else
    fp = fopen(path, "w");
#endif
    if (!fp) {
        OSH_DIAG_ERRORF(diag, "profile: cannot open '%s' for writing", path);
        return OSH_EIO;
    }

    fprintf(fp,
            "{\"schema\":1,\"version\":\"%s\",\"compiler\":\"%s\","
            "\"nstat\":%zu,\"rndseed\":%d,\"rndoffset\":%d,"
            "\"setup_parse_s\":%.9g,\"setup_compile_s\":%.9g,"
            "\"run_s\":%.9g,\"transport_s\":%.9g,\"save_s\":%.9g,"
            "\"prim_per_s\":%.9g,"
            "\"phases\":{\"fill_s\":%.9g,\"zone_ref_s\":%.9g,\"distance_s\":%.9g,"
            "\"step_s\":%.9g,\"compact_s\":%.9g,\"sum_s\":%.9g},"
            "\"counters\":{\"steps\":%llu,\"steps_per_primary\":%.9g,"
            "\"iterations\":%llu,\"nuclear_events\":%llu,\"secondaries\":%llu,"
            "\"neutrons_banked\":%llu,\"fragments_banked\":%llu,"
            "\"ion_secondaries_dropped\":%llu}}\n",
            osh_version_string(),
            run_compiler_string(),
            nstat,
            rndseed,
            rndoffset,
            parse_s,
            compile_s,
            run_s,
            prof->transport_s,
            save_s,
            prim_per_s,
            prof->phase_fill_s,
            prof->phase_zone_ref_s,
            prof->phase_distance_s,
            prof->phase_step_s,
            prof->phase_compact_s,
            phase_sum_s,
            prof->steps,
            steps_per_primary,
            prof->iterations,
            prof->nuclear_events,
            prof->secondaries,
            prof->neutrons_banked,
            prof->fragments_banked,
            prof->ion_secondaries_dropped);

    if (fclose(fp) != 0) {
        OSH_DIAG_ERRORF(diag, "profile: failed to finalize '%s'", path);
        return OSH_EIO;
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
