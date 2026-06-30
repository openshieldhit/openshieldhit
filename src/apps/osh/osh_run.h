#ifndef OSH_RUN_H
#define OSH_RUN_H

#include <stdio.h>

#include "openshieldhit/diag.h"
#include "openshieldhit/status.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Options controlling a single OSH simulation run.
 *
 * @details
 * All path pointers are borrowed (not owned); the caller keeps them alive
 * for the duration of the osh_run() call.  NULL values fall back to
 * defaults (workdir defaults to ".", other paths to workdir/<name>.dat).
 */
struct osh_run_options {
    char const *workdir;              /**< Working directory for resolving default file names. */
    char const *out_dir;              /**< Output directory; NULL falls back to workdir. */
    char const *geo_path;             /**< Explicit geometry file path; NULL → workdir/geo.dat. */
    char const *beam_path;            /**< Explicit beam file path;     NULL → workdir/beam.dat. */
    char const *mat_path;             /**< Explicit material file path;  NULL → workdir/mat.dat. */
    char const *detect_path;          /**< Explicit detect file path;    NULL → workdir/detect.dat. */
    unsigned long long nstat;         /**< Primary history count override; used only when has_nstat != 0. */
    int has_nstat;                    /**< 1 if nstat should override the beam file value. */
    unsigned long long seed_offset;   /**< RNG seed offset override; used only when has_seed_offset != 0. */
    int has_seed_offset;              /**< 1 if seed_offset should override the beam file value. */
    int validate_only;                /**< 1 = validate inputs then exit without running transport. */
    unsigned long long pool_capacity; /**< Transport pool capacity override; used only when has_pool_capacity != 0. */
    int has_pool_capacity;            /**< 1 if pool_capacity should override the compiled default. */
    char const *mem_budget;           /**< Memory-budget override string (e.g. "8GB", "80%"); NULL = default policy. */
    char const *profile_path;         /**< Profile JSON output path; NULL disables profiling. */
    double max_time_s;                /**< Wall-time budget [s]; used only when has_max_time != 0. */
    int has_max_time;                 /**< 1 if max_time_s should override the beam file MAXTIME card. */
    int (*should_stop)(void *user);   /**< Borrowed graceful-stop callback (e.g. SIGINT-backed); NULL = none. */
    void *should_stop_user;           /**< Opaque context passed to @ref should_stop. */
    struct osh_diag_sink const *diag; /**< Borrowed diagnostics sink for simulation/transport messages. */
};

/**
 * @brief Run a complete OSH simulation.
 *
 * @details
 * Loads all inputs, runs the transport loop, postprocesses scoring, and
 * saves outputs.  Informational messages are written to @p out; error
 * messages are written to @p err.  Either may be NULL to suppress output.
 *
 * The caller is responsible for initialising the logger before calling
 * this function.
 *
 * @param[in] opt  Run options (must not be NULL).
 * @param[in] out  Stream for progress output, or NULL.
 * @param[in] err  Stream for error messages, or NULL.
 *
 * @returns OSH_OK on success, or an OSH_E* code on failure.
 */
enum osh_status osh_run(struct osh_run_options const *opt, FILE *out, FILE *err);

#ifdef __cplusplus
}
#endif

#endif /* OSH_RUN_H */
