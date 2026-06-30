#ifndef OPENSHIELDHIT_SIMULATION_H
#define OPENSHIELDHIT_SIMULATION_H

#include <signal.h> /* sig_atomic_t */

#include "openshieldhit/beam.h"
#include "openshieldhit/geometry.h"
#include "openshieldhit/material.h"
#include "openshieldhit/scoring.h"
#include "openshieldhit/status.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Opaque handle for a compiled, ready-to-run simulation.
 *
 * @details
 * Holds all four runtime representations compiled from the cold workspaces.
 * The cold workspaces are borrowed (not owned); the caller must keep them
 * alive for the lifetime of the simulation object.
 */
struct osh_simulation;
struct osh_results;

/**
 * @brief Compile four cold workspaces into a simulation ready to run.
 *
 * @details
 * Performs zone-to-material index resolution, compiles geometry and scoring
 * runtimes, compiles transport tables, and initialises the beam source.
 *
 * @param[in]  beam     Prepared beam workspace.
 * @param[in]  geo      Prepared geometry workspace.
 * @param[in]  mat      Finalized material workspace.
 * @param[in]  scoring  Parsed scoring workspace.
 * @param[in]  diag     Borrowed diagnostics sink for run/save messages; NULL for silent.
 * @param[out] sim_out  Receives the new simulation handle on success.
 *
 * @returns OSH_OK on success, or an error code.
 */
enum osh_status osh_simulation_create(struct osh_beam_workspace *beam,
                                      struct osh_geometry_workspace *geo,
                                      struct osh_material_workspace *mat,
                                      struct osh_scoring_workspace *scoring,
                                      struct osh_diag_sink const *diag,
                                      struct osh_simulation **sim_out);

/**
 * @brief Transport-phase wall-clock timers and event counters for one run.
 *
 * @details
 * Filled by osh_simulation_get_profile() after a completed run when profiling
 * was enabled with osh_simulation_set_profiling().  The five phase timers
 * decompose transport_s; their sum is slightly below transport_s because
 * progress reporting and loop bookkeeping are not attributed to any phase.
 *
 * Profiling reads only the monotonic clock and pre-existing counters; it never
 * touches the RNG streams or physics state, so profiled runs produce
 * bit-identical scoring output to unprofiled ones.
 */
struct osh_simulation_profile {
    double transport_s;                  /**< Total transport wall time [s]. */
    double phase_fill_s;                 /**< Pool refill from the beam source [s]. */
    double phase_zone_ref_s;             /**< Batched zone-ref lookup [s]. */
    double phase_distance_s;             /**< Batched boundary-distance query [s]. */
    double phase_step_s;                 /**< Per-particle physics step loop [s]. */
    double phase_compact_s;              /**< Dead-slot pool compaction [s]. */
    unsigned long long steps;            /**< Total transport steps taken. */
    unsigned long long iterations;       /**< Wavefront loop iterations. */
    unsigned long long nuclear_events;   /**< Nuclear interactions sampled. */
    unsigned long long secondaries;      /**< Secondaries produced by nuclear events. */
    unsigned long long neutrons_banked;  /**< Neutrons routed to the (untransported) neutron pool. */
    unsigned long long fragments_banked; /**< Residual fragments banked for future breakup. */
};

/**
 * @brief Enable or disable transport profiling for subsequent runs.
 *
 * @details
 * Must be called after osh_simulation_create() and before
 * osh_simulation_run().  When enabled, the transport loop accumulates phase
 * timers and counters retrievable via osh_simulation_get_profile().  The
 * off-path cost is a single pointer test per loop phase; results are
 * bit-identical either way.
 *
 * @param[in] sim     Simulation handle created by osh_simulation_create().
 * @param[in] enable  Non-zero to enable, zero to disable.
 *
 * @returns OSH_OK on success, OSH_EINVAL when @p sim is NULL.
 */
enum osh_status osh_simulation_set_profiling(struct osh_simulation *sim, int enable);

/**
 * @brief Set the live-history pool capacity for the next run.
 *
 * @details
 * Must be called after osh_simulation_create() and before
 * osh_simulation_run().  The pool capacity is the number of particle
 * histories transported simultaneously; it trades cache footprint (small)
 * against batch parallelism (large).  It is a pure performance knob: because
 * each history owns an independent RNG stream keyed by its global index, every
 * history consumes the same random draws at any capacity, so scored results
 * are independent of the capacity chosen up to floating-point summation order
 * in the shared scoring accumulators.
 *
 * @param[in] sim       Simulation handle created by osh_simulation_create().
 * @param[in] capacity  Number of simultaneous histories; 0 selects the
 *                      compiled default.  Values above nstat are clamped.
 *
 * @returns OSH_OK on success, OSH_EINVAL when @p sim is NULL, OSH_ENOMEM
 *          when pool reallocation fails (simulation is then unusable and must
 *          be freed with osh_simulation_free()).
 */
enum osh_status osh_simulation_set_pool_capacity(struct osh_simulation *sim, size_t capacity);

/**
 * @brief Install a clean-stop / wall-time budget policy for the next run.
 *
 * @details
 * Must be called after osh_simulation_create() and before osh_simulation_run().
 * The control interface is deliberately a *flag*, not an OS signal, so it is
 * portable: the caller raises @p stop_flag from a POSIX signal handler, a
 * Windows console handler, a watchdog thread, or a WASM message, and the run
 * stops cleanly at the next safe point.
 *
 * "Cleanly" means the transport stops *injecting* new primaries but lets every
 * in-flight history finish, then drains all banked secondary families, so the
 * partial result is family-exact for exactly the primaries that completed.
 * osh_results_completed_nstat() then reports that true count and every output is
 * normalised by it.  Passing wall_budget_s <= 0 and stop_flag == NULL disables
 * the policy (identical to not calling this function).
 *
 * @param[in] sim            Simulation handle created by osh_simulation_create().
 * @param[in] wall_budget_s  Wall-clock budget [s]; <= 0 means unlimited.
 * @param[in] stop_flag      Borrowed flag; raised → stop cleanly.  NULL = none.
 *
 * @returns OSH_OK on success, OSH_EINVAL when @p sim is NULL.
 */
enum osh_status
osh_simulation_set_run_control(struct osh_simulation *sim, double wall_budget_s, sig_atomic_t volatile *stop_flag);

/**
 * @brief Retrieve the transport profile of the last completed run.
 *
 * @param[in]  sim  Simulation handle created by osh_simulation_create().
 * @param[out] out  Receives a copy of the profile.
 *
 * @returns OSH_OK on success, OSH_EINVAL on NULL arguments, OSH_ESTATE when
 *          profiling was not enabled before the run.
 */
enum osh_status osh_simulation_get_profile(struct osh_simulation const *sim, struct osh_simulation_profile *out);

/**
 * @brief Run the simulation transport and finalize scoring accumulators.
 *
 * @details
 * Drives the transport loop and postprocesses scoring accumulators. The actual
 * number of simulated primaries is stored internally and used by
 * osh_simulation_save(). Saving is a separate explicit step.
 *
 * @param[in] sim  Simulation handle created by osh_simulation_create().
 *
 * @returns OSH_OK on success, or an error code.
 */
enum osh_status osh_simulation_run(struct osh_simulation *sim);

/**
 * @brief Borrow a read-only handle to the current compiled scoring results.
 *
 * @details
 * The returned handle is owned by @p sim and remains valid until the
 * simulation is freed. Its contents reflect the current scoring runtime state:
 * immediately after osh_simulation_create() that means zeroed accumulators;
 * after osh_simulation_run() it means finalized postprocessed results.
 *
 * @param[in]  sim  Simulation handle created by osh_simulation_create().
 * @param[out] out  Receives the borrowed results handle on success.
 *
 * @returns OSH_OK on success, or an error code.
 */
enum osh_status osh_simulation_get_results(struct osh_simulation const *sim, struct osh_results const **out);

/**
 * @brief Return the requested primary count associated with a results handle.
 *
 * @details
 * This is the target number of primaries requested when the simulation was
 * created. It may differ from the completed primary count when future chunked
 * or time-limited runs are introduced.
 *
 * @param[in] results  Results handle returned by osh_simulation_get_results().
 *
 * @returns Requested primary count, or 0 if @p results is NULL.
 */
unsigned long long osh_results_requested_nstat(struct osh_results const *results);

/**
 * @brief Return the completed primary count currently represented by results.
 *
 * @details
 * Before the first completed run this returns 0. After a successful run it
 * returns the actual number of transported primaries represented by the
 * current results snapshot.
 *
 * @param[in] results  Results handle returned by osh_simulation_get_results().
 *
 * @returns Completed primary count, or 0 if @p results is NULL.
 */
unsigned long long osh_results_completed_nstat(struct osh_results const *results);

/**
 * @brief Return whether the results handle contains completed run data.
 *
 * @details
 * This can be used by callers to distinguish "simulation created but not yet
 * run" from a completed run that happened to transport zero primaries.
 *
 * @param[in] results  Results handle returned by osh_simulation_get_results().
 *
 * @returns 1 if completed run data are available, 0 otherwise.
 */
int osh_results_has_completed_run(struct osh_results const *results);

/**
 * @brief Save all scoring outputs for a finished simulation.
 *
 * @details
 * Iterates all outputs defined in the scoring workspace and writes each one in
 * its configured format (BDO or ASCII). The completed primary count stored by
 * osh_simulation_run() is embedded in BDO files and used to normalise ASCII
 * output per primary.
 *
 * Saving before a completed run is invalid and returns OSH_ESTATE.
 *
 * Output file paths come from the scoring workspace and are expected to be
 * fully resolved by the calling application before the simulation is saved.
 *
 * @param[in] sim  Simulation handle created by osh_simulation_create().
 *
 * @returns OSH_OK on success, OSH_ESTATE if called before a completed run,
 *          OSH_ENOTSUP if a configured format is not supported, or another
 *          OSH_E* on I/O error.
 */
enum osh_status osh_simulation_save(struct osh_simulation const *sim);

/**
 * @brief Release the simulation.
 *
 * @details
 * Frees all runtime resources owned by the simulation.  The cold workspaces
 * passed to osh_simulation_create() are not freed.
 *
 * @param[in] sim  May be NULL (no-op).
 *
 * @returns OSH_OK always.
 */
enum osh_status osh_simulation_free(struct osh_simulation *sim);

#ifdef __cplusplus
}
#endif

#endif /* OPENSHIELDHIT_SIMULATION_H */
