#ifndef OPENSHIELDHIT_SIMULATION_H
#define OPENSHIELDHIT_SIMULATION_H

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
    double transport_s;                         /**< Total transport wall time [s]. */
    double phase_fill_s;                        /**< Pool refill from the beam source [s]. */
    double phase_zone_ref_s;                    /**< Batched zone-ref lookup [s]. */
    double phase_distance_s;                    /**< Batched boundary-distance query [s]. */
    double phase_step_s;                        /**< Per-particle physics step loop [s]. */
    double phase_compact_s;                     /**< Dead-slot pool compaction [s]. */
    unsigned long long steps;                   /**< Total transport steps taken. */
    unsigned long long iterations;              /**< Wavefront loop iterations. */
    unsigned long long nuclear_events;          /**< Nuclear interactions sampled. */
    unsigned long long secondaries;             /**< Secondaries produced by nuclear events. */
    unsigned long long neutrons_banked;         /**< Neutrons routed to the (untransported) neutron pool. */
    unsigned long long fragments_banked;        /**< Residual fragments banked for future breakup. */
    unsigned long long ion_secondaries_dropped; /**< Ion secondaries lost to pool overflow; 0 in the default
                                                     configuration, which reserves secondary headroom (issue #213). */
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
 * The control interface is deliberately a *callback*, not an OS signal or a
 * typed flag, so it is portable: the library only asks @p should_stop(@p user)
 * at safe points and never names a signal, thread, or browser primitive.  The
 * caller backs that answer however it likes — a POSIX signal flag, a Windows
 * console handler, a watchdog thread, a GUI button, or a WASM message — and the
 * run stops cleanly the next time the callback returns non-zero.  @p user is an
 * opaque context handed straight back to @p should_stop (the usual C "callback +
 * context" idiom; mirrors the scoring sink).
 *
 * "Cleanly" means the transport stops *injecting* new primaries but lets every
 * in-flight history finish, then drains all banked secondary families, so the
 * partial result is family-exact for exactly the primaries that completed.
 * osh_results_completed_nstat() then reports that true count and every output is
 * normalised by it.  Passing wall_budget_s <= 0 and should_stop == NULL disables
 * the policy (identical to not calling this function).
 *
 * @p should_stop must be cheap and non-blocking: it is polled once per wavefront
 * safe point, not once per history.
 *
 * @param[in] sim            Simulation handle created by osh_simulation_create().
 * @param[in] wall_budget_s  Wall-clock budget [s]; <= 0 means unlimited.
 * @param[in] should_stop    Borrowed callback; returns non-zero to stop cleanly.
 *                           NULL = no external stop source.
 * @param[in] user           Opaque context passed to @p should_stop; may be NULL.
 *
 * @returns OSH_OK on success, OSH_EINVAL when @p sim is NULL.
 */
enum osh_status osh_simulation_set_run_control(struct osh_simulation *sim,
                                               double wall_budget_s,
                                               int (*should_stop)(void *user),
                                               void *user);

/**
 * @brief Configure family-complete checkpoint batching for the next run (issue #195).
 *
 * @details
 * Must be called after osh_simulation_create() and before osh_simulation_run().
 * Selects how often the run is brought to a quiescent, *family-complete*
 * checkpoint — a point where every transport family (ions, then the
 * neutrons/fragments they banked, …) has been drained into scoring, so a result
 * observed there is physically complete rather than an ion-only fraction.
 *
 * @p every_primaries is a **count cadence**: the run is transported in
 * family-complete batches of up to @p every_primaries primaries each, and every
 * batch boundary is a checkpoint.  Count cadence is deterministic and
 * order-independent (each history's RNG stream is a pure function of its global
 * index), which makes it the reproducible cadence for tests and CI.
 *
 *   - @p every_primaries == 0 → **FINAL-ONLY** (default): one batch of K = nstat,
 *     the fastest path, byte-for-byte identical to a run that never calls this.
 *   - @p every_primaries  > 0 → **LIVE**: family-complete batches of that size.
 *     Scored output matches the final-only result up to floating-point reduction
 *     order.
 *
 * This lands the batch-aware seam and the quiescence guarantee only; the machinery
 * that hangs off a checkpoint — periodic file dumps and the time cadence (#193),
 * variance-batch folding (#169), and per-worker accumulator merges (#161) — is
 * added by those follow-ups.
 *
 * @param[in] sim              Simulation handle created by osh_simulation_create().
 * @param[in] every_primaries  Count cadence in primaries; 0 = final-only (default).
 *
 * @returns OSH_OK on success; OSH_EINVAL when @p sim is NULL, or when
 *          @p every_primaries exceeds what the platform's @c size_t can hold
 *          (only possible where @c size_t is narrower than @c unsigned @c long
 *          @c long, e.g. some 32-bit builds).
 */
enum osh_status osh_simulation_set_checkpoint_policy(struct osh_simulation *sim, unsigned long long every_primaries);

/**
 * @brief Configure the sequential score-replica diagnostic harness (issue #230).
 *
 * @details
 * Must be called after osh_simulation_create() and before osh_simulation_run().
 * Splits the run's [0, nstat) histories into @p replicas contiguous sub-ranges,
 * transports each one **sequentially** (no threads) into its own private
 * accumulator set, then merges all of them into the master before the normal
 * postprocess + save.  It is a correctness- and profiling-harness: the first thing
 * that exercises the private-accumulator + merge reduce every parallel backend
 * (threads, MPI, WASM workers) will depend on — with zero concurrency risk.  It
 * does **not** speed anything up.
 *
 * Reproducibility (the deterministic core of the parallel contract, issue #168):
 *   - @p replicas == 0 (default) → the shared-master fast path, byte-for-byte the
 *     un-replica'd run.
 *   - @p replicas == 1 → one private set merged into an empty master: the same
 *     summation order as serial, so **bit-identical** to it.
 *   - @p replicas  > 1 → identical per-history physics (each history's RNG stream
 *     is a pure function of its global index); only cross-partition floating-point
 *     summation order differs, so the result matches serial within tolerance, and
 *     the **same** @p replicas run twice is bit-identical.
 *
 * @param[in] sim       Simulation handle created by osh_simulation_create().
 * @param[in] replicas  Number of sequential replicas; 0 disables the harness.
 *                      Must be <= the run's nstat (each replica needs at least one
 *                      history) or the call is rejected.
 *
 * @returns OSH_OK on success; OSH_EINVAL when @p sim is NULL or @p replicas
 *          exceeds nstat.
 */
enum osh_status osh_simulation_set_score_replicas(struct osh_simulation *sim, size_t replicas);

/**
 * @brief Configure periodic and on-demand partial-result dumps (issue #193).
 *
 * @details
 * Must be called after osh_simulation_create() and before osh_simulation_run().
 * Arranges for the run to write a *non-destructive* snapshot of the current
 * scoring result to the same output files as the final save, refining them on
 * disk as the run progresses.  Each dump is taken at a family-complete checkpoint
 * (issue #195), so the partial result written is physically **exact** — every
 * secondary family the completed primaries banked has been drained into scoring —
 * and the live accumulators are left byte-identical, so the run continues
 * unperturbed.
 *
 * Three independent triggers, evaluated at each checkpoint boundary:
 *   - @p dump_every_s > 0 — a **time** cadence: dump roughly every this many
 *     wall-clock seconds.  The production cadence; its per-run overhead is bounded
 *     independent of core count.  Non-deterministic by nature.
 *   - @p dump_every_primaries > 0 — a **count** cadence: dump every this many
 *     completed primaries.  Deterministic and reproducible (the natural cadence
 *     for tests); equivalent to the beam.dat @c NSTAT save step / @c nsave.
 *   - @p should_dump — an **on-demand** callback (e.g. POSIX @c SIGUSR1) polled at
 *     each checkpoint; expected to be edge-triggered (read-and-clear) so one
 *     external request yields one dump.  NULL for none.
 *
 * Setting a time or count cadence also puts the run in LIVE checkpoint-batching
 * mode (the dump cadence *is* the checkpoint cadence — "a periodic dump and a
 * parallel checkpoint are the same operation", issue #170), overriding any prior
 * osh_simulation_set_checkpoint_policy().  Passing all three off (0, 0, NULL)
 * disables dumping and restores the final-only fast path.
 *
 * @note @p should_dump is observed only at checkpoint boundaries, and those exist
 * only when a **cadence** is also set — the final boundary is deliberately skipped
 * (the run's own end-of-run save already writes the complete result).  An
 * on-demand trigger with no cadence therefore has **no effect**: pair @c SIGUSR1
 * with a @p dump_every_s / @p dump_every_primaries cadence, which both enables
 * dumps and bounds how soon the on-demand request is serviced.
 *
 * A dump is a preview, never the run's product: a write or allocation failure is
 * logged and the run continues to its exact final save rather than aborting.
 *
 * @param[in] sim                  Simulation handle created by osh_simulation_create().
 * @param[in] dump_every_s         Time cadence [s]; <= 0 disables the time trigger.
 * @param[in] dump_every_primaries Count cadence in primaries; 0 disables the count trigger.
 * @param[in] should_dump          Borrowed on-demand callback; NULL = none.  Observed
 *                                 only at checkpoints, so it needs a cadence to fire.
 * @param[in] user                 Opaque context passed to @p should_dump; may be NULL.
 *
 * @returns OSH_OK on success; OSH_EINVAL when @p sim is NULL, or when
 *          @p dump_every_primaries exceeds what the platform's @c size_t can hold.
 */
enum osh_status osh_simulation_set_dump_control(struct osh_simulation *sim,
                                                double dump_every_s,
                                                unsigned long long dump_every_primaries,
                                                int (*should_dump)(void *user),
                                                void *user);

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
