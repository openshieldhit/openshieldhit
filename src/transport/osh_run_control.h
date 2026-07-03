#ifndef OSH_RUN_CONTROL_H
#define OSH_RUN_CONTROL_H

#include <stddef.h> /* size_t */

#ifdef __cplusplus
extern "C" {
#endif

/* Pointer-only references (G1/G2): the run-control block carries a dump
 * destination and its reusable postprocess scratch without transport/ having to
 * include the scoring save/snapshot headers.  Only the .c files that actually
 * fire a dump pull in scoring/runtime/osh_scoring_snapshot.h. */
struct osh_scoring_sink;
struct osh_scoring_shadow;

/**
 * @brief Portable run-control transport — issue #170 / #192 (PR-2 of 3).
 *
 * @details
 * Carries the *intent* knobs that let a driver stop a run cleanly or, later,
 * dump periodic previews — without the transport layer ever learning about
 * files, signals, or the host OS.  The external stop source is a generic
 * *callback* (goal G3 — "how the stop is requested is platform glue"): the
 * library only asks @c should_stop(user) at safe points and never names a
 * signal, thread, or browser primitive.  The app supplies whatever backs that
 * answer — a POSIX signal flag, a Windows console handler, a watchdog thread, a
 * GUI button, or a WASM @c postMessage — so the same transport code is portable.
 * This mirrors the @c fn + @c void* idiom already used by the scoring sink (G1).
 *
 * @b PR-3 scope (issue #193): the dump-cadence scalars (@ref dump_every_s,
 * @ref dump_every_primaries) and the symmetric on-demand @ref should_dump
 * callback (G3) now drive @ref run_ctl_should_dump, and the dump destination
 * (@ref dump_sink + reusable @ref dump_shadow scratch, G1) plus the output
 * selector (@ref dump_outputs, G2) let a driver fire a non-destructive snapshot
 * without the transport layer knowing about files, signals, or the host OS.
 *
 * @b Quiescence (issue #195): both a clean stop and a periodic dump are only
 * *physically meaningful* at a family-complete point.  This module never
 * short-circuits the family scheduler; it only tells the *primary* (ion) loop to
 * stop injecting new primaries and drain what is in flight.  The scheduler then
 * runs the remaining family passes (neutrons, …) over the secondaries those
 * completed primaries banked, so a partial result is *family-exact* for exactly
 * the primaries that finished — never an ion-only fraction under a full
 * @c completed_nstat.  Consequently the driver evaluates @ref run_ctl_should_dump
 * at the *checkpoint boundary* (family-quiescent), not the inner ion safe point,
 * so every dump this module drives is EXACT.
 */
struct osh_run_control {
    double wall_budget_s;           /**< Wall-clock budget [s]; 0 = unlimited. */
    double dump_every_s;            /**< Wall-time dump cadence [s]; 0 = off. */
    size_t dump_every_primaries;    /**< Primary-count dump cadence; 0 = off. */
    int (*should_stop)(void *user); /**< Borrowed; returns non-zero to stop cleanly.  NULL = never stop
                                         early.  Must be cheap and non-blocking — called at safe points. */
    void *should_stop_user;         /**< Opaque context handed back to @ref should_stop. */
    int (*should_dump)(void *user); /**< Borrowed on-demand trigger (G3); returns non-zero to request a
                                         one-off dump (e.g. POSIX SIGUSR1).  NULL = no on-demand dumps.
                                         Expected to be edge-triggered (read-and-clear) so each external
                                         signal fires exactly one dump.  Cheap and non-blocking. */
    void *should_dump_user;         /**< Opaque context handed back to @ref should_dump. */

    /* Dump destination (borrowed; set by the driver, not the transport).  All
     * three are needed together for a dump to fire; any NULL disables dumping. */
    struct osh_scoring_sink const *dump_sink; /**< Where a snapshot is written (G1: file / buffer / …). */
    struct osh_scoring_shadow *dump_shadow;   /**< Reusable non-destructive postprocess scratch bound to the
                                                   live scoring runtime; mutated per dump, never the live acc. */
    size_t const *dump_outputs;               /**< Output selector (G2): indices to dump; NULL = all outputs. */
    size_t dump_noutputs;                     /**< Number of entries in @ref dump_outputs. */

    /* Internal run-lifetime state (set by run_ctl_start, not by callers).  The
     * dump bookkeeping is kept in run-relative seconds (elapsed since start), the
     * same clock @ref run_ctl_should_dump is passed, so a dump cadence needs no
     * absolute-time baseline of its own. */
    double t_start;             /**< Monotonic timestamp of run start [s]. */
    double last_dump_s;         /**< Elapsed seconds at the last dump. */
    size_t last_dump_primaries; /**< Completed-primary count at the last dump. */
};

/**
 * @brief Zero-initialise a run-control block to "no limits, no early stop".
 *
 * @param[out] ctl  Block to initialise; no-op when NULL.
 */
void osh_run_control_init(struct osh_run_control *ctl);

/**
 * @brief Mark the start of the run (records the monotonic baseline @p t_now).
 *
 * @details
 * The wall-time budget is measured from this baseline, so the driver should
 * call it once, immediately before transport begins, with a timestamp from
 * @c osh_monotonic_seconds().  Resets the dump bookkeeping to the same instant.
 *
 * @param[in,out] ctl    Block to arm; no-op when NULL.
 * @param[in]     t_now  Monotonic start timestamp [s].
 */
void osh_run_control_start(struct osh_run_control *ctl, double t_now);

/**
 * @brief Decide whether the run should stop injecting new primaries.
 *
 * @details
 * @c should_stop = (wall_budget_s > 0 && elapsed >= wall_budget_s)
 *                  || (should_stop && should_stop(should_stop_user)).
 *
 * A NULL @p ctl never stops (the un-controlled baseline).  @p completed is
 * accepted for a stable signature and future count-based policies; it is not
 * used by the current rule.
 *
 * @param[in] ctl        Run-control block, or NULL.
 * @param[in] elapsed    Seconds since @ref osh_run_control_start.
 * @param[in] completed  Primaries completed so far (currently unused).
 *
 * @returns 1 to stop cleanly, 0 to keep going.
 */
int run_ctl_should_stop(const struct osh_run_control *ctl, double elapsed, size_t completed);

/**
 * @brief Decide whether a partial-result dump is due, and record that it fired.
 *
 * @details
 * Returns 1 when any of the three triggers is active, then advances the dump
 * bookkeeping so the next call measures from *this* dump:
 *   - **on-demand**: @c should_dump(@c should_dump_user) returns non-zero (a
 *     SIGUSR1-style one-off request; the callback is expected to read-and-clear);
 *   - **time cadence**: @c dump_every_s elapsed since the last dump
 *     (@p elapsed − @c last_dump_s ≥ @c dump_every_s);
 *   - **count cadence**: @c dump_every_primaries reached since the last dump
 *     (@p completed − @c last_dump_primaries ≥ @c dump_every_primaries).
 *
 * A NULL @p ctl never dumps.  The caller is expected to evaluate this at a
 * family-complete checkpoint boundary (issue #195) so the dump it triggers is
 * physically exact; the function itself is pure bookkeeping and does not enforce
 * where it is called.  Unlike @ref run_ctl_should_stop this takes a **mutable**
 * @p ctl because a positive answer updates @c last_dump_s / @c last_dump_primaries.
 *
 * @param[in,out] ctl        Run-control block, or NULL.
 * @param[in]     elapsed    Seconds since @ref osh_run_control_start.
 * @param[in]     completed  Primaries completed so far (family-exact at a checkpoint).
 *
 * @returns 1 when a dump should fire now, 0 otherwise.
 */
int run_ctl_should_dump(struct osh_run_control *ctl, double elapsed, size_t completed);

/**
 * @brief Is a *scheduled* (cadence-driven) dump configured on this control?
 *
 * @details
 * True when a time or count cadence is set — i.e. dumps will fire on a schedule
 * regardless of any external signal.  This is the discriminator for the shadow
 * memory budget-reservation rule (issue #193): a scheduled dump budgets its
 * shadow up front (it *will* happen, so its cost is accounted before the run
 * rather than discovered mid-run), whereas an on-demand-only dump (@ref
 * should_dump set but no cadence) allocates lazily and is fail-soft.  The shadow
 * itself still allocates lazily at the first dump, so the reservation lowers risk
 * rather than being an absolute guarantee.  An on-demand callback alone does
 * **not** count as scheduled.
 *
 * @param[in] ctl  Run-control block, or NULL.
 * @returns 1 when @c dump_every_s > 0 or @c dump_every_primaries > 0, else 0.
 */
int run_ctl_has_scheduled_dump(struct osh_run_control const *ctl);

#ifdef __cplusplus
}
#endif

#endif /* OSH_RUN_CONTROL_H */
