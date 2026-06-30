#ifndef OSH_RUN_CONTROL_H
#define OSH_RUN_CONTROL_H

#include <signal.h> /* sig_atomic_t */
#include <stddef.h> /* size_t */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Portable run-control transport — issue #170 / #192 (PR-2 of 3).
 *
 * @details
 * Carries the *intent* knobs that let a driver stop a run cleanly or, later,
 * dump periodic previews — without the transport layer ever learning about
 * files, signals, or the host OS.  The control interface is deliberately
 * *flags*, not signals (goal G3): a borrowed @c stop_flag the app raises
 * however it likes (POSIX @c sigaction, a Windows console handler, a WASM
 * message, a watchdog thread), so the same transport code is portable.
 *
 * @b PR-2 scope (this struct, what is consumed now): @ref wall_budget_s and
 * @ref stop_flag drive @ref run_ctl_should_stop.  The dump-cadence fields are
 * declared here so the struct shape is stable, but they are *not* read by the
 * transport yet — periodic / on-demand dumps land in PR-3 (issue #193), gated
 * on the checkpoint/quiescence contract worked out in issue #195.
 *
 * @b Quiescence (issue #195): a clean stop is only *physically meaningful* at a
 * family-complete point.  This module never short-circuits the family
 * scheduler; it only tells the *primary* (ion) loop to stop injecting new
 * primaries and drain what is in flight.  The scheduler then runs the remaining
 * family passes (neutrons, …) over the secondaries those completed primaries
 * banked, so the partial result is *family-exact* for exactly the primaries
 * that finished — never an ion-only fraction under a full @c completed_nstat.
 */
struct osh_run_control {
    double wall_budget_s;             /**< Wall-clock budget [s]; 0 = unlimited. */
    double dump_every_s;              /**< Wall-time dump cadence [s]; 0 = off (consumed in PR-3). */
    size_t dump_every_primaries;      /**< Primary-count dump cadence; 0 = off (consumed in PR-3). */
    sig_atomic_t volatile *stop_flag; /**< Borrowed; raised → stop cleanly.  NULL = never stop early. */
    sig_atomic_t volatile *dump_flag; /**< Borrowed; on-demand dump request.  NULL = none (PR-3). */

    /* Internal run-lifetime state (set by run_ctl_start, not by callers). */
    double t_start;             /**< Monotonic timestamp of run start [s]. */
    double last_dump_s;         /**< Monotonic timestamp of the last dump [s] (PR-3). */
    size_t last_dump_primaries; /**< Completed-primary count at the last dump (PR-3). */
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
 *                  || (stop_flag && *stop_flag).
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

#ifdef __cplusplus
}
#endif

#endif /* OSH_RUN_CONTROL_H */
