#ifndef OSH_RUN_CONTROL_H
#define OSH_RUN_CONTROL_H

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
 * files, signals, or the host OS.  The external stop source is a generic
 * *callback* (goal G3 — "how the stop is requested is platform glue"): the
 * library only asks @c should_stop(user) at safe points and never names a
 * signal, thread, or browser primitive.  The app supplies whatever backs that
 * answer — a POSIX signal flag, a Windows console handler, a watchdog thread, a
 * GUI button, or a WASM @c postMessage — so the same transport code is portable.
 * This mirrors the @c fn + @c void* idiom already used by the scoring sink (G1).
 *
 * @b PR-2 scope (this struct, what is consumed now): @ref wall_budget_s and
 * @ref should_stop drive @ref run_ctl_should_stop.  The dump-cadence scalars are
 * declared here so the struct shape is stable, but they are *not* read by the
 * transport yet — periodic / on-demand dumps land in PR-3 (issue #193), which
 * will add a symmetric @c should_dump callback, gated on the checkpoint /
 * quiescence contract worked out in issue #195.
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
    double wall_budget_s;           /**< Wall-clock budget [s]; 0 = unlimited. */
    double dump_every_s;            /**< Wall-time dump cadence [s]; 0 = off (consumed in PR-3). */
    size_t dump_every_primaries;    /**< Primary-count dump cadence; 0 = off (consumed in PR-3). */
    int (*should_stop)(void *user); /**< Borrowed; returns non-zero to stop cleanly.  NULL = never stop
                                         early.  Must be cheap and non-blocking — called at safe points. */
    void *should_stop_user;         /**< Opaque context handed back to @ref should_stop. */

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

#ifdef __cplusplus
}
#endif

#endif /* OSH_RUN_CONTROL_H */
