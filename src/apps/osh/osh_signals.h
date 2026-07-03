#ifndef OSH_SIGNALS_H
#define OSH_SIGNALS_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Install a graceful-stop handler for interactive interruption.
 *
 * @details
 * Maps an interrupt request from the OS onto a single process-wide stop flag,
 * polled via @ref osh_signals_should_stop.  On POSIX this catches @c SIGINT via
 * @c sigaction; on Windows it registers a @c SetConsoleCtrlHandler so Ctrl-C
 * (and console close/logoff) set the same flag.  On platforms without either
 * mechanism the call is a no-op and the flag simply stays clear, so the run
 * behaves exactly as if no handler were installed.
 *
 * The handler does nothing but raise the flag — no I/O, no allocation — and the
 * flag is updated with the right concurrency primitive for each platform (a
 * @c sig_atomic_t on POSIX, where the handler runs in the interrupted thread;
 * an interlocked @c LONG on Windows, where the console handler runs on its own
 * thread).  The transport layer polls the flag at its safe points and stops
 * cleanly (in-flight histories finish, secondaries drain), never mid-history.
 *
 * Idempotent; safe to call once at startup.
 */
void osh_signals_install_stop(void);

/**
 * @brief Run-control stop callback backed by the graceful-stop flag.
 *
 * @details
 * Ready-made adapter matching @c osh_simulation_set_run_control's
 * @c should_stop signature: returns non-zero once the installed handler has
 * raised the flag.  This is the seam that keeps every signal/OS detail in the
 * app layer — the library only ever calls this through a plain function pointer
 * and never sees the underlying flag type.  @p user is ignored (the flag is a
 * process singleton), so pass NULL.
 *
 * @param[in] user  Unused; pass NULL.
 * @returns 1 if a graceful stop has been requested, 0 otherwise.
 */
int osh_signals_should_stop(void *user);

/**
 * @brief Clear the graceful-stop flag.
 *
 * @details
 * Resets the request so @ref osh_signals_should_stop reports 0 again.  The
 * application never needs this (a real run sets the flag once and exits); it
 * exists so tests can exercise the set/clear cycle without reaching into the
 * platform-specific flag directly.
 */
void osh_signals_reset_stop(void);

/**
 * @brief Install an on-demand dump handler for interactive, mid-run snapshots.
 *
 * @details
 * On POSIX this catches @c SIGUSR1 via @c sigaction and raises a one-shot dump
 * flag, polled and cleared by @ref osh_signals_should_dump.  The handler only
 * raises the flag — no I/O, no allocation — and the run loop performs the actual
 * dump at its next family-complete checkpoint (issue #193/#195), so the partial
 * result written is physically exact and the live accumulators are untouched.
 *
 * There is no Windows equivalent to @c SIGUSR1, so on Windows (and any platform
 * without the signal) this is a no-op and @ref osh_signals_should_dump always
 * reports 0; scheduled @c --dump-every[-primaries] dumps still work everywhere.
 *
 * Idempotent; safe to call once at startup.
 */
void osh_signals_install_dump(void);

/**
 * @brief Run-control dump callback backed by the on-demand dump flag.
 *
 * @details
 * Ready-made adapter matching @c osh_simulation_set_dump_control's
 * @c should_dump signature.  Unlike @ref osh_signals_should_stop this is
 * **edge-triggered**: it reads *and clears* the flag, returning non-zero exactly
 * once per @c SIGUSR1 so a single signal produces a single dump rather than a
 * dump at every subsequent checkpoint.  @p user is ignored (the flag is a
 * process singleton); pass NULL.
 *
 * @param[in] user  Unused; pass NULL.
 * @returns 1 if a dump was requested since the last call, 0 otherwise.
 */
int osh_signals_should_dump(void *user);

/**
 * @brief Clear the on-demand dump flag.
 *
 * @details
 * Resets any pending request so @ref osh_signals_should_dump reports 0.  Exists
 * for tests; a real run relies on the read-and-clear behaviour of the poll.
 */
void osh_signals_reset_dump(void);

#ifdef __cplusplus
}
#endif

#endif /* OSH_SIGNALS_H */
