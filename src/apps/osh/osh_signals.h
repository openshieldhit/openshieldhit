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

#ifdef __cplusplus
}
#endif

#endif /* OSH_SIGNALS_H */
