#ifndef OSH_SIGNALS_H
#define OSH_SIGNALS_H

#include <signal.h> /* sig_atomic_t */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Install a graceful-stop handler for interactive interruption.
 *
 * @details
 * Maps an interrupt request from the OS onto a single async-signal-safe flag
 * (see @ref osh_signals_stop_flag).  On POSIX this catches @c SIGINT via
 * @c sigaction; on Windows it registers a @c SetConsoleCtrlHandler so Ctrl-C
 * (and console close/logoff) set the same flag.  On platforms without either
 * mechanism the call is a no-op and the flag simply stays clear, so the run
 * behaves exactly as if no handler were installed.
 *
 * The handler does nothing but raise the flag — no I/O, no allocation — which
 * keeps it async-signal-safe.  The transport layer polls the flag at its safe
 * points and stops cleanly (in-flight histories finish, secondaries drain),
 * rather than aborting mid-history.
 *
 * Idempotent; safe to call once at startup.
 */
void osh_signals_install_stop(void);

/**
 * @brief Borrow the process-wide graceful-stop flag.
 *
 * @details
 * Returns a pointer to the flag raised by the installed handler.  Pass it to
 * @c osh_simulation_set_run_control() so the run can observe interrupt
 * requests.  The pointer is always valid (it refers to a static); it just never
 * changes from 0 when no handler was installed.
 *
 * @returns Address of the volatile stop flag (non-NULL).
 */
sig_atomic_t volatile *osh_signals_stop_flag(void);

#ifdef __cplusplus
}
#endif

#endif /* OSH_SIGNALS_H */
