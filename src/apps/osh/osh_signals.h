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
 * @brief Run-control stop callback backed by the graceful-stop flag.
 *
 * @details
 * Ready-made adapter matching @c osh_simulation_set_run_control's
 * @c should_stop signature: returns non-zero once the installed handler has
 * raised the flag.  This is the seam that keeps every signal/OS detail in the
 * app layer — the library only ever calls this through a plain function pointer
 * and never sees @c sig_atomic_t.  @p user is ignored (the flag is a process
 * singleton), so pass NULL.
 *
 * @param[in] user  Unused; pass NULL.
 * @returns 1 if a graceful stop has been requested, 0 otherwise.
 */
int osh_signals_should_stop(void *user);

/**
 * @brief Borrow the process-wide graceful-stop flag.
 *
 * @details
 * Returns a pointer to the flag raised by the installed handler.  Most callers
 * want @ref osh_signals_should_stop instead; this raw accessor exists for
 * tests and introspection.  The pointer is always valid (it refers to a
 * static); it just never changes from 0 when no handler was installed.
 *
 * @returns Address of the volatile stop flag (non-NULL).
 */
sig_atomic_t volatile *osh_signals_stop_flag(void);

#ifdef __cplusplus
}
#endif

#endif /* OSH_SIGNALS_H */
