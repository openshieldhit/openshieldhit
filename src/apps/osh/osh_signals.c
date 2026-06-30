#include "apps/osh/osh_signals.h"

#include <stddef.h> /* NULL */

/*
 * Graceful-stop flag.  The OS handler only ever raises it; the run loop polls it
 * via osh_signals_should_stop() at safe points.  The flag's type and the way it
 * is read/written differ by platform because the two handler mechanisms have
 * different concurrency models — see each branch.
 */

#if defined(_WIN32)

#include <windows.h>

/* A console control handler runs on a *separate OS thread*, not the main
 * thread, so a plain `sig_atomic_t` (which the C standard only blesses for
 * signal-handler <-> main-flow concurrency) is not enough.  Use a LONG written
 * and read through Interlocked* so the handler-thread store is atomic and
 * ordered with respect to the polling thread. */
static LONG volatile g_stop = 0;

static BOOL WINAPI on_console_ctrl(DWORD ctrl_type) {
    switch (ctrl_type) {
    case CTRL_C_EVENT:
    case CTRL_BREAK_EVENT:
    case CTRL_CLOSE_EVENT:
    case CTRL_LOGOFF_EVENT:
    case CTRL_SHUTDOWN_EVENT:
        InterlockedExchange(&g_stop, 1);
        return TRUE; /* handled: do not run the default terminator */
    default:
        return FALSE;
    }
}

void osh_signals_install_stop(void) {
    /* Returns 0 on failure; we deliberately ignore it — a missing handler just
     * means Ctrl-C keeps its default behaviour, never a crash. */
    (void) SetConsoleCtrlHandler(on_console_ctrl, TRUE);
}

int osh_signals_should_stop(void *user) {
    (void) user; /* the flag is a process singleton; no per-call context needed */
    /* Read with a full barrier so the handler-thread write is observed. */
    return InterlockedCompareExchange(&g_stop, 0, 0) != 0 ? 1 : 0;
}

void osh_signals_reset_stop(void) {
    InterlockedExchange(&g_stop, 0);
}

#else /* POSIX */

#include <signal.h>
#include <string.h> /* memset */

/* The only writer is the SIGINT handler, which runs in the context of the
 * interrupted thread (not a separate one), so `volatile sig_atomic_t` — the
 * type the C standard guarantees safe between a signal handler and the main
 * flow — is exactly right, and writing the flag is async-signal-safe. */
static sig_atomic_t volatile g_stop = 0;

static void on_sigint(int sig) {
    (void) sig;
    g_stop = 1; /* async-signal-safe: flag only */
}

void osh_signals_install_stop(void) {
    struct sigaction sa;
    /* Zero every field first, including implementation-defined ones such as
     * sa_restorer and any padding, so nothing uninitialised reaches the libc
     * wrapper; sa_handler/sa_mask/sa_flags are then set explicitly below. */
    memset(&sa, 0, sizeof(sa));
    sa.sa_handler = on_sigint;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = 0; /* no SA_RESTART: let blocking calls return so the flag is seen promptly */
    (void) sigaction(SIGINT, &sa, NULL);
}

int osh_signals_should_stop(void *user) {
    (void) user; /* the flag is a process singleton; no per-call context needed */
    return g_stop ? 1 : 0;
}

void osh_signals_reset_stop(void) {
    g_stop = 0;
}

#endif
