#include "apps/osh/osh_signals.h"

#include <stddef.h> /* NULL */
#include <string.h> /* memset */

/*
 * Graceful-stop flag.  The only thing the OS handler touches, and it only ever
 * writes 1 — that keeps the handler async-signal-safe (no I/O, no allocation,
 * no library calls).  volatile sig_atomic_t is the one type the C standard
 * guarantees can be written in a signal handler and read elsewhere without a
 * data race.
 */
static sig_atomic_t volatile g_stop = 0;

#if defined(_WIN32)

#include <windows.h>

static BOOL WINAPI on_console_ctrl(DWORD ctrl_type) {
    switch (ctrl_type) {
    case CTRL_C_EVENT:
    case CTRL_BREAK_EVENT:
    case CTRL_CLOSE_EVENT:
    case CTRL_LOGOFF_EVENT:
    case CTRL_SHUTDOWN_EVENT:
        g_stop = 1;
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

#else /* POSIX */

#include <signal.h>

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

#endif

sig_atomic_t volatile *osh_signals_stop_flag(void) {
    return &g_stop;
}
