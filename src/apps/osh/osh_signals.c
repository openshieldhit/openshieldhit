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

/* Windows has no SIGUSR1 equivalent, so on-demand dumps are unavailable here.
 * Keep the symbols so main.c and tests link identically on every platform: the
 * install is a no-op and the poll always reports "no dump requested".  Scheduled
 * --dump-every[-primaries] dumps are unaffected — they never use these. */
void osh_signals_install_dump(void) {
}

int osh_signals_should_dump(void *user) {
    (void) user;
    return 0;
}

void osh_signals_reset_dump(void) {
}

#else /* POSIX */

#include <signal.h>
#include <string.h> /* memset */

/* The only writer is the SIGINT handler, which runs in the context of the
 * interrupted thread (not a separate one), so `volatile sig_atomic_t` — the
 * type the C standard guarantees safe between a signal handler and the main
 * flow — is exactly right, and writing the flag is async-signal-safe. */
static sig_atomic_t volatile g_stop = 0;

/* On-demand dump flag, raised by SIGUSR1 and consumed (read-and-cleared) by the
 * run loop.  Same concurrency model as g_stop: the handler runs in the
 * interrupted thread, so volatile sig_atomic_t is the right, async-signal-safe
 * type.  Best-effort edge trigger, unlike the level-held g_stop: one or more
 * SIGUSR1s pending before a checkpoint coalesce into a single dump (standard
 * POSIX signals coalesce, so N rapid signals need not yield N dumps). */
static sig_atomic_t volatile g_dump = 0;

static void on_sigint(int sig) {
    (void) sig;
    g_stop = 1; /* async-signal-safe: flag only */
}

static void on_sigusr1(int sig) {
    (void) sig;
    g_dump = 1; /* async-signal-safe: flag only; dump happens at the next checkpoint */
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

void osh_signals_install_dump(void) {
    struct sigaction sa;
    /* Zero every field first (see osh_signals_install_stop for why), then set the
     * handler explicitly.  SIGUSR1 is the conventional user-defined signal for a
     * dump-and-continue request. */
    memset(&sa, 0, sizeof(sa));
    sa.sa_handler = on_sigusr1;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = 0; /* no SA_RESTART: let a blocking call return so the flag is seen promptly */
    (void) sigaction(SIGUSR1, &sa, NULL);
}

int osh_signals_should_dump(void *user) {
    (void) user; /* the flag is a process singleton; no per-call context needed */
    /* Best-effort edge trigger: read and clear so pending SIGUSR1 requests produce
     * one dump at this poll.  This is not a lossless one-dump-per-signal guarantee:
     * signals coalesce, and a SIGUSR1 landing between the read and the clear is
     * dropped rather than deferred.  That is acceptable for a dump-and-continue
     * hint — the user simply re-signals, and a scheduled cadence bounds latency. */
    if (g_dump) {
        g_dump = 0;
        return 1;
    }
    return 0;
}

void osh_signals_reset_dump(void) {
    g_dump = 0;
}

#endif
