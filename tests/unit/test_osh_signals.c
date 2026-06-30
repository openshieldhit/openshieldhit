#include "apps/osh/osh_signals.h"
#include "test_assert.h"

#if !defined(_WIN32)
#include <signal.h> /* raise, SIGINT */
#endif

/* The stop request starts (and resets) clear; the adapter reports it. */
static void test_stop_starts_clear(void) {
    osh_signals_reset_stop();
    ASSERT_TRUE(osh_signals_should_stop(NULL) == 0);
}

/*
 * Installing the handler must be safe and must not raise the flag.  On POSIX a
 * delivered SIGINT must then set it (exercising on_sigint), and a reset clears
 * it again.  On Windows the console Ctrl-C handler runs on its own thread and
 * cannot be triggered portably from a unit test, so there we only assert the
 * install/poll/reset contract, which is identical across platforms.
 */
static void test_install_and_signal(void) {
    osh_signals_reset_stop();
    osh_signals_install_stop();
    ASSERT_TRUE(osh_signals_should_stop(NULL) == 0); /* installing must not raise the flag */

#if !defined(_WIN32)
    raise(SIGINT);                                   /* handler is installed, so this does not terminate */
    ASSERT_TRUE(osh_signals_should_stop(NULL) == 1); /* the handler raised the flag */
    osh_signals_reset_stop();
    ASSERT_TRUE(osh_signals_should_stop(NULL) == 0); /* reset clears it again */
#endif
}

int main(void) {
    test_stop_starts_clear();
    test_install_and_signal();
    return 0;
}
