#include <signal.h>

#include "apps/osh/osh_signals.h"
#include "test_assert.h"

/* The stop-flag accessor always returns a usable, initially-clear flag. */
static void test_stop_flag_accessor(void) {
    sig_atomic_t volatile *flag = osh_signals_stop_flag();
    ASSERT_TRUE(flag != NULL);
}

/*
 * Installing the handler must be safe, and on POSIX a delivered SIGINT must set
 * the stop flag (exercising on_sigint).  On Windows the console Ctrl-C handler
 * cannot be triggered portably from a unit test, so there we only assert that
 * install is a safe no-op and the flag stays readable; the flag contract itself
 * is what callers rely on and it is identical across platforms.
 */
static void test_install_and_signal(void) {
    sig_atomic_t volatile *flag = osh_signals_stop_flag();

    *flag = 0;
    osh_signals_install_stop();
    ASSERT_TRUE(*flag == 0); /* installing must not raise the flag */

#if !defined(_WIN32)
    raise(SIGINT);           /* handler is installed, so this does not terminate */
    ASSERT_TRUE(*flag == 1); /* the handler raised the flag */
    *flag = 0;               /* leave it clear for any later use */
#endif
}

int main(void) {
    test_stop_flag_accessor();
    test_install_and_signal();
    return 0;
}
