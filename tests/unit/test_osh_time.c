#include "common/osh_time.h"
#include "test_assert.h"

static void test_monotonic_seconds_is_monotonic(void) {
    double t0;
    double t1;
    int i;

    t0 = osh_monotonic_seconds();
    ASSERT_TRUE(t0 > 0.0);

    /* Successive reads never go backwards, and the clock visibly advances
     * within a bounded busy-wait (sub-second granularity is required for the
     * transport phase timers to be useful). */
    for (i = 0; i < 1000; ++i) {
        t1 = osh_monotonic_seconds();
        ASSERT_TRUE(t1 >= t0);
        t0 = t1;
    }

    t0 = osh_monotonic_seconds();
    do {
        t1 = osh_monotonic_seconds();
        ASSERT_TRUE(t1 >= t0);
    } while (t1 - t0 < 1.0e-6);
    ASSERT_TRUE(t1 - t0 < 60.0);
}

int main(void) {
    test_monotonic_seconds_is_monotonic();
    return 0;
}
