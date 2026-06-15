#include "common/osh_time.h"

#include <time.h>

/* Fallback chain: CLOCK_MONOTONIC (POSIX) -> timespec_get(TIME_UTC) (C11) ->
 * time().  The latter two are wall-clock sources and can move backwards on
 * NTP adjustments, producing negative phase durations in the profiling output.
 * This is acceptable for a best-effort perf timer; correctness is not
 * affected.  A platform-specific monotonic source (e.g. QueryPerformanceCounter
 * on Windows) is not worth the added complexity here. */
double osh_monotonic_seconds(void) {
#if defined(CLOCK_MONOTONIC)
    {
        struct timespec ts;
        if (clock_gettime(CLOCK_MONOTONIC, &ts) == 0) {
            return (double) ts.tv_sec + 1.0e-9 * (double) ts.tv_nsec;
        }
    }
#endif

#if defined(TIME_UTC)
    {
        struct timespec ts;
        if (timespec_get(&ts, TIME_UTC) == TIME_UTC) {
            return (double) ts.tv_sec + 1.0e-9 * (double) ts.tv_nsec;
        }
    }
#endif

    return (double) time(NULL);
}
