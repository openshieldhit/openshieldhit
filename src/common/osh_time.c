#include "common/osh_time.h"

#include <time.h>

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
