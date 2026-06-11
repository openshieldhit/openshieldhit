#ifndef OSH_TIME_H
#define OSH_TIME_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Read a monotonic wall-clock timestamp in seconds.
 *
 * @details
 * Uses CLOCK_MONOTONIC where available, falling back to timespec_get() and
 * finally time().  Only differences between two return values are meaningful;
 * the epoch is unspecified.  Used for progress reporting and the optional
 * profiling phase timers.
 *
 * @returns Monotonic timestamp [s].
 */
double osh_monotonic_seconds(void);

#ifdef __cplusplus
}
#endif

#endif /* OSH_TIME_H */
