#ifndef OSH_DURATION_H
#define OSH_DURATION_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Parse a human-friendly duration string into seconds.
 *
 * @details
 * Shared by the CLI (`--max-time`) and the beam.dat `MAXTIME` keyword so both
 * accept the same syntax.  The grammar is a non-negative decimal number with an
 * optional single-letter unit suffix:
 *
 *   - `s` — seconds  (also the default when no suffix is given)
 *   - `m` — minutes
 *   - `h` — hours
 *
 * Examples: `"30"` → 30, `"30s"` → 30, `"30m"` → 1800, `"1.5h"` → 5400.
 * The suffix is case-insensitive.  Leading/trailing surrounding whitespace is
 * not accepted (callers pass a single trimmed token), and any trailing
 * characters beyond the optional unit are rejected.
 *
 * @param[in]  s            Null-terminated duration token.
 * @param[out] seconds_out  Receives the parsed value in seconds on success.
 *
 * @returns 1 on success, 0 if @p s is NULL/empty, negative, malformed, or has
 *          an unknown suffix or trailing junk.
 */
int osh_parse_duration(char const *s, double *seconds_out);

#ifdef __cplusplus
}
#endif

#endif /* OSH_DURATION_H */
