#include "common/osh_duration.h"

#include <errno.h>
#include <stdlib.h>

int osh_parse_duration(char const *s, double *seconds_out) {
    double value;
    double scale;
    char *end = NULL;

    if (!s || !*s || !seconds_out) {
        return 0;
    }
    /* Reject leading whitespace / sign so the grammar matches the CLI's other
     * numeric options (strtod would otherwise accept " 30", "+30", "-30"). */
    if (s[0] != '.' && (s[0] < '0' || s[0] > '9')) {
        return 0;
    }

    errno = 0;
    value = strtod(s, &end);
    if (errno != 0 || end == s || value < 0.0) {
        return 0;
    }

    /* Optional single-letter unit suffix; default (no suffix) is seconds. */
    scale = 1.0;
    if (*end != '\0') {
        switch (*end) {
        case 's':
        case 'S':
            scale = 1.0;
            break;
        case 'm':
        case 'M':
            scale = 60.0;
            break;
        case 'h':
        case 'H':
            scale = 3600.0;
            break;
        default:
            return 0;
        }
        ++end;
        if (*end != '\0') {
            return 0; /* trailing junk after the unit */
        }
    }

    *seconds_out = value * scale;
    return 1;
}
