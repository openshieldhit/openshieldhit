#include "apps/osh/osh_membudget.h"

#include <ctype.h>
#include <stdlib.h>
#include <string.h>

/*
 * Policy constants — all #defined so they appear in user-facing error messages
 * and can be patched at compile time for specialised deployments.
 */

/* Smallest absolute reserve we always keep below the detected memory limit.
 * Written as a product of visible factors so the value is self-documenting. */
#define OSH_MEMBUDGET_RESERVE_MIN ((uint64_t) 256 * 1024 * 1024) /* 256 MiB */

/* Fraction of the memory limit to offer as the preliminary budget (before
 * applying the reserve floor).  0.80 means "use at most 80 % of available". */
#define OSH_MEMBUDGET_FRACTION 0.80

/* Lower bound on the reserve expressed as a fraction of *total* RAM.
 * On machines with little RAM (e.g. 1 GiB), 5 % keeps the reserve larger
 * than 50 MiB even when OSH_MEMBUDGET_RESERVE_MIN would be the binding limit. */
#define OSH_MEMBUDGET_RESERVE_FRACTION 0.05

/*
 * unit_multiplier — map a binary-size suffix to its byte multiplier.
 *
 * Accepts a suffix string that immediately follows the numeric part of a size
 * token.  Leading whitespace is consumed.  The token is matched case-insensitively
 * and may be followed only by whitespace; any extra non-whitespace characters
 * cause a 0 return so the caller can report an unrecognised unit.
 *
 * Recognised suffixes and their multipliers (all 1024-based / binary):
 *
 *   ""  / "B"             →          1  (plain bytes)
 *   "K" / "KB" / "KiB"   →      1 024  (kibibytes)
 *   "M" / "MB" / "MiB"   →  1 048 576  (mebibytes)
 *   "G" / "GB" / "GiB"   →  1 073 741 824            (gibibytes)
 *   "T" / "TB" / "TiB"   →  1 099 511 627 776        (tebibytes)
 *   "P" / "PB" / "PiB"   →  1 125 899 906 842 624    (pebibytes)
 *
 * "KB" and "KiB" are treated as synonyms — both mean 1024 bytes.
 * Returns 0 for any unrecognised or malformed suffix.
 */
static uint64_t unit_multiplier(char const *suffix) {
    char u[4]; /* longest valid token is "KiB\0" = 4 chars */
    size_t n = 0u;

    /* Skip optional whitespace between the number and the unit token. */
    while (*suffix == ' ' || *suffix == '\t') {
        ++suffix;
    }

    /* Copy the unit token into a local buffer, upper-casing as we go.
     * Stop at whitespace or end-of-string; reject tokens longer than 3 chars. */
    while (suffix[n] != '\0' && suffix[n] != ' ' && suffix[n] != '\t') {
        if (n >= sizeof(u) - 1u) {
            return 0u; /* longer than any valid unit ("PiB" is the longest = 3 chars) */
        }
        u[n] = (char) toupper((unsigned char) suffix[n]);
        ++n;
    }
    u[n] = '\0';

    /* Anything after the token must be trailing whitespace only.
     * e.g. "1GB " is OK but "1GB extra" is not. */
    {
        char const *tail = suffix + n;
        while (*tail == ' ' || *tail == '\t') {
            ++tail;
        }
        if (*tail != '\0') {
            return 0u;
        }
    }

    /* Match the normalised token. Empty string or bare "B" means plain bytes. */
    if (u[0] == '\0' || strcmp(u, "B") == 0) {
        return 1u;
    }
    if (strcmp(u, "K") == 0 || strcmp(u, "KB") == 0 || strcmp(u, "KIB") == 0) {
        return (uint64_t) 1u << 10; /* 1 024 */
    }
    if (strcmp(u, "M") == 0 || strcmp(u, "MB") == 0 || strcmp(u, "MIB") == 0) {
        return (uint64_t) 1u << 20; /* 1 048 576 */
    }
    if (strcmp(u, "G") == 0 || strcmp(u, "GB") == 0 || strcmp(u, "GIB") == 0) {
        return (uint64_t) 1u << 30; /* 1 073 741 824 */
    }
    if (strcmp(u, "T") == 0 || strcmp(u, "TB") == 0 || strcmp(u, "TIB") == 0) {
        return (uint64_t) 1u << 40; /* 1 099 511 627 776 */
    }
    if (strcmp(u, "P") == 0 || strcmp(u, "PB") == 0 || strcmp(u, "PIB") == 0) {
        return (uint64_t) 1u << 50; /* 1 125 899 906 842 624 */
    }
    return 0u; /* unrecognised */
}

int osh_membudget_parse(char const *text, uint64_t base, uint64_t *out_bytes) {
    char const *p;
    char *end = NULL;
    double value;

    if (!text || !out_bytes) {
        return 0;
    }

    /* Skip leading whitespace; reject immediately if the string is empty or
     * starts with a sign character (budgets are non-negative magnitudes). */
    p = text;
    while (*p == ' ' || *p == '\t') {
        ++p;
    }
    if (*p == '\0' || *p == '-' || *p == '+') {
        return 0;
    }

    /* Parse the leading numeric part.  strtod handles integer and fractional
     * forms ("512", "1.5", "62.5") and stops at the first non-numeric char. */
    value = strtod(p, &end);
    if (end == p || value < 0.0) {
        return 0; /* no number parsed, or negative */
    }

    /* Skip optional whitespace between the number and its unit/percent marker
     * so "4 GiB" and "4GiB" are both accepted. */
    while (*end == ' ' || *end == '\t') {
        ++end;
    }

    /* --- Percentage form: "80%" means 80 % of base ---
     * base == 0 is rejected because the result would always be 0, which would
     * silently disable the budget instead of refusing an unusable input. */
    if (*end == '%') {
        char const *after = end + 1;
        double bytes_d;

        /* Consume any trailing whitespace after '%'; any other character is junk. */
        while (*after == ' ' || *after == '\t') {
            ++after;
        }
        if (*after != '\0' || base == 0u) {
            return 0; /* trailing junk after % or unknown base */
        }

        bytes_d = (value / 100.0) * (double) base;
        /* Guard against NaN, negative results, and values too large for uint64_t. */
        if (!(bytes_d >= 0.0) || bytes_d > (double) UINT64_MAX) {
            return 0;
        }
        *out_bytes = (uint64_t) bytes_d;
        return 1;
    }

    /* --- Unit suffix form: "512MB", "4GiB", "2k", "1048576B", "" (plain bytes) ---
     * unit_multiplier() handles the suffix and rejects trailing junk. */
    {
        uint64_t const mult = unit_multiplier(end);
        double bytes_d;

        if (mult == 0u) {
            return 0; /* unrecognised unit or trailing garbage */
        }

        bytes_d = value * (double) mult;
        /* Guard against NaN, negative results, and values too large for uint64_t. */
        if (!(bytes_d >= 0.0) || bytes_d > (double) UINT64_MAX) {
            return 0;
        }
        *out_bytes = (uint64_t) bytes_d;
        return 1;
    }
}

uint64_t osh_membudget_default(struct osh_sysinfo const *info) {
    uint64_t limit;        /* The memory figure we base the budget on. */
    uint64_t reserve;      /* Bytes kept back as a safety margin. */
    uint64_t frac_cap;     /* Budget from the OSH_MEMBUDGET_FRACTION percentage. */
    uint64_t headroom_cap; /* Budget from the reserve-floor approach. */

    if (!info) {
        return 0u;
    }

    /* Choose the memory baseline: prefer currently *available* RAM because it
     * already accounts for other running processes.  Fall back to total RAM on
     * platforms that cannot report available RAM (e.g. some WSL1 versions or
     * kernels without /proc/meminfo MemAvailable).  0 means unknown → no budget. */
    if (info->ram_available_bytes > 0u) {
        limit = info->ram_available_bytes;
    } else {
        limit = info->ram_total_bytes;
    }
    if (limit == 0u) {
        return 0u;
    }

    /* Compute the reserve floor: the larger of the hard minimum (256 MiB) and
     * 5 % of *total* RAM.  Using total (not available) keeps the reserve stable
     * and independent of transient allocations by other processes.
     * On a 32 GiB machine: 5 % = 1.6 GiB > 256 MiB → reserve = 1.6 GiB.
     * On a  2 GiB machine: 5 % = 100 MiB < 256 MiB → reserve = 256 MiB. */
    reserve = OSH_MEMBUDGET_RESERVE_MIN;
    if (info->ram_total_bytes > 0u) {
        uint64_t const pct_reserve = (uint64_t) (OSH_MEMBUDGET_RESERVE_FRACTION * (double) info->ram_total_bytes);
        if (pct_reserve > reserve) {
            reserve = pct_reserve;
        }
    }

    /* Two independent caps; we return whichever is *smaller* so both constraints
     * are always satisfied simultaneously:
     *
     *   frac_cap     = 80 % of limit
     *                  → never use more than 80 % of (available/total) RAM.
     *
     *   headroom_cap = limit - reserve
     *                  → always leave at least `reserve` bytes untouched for the
     *                    OS, other processes, and OpenShieldHIT's own (small)
     *                    non-scoring allocations.
     *
     * If limit ≤ reserve the machine is already tight; return 0 to skip
     * budget enforcement rather than returning a nonsensical negative-ish value. */
    frac_cap = (uint64_t) (OSH_MEMBUDGET_FRACTION * (double) limit);
    headroom_cap = (limit > reserve) ? (limit - reserve) : 0u;

    return (frac_cap < headroom_cap) ? frac_cap : headroom_cap;
}
