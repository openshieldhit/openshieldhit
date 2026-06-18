#include "apps/osh/osh_membudget.h"

#include <ctype.h>
#include <stdlib.h>
#include <string.h>

/* Smallest absolute reserve we always leave below the detected limit. */
#define OSH_MEMBUDGET_RESERVE_MIN ((uint64_t) 256u << 20) /* 256 MiB */
/* Fraction of the limit used as the budget before applying the reserve floor. */
#define OSH_MEMBUDGET_FRACTION 0.80
/* Fraction of total RAM that also bounds the reserve from below. */
#define OSH_MEMBUDGET_RESERVE_FRACTION 0.05

/* Return the byte multiplier for a (case-insensitive) unit suffix, or 0 if the
 * suffix is not a recognised unit.  "", "B" → 1; "K"/"KB"/"KiB" → 1024; etc.
 * Units are binary (1024-based): "KB" and "KiB" are treated identically. */
static uint64_t unit_multiplier(char const *suffix) {
    char u[4];
    size_t n = 0u;

    while (*suffix == ' ' || *suffix == '\t') {
        ++suffix;
    }
    /* Copy the unit token, stopping at whitespace or end. */
    while (suffix[n] != '\0' && suffix[n] != ' ' && suffix[n] != '\t') {
        if (n >= sizeof(u) - 1u) {
            return 0u; /* longer than any valid unit */
        }
        u[n] = (char) toupper((unsigned char) suffix[n]);
        ++n;
    }
    u[n] = '\0';
    /* Anything after the token must be trailing whitespace only. */
    {
        char const *tail = suffix + n;
        while (*tail == ' ' || *tail == '\t') {
            ++tail;
        }
        if (*tail != '\0') {
            return 0u;
        }
    }

    if (u[0] == '\0' || strcmp(u, "B") == 0) {
        return 1u;
    }
    if (strcmp(u, "K") == 0 || strcmp(u, "KB") == 0 || strcmp(u, "KIB") == 0) {
        return (uint64_t) 1u << 10;
    }
    if (strcmp(u, "M") == 0 || strcmp(u, "MB") == 0 || strcmp(u, "MIB") == 0) {
        return (uint64_t) 1u << 20;
    }
    if (strcmp(u, "G") == 0 || strcmp(u, "GB") == 0 || strcmp(u, "GIB") == 0) {
        return (uint64_t) 1u << 30;
    }
    if (strcmp(u, "T") == 0 || strcmp(u, "TB") == 0 || strcmp(u, "TIB") == 0) {
        return (uint64_t) 1u << 40;
    }
    if (strcmp(u, "P") == 0 || strcmp(u, "PB") == 0 || strcmp(u, "PIB") == 0) {
        return (uint64_t) 1u << 50;
    }
    return 0u;
}

int osh_membudget_parse(char const *text, uint64_t base, uint64_t *out_bytes) {
    char const *p;
    char *end = NULL;
    double value;

    if (!text || !out_bytes) {
        return 0;
    }

    p = text;
    while (*p == ' ' || *p == '\t') {
        ++p;
    }
    if (*p == '\0' || *p == '-' || *p == '+') {
        return 0; /* empty or signed: reject (budgets are non-negative magnitudes) */
    }

    value = strtod(p, &end);
    if (end == p || value < 0.0) {
        return 0; /* no number parsed, or negative */
    }

    /* Skip spaces between the number and its unit/percent marker. */
    while (*end == ' ' || *end == '\t') {
        ++end;
    }

    if (*end == '%') {
        char const *after = end + 1;
        while (*after == ' ' || *after == '\t') {
            ++after;
        }
        if (*after != '\0') {
            return 0; /* trailing junk after % */
        }
        *out_bytes = (uint64_t) ((value / 100.0) * (double) base);
        return 1;
    }

    {
        uint64_t const mult = unit_multiplier(end);
        if (mult == 0u) {
            return 0; /* unrecognised unit */
        }
        *out_bytes = (uint64_t) (value * (double) mult);
        return 1;
    }
}

uint64_t osh_membudget_default(struct osh_sysinfo const *info) {
    uint64_t limit;
    uint64_t reserve;
    uint64_t headroom_cap;
    uint64_t frac_cap;

    if (!info) {
        return 0u;
    }

    /* Prefer available RAM; fall back to total. 0 ⇒ unknown ⇒ no budget. */
    limit = (info->ram_available_bytes > 0u) ? info->ram_available_bytes : info->ram_total_bytes;
    if (limit == 0u) {
        return 0u;
    }

    reserve = OSH_MEMBUDGET_RESERVE_MIN;
    if (info->ram_total_bytes > 0u) {
        uint64_t const frac_reserve = (uint64_t) (OSH_MEMBUDGET_RESERVE_FRACTION * (double) info->ram_total_bytes);
        if (frac_reserve > reserve) {
            reserve = frac_reserve;
        }
    }

    frac_cap = (uint64_t) (OSH_MEMBUDGET_FRACTION * (double) limit);
    headroom_cap = (limit > reserve) ? (limit - reserve) : 0u;

    return (frac_cap < headroom_cap) ? frac_cap : headroom_cap;
}
