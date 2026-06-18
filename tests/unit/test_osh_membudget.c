/*
 * Unit tests for the memory-budget parser and default policy (osh_membudget).
 * Pure arithmetic, so these run identically on every OS.
 *
 * Coverage:
 *   test_parse_plain_bytes      — bare integer strings (no unit suffix)
 *   test_parse_all_units        — all recognised unit suffixes (K/KB/KiB … P/PB/PiB)
 *   test_parse_spacing          — leading/trailing/internal whitespace
 *   test_parse_percent          — "N%" form, including base-zero rejection
 *   test_parse_rejects_garbage  — malformed inputs that must return 0
 *   test_default_policy         — osh_membudget_default() policy invariants
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "apps/osh/osh_membudget.h"
#include "test_assert.h"

/* Byte-count helpers — spelled out as products so the intent is obvious. */
#define KIB ((uint64_t) 1024)
#define MIB ((uint64_t) 1024 * 1024)
#define GIB ((uint64_t) 1024 * 1024 * 1024)
#define TIB ((uint64_t) 1024 * 1024 * 1024 * 1024)
#define PIB ((uint64_t) 1024 * 1024 * 1024 * 1024 * 1024)

/* ---- Plain byte strings --------------------------------------------------- */

static void test_parse_plain_bytes(void) {
    uint64_t v = 0u;
    /* Exact integer byte counts, no unit suffix. */
    ASSERT_TRUE(osh_membudget_parse("0", 0u, &v) && v == 0u);
    ASSERT_TRUE(osh_membudget_parse("1", 0u, &v) && v == 1u);
    ASSERT_TRUE(osh_membudget_parse("1048576", 0u, &v) && v == 1048576u);
    ASSERT_TRUE(osh_membudget_parse("4294967295", 0u, &v) && v == 4294967295u); /* 2^32-1 */
}

/* ---- Unit suffix variants ------------------------------------------------- */

/*
 * unit_multiplier() is a private static function; test it indirectly through
 * osh_membudget_parse().  Every suffix variant (short, long, SI-looking, IEC)
 * must produce the same binary (1024-based) result.
 */
static void test_parse_all_units(void) {
    uint64_t v = 0u;

    /* Plain bytes — "B" suffix or no suffix at all. */
    ASSERT_TRUE(osh_membudget_parse("1B", 0u, &v) && v == 1u);
    ASSERT_TRUE(osh_membudget_parse("512B", 0u, &v) && v == 512u);

    /* Kibibytes: K, KB, KiB (all lowercase/mixed-case too). */
    ASSERT_TRUE(osh_membudget_parse("1K", 0u, &v) && v == KIB);
    ASSERT_TRUE(osh_membudget_parse("1k", 0u, &v) && v == KIB);
    ASSERT_TRUE(osh_membudget_parse("1KB", 0u, &v) && v == KIB);
    ASSERT_TRUE(osh_membudget_parse("1kb", 0u, &v) && v == KIB);
    ASSERT_TRUE(osh_membudget_parse("1KiB", 0u, &v) && v == KIB);
    ASSERT_TRUE(osh_membudget_parse("1kib", 0u, &v) && v == KIB);

    /* Mebibytes: M, MB, MiB. */
    ASSERT_TRUE(osh_membudget_parse("1M", 0u, &v) && v == MIB);
    ASSERT_TRUE(osh_membudget_parse("2M", 0u, &v) && v == 2u * MIB);
    ASSERT_TRUE(osh_membudget_parse("1MB", 0u, &v) && v == MIB);
    ASSERT_TRUE(osh_membudget_parse("1MiB", 0u, &v) && v == MIB);
    ASSERT_TRUE(osh_membudget_parse("256MB", 0u, &v) && v == 256u * MIB);

    /* Gibibytes: G, GB, GiB. */
    ASSERT_TRUE(osh_membudget_parse("1G", 0u, &v) && v == GIB);
    ASSERT_TRUE(osh_membudget_parse("8GB", 0u, &v) && v == 8u * GIB);
    ASSERT_TRUE(osh_membudget_parse("1GiB", 0u, &v) && v == GIB);
    ASSERT_TRUE(osh_membudget_parse("1.5GiB", 0u, &v) && v == (uint64_t) (1.5 * (double) GIB));

    /* Tebibytes: T, TB, TiB. */
    ASSERT_TRUE(osh_membudget_parse("1T", 0u, &v) && v == TIB);
    ASSERT_TRUE(osh_membudget_parse("1TB", 0u, &v) && v == TIB);
    ASSERT_TRUE(osh_membudget_parse("1TiB", 0u, &v) && v == TIB);
    ASSERT_TRUE(osh_membudget_parse("2TB", 0u, &v) && v == 2u * TIB);

    /* Pebibytes: P, PB, PiB. */
    ASSERT_TRUE(osh_membudget_parse("1P", 0u, &v) && v == PIB);
    ASSERT_TRUE(osh_membudget_parse("1PB", 0u, &v) && v == PIB);
    ASSERT_TRUE(osh_membudget_parse("1PiB", 0u, &v) && v == PIB);
}

/* ---- Whitespace tolerance ------------------------------------------------- */

static void test_parse_spacing(void) {
    uint64_t v = 0u;
    /* Leading/trailing spaces around the whole string. */
    ASSERT_TRUE(osh_membudget_parse("  512MB", 0u, &v) && v == 512u * MIB);
    ASSERT_TRUE(osh_membudget_parse("512MB  ", 0u, &v) && v == 512u * MIB);
    ASSERT_TRUE(osh_membudget_parse("  512MB  ", 0u, &v) && v == 512u * MIB);
    /* Space between number and unit: "4 GiB" → same as "4GiB". */
    ASSERT_TRUE(osh_membudget_parse("4 GiB", 0u, &v) && v == 4u * GIB);
    ASSERT_TRUE(osh_membudget_parse("  4 GiB ", 0u, &v) && v == 4u * GIB);
    /* Plain byte string with spaces. */
    ASSERT_TRUE(osh_membudget_parse("  1048576  ", 0u, &v) && v == 1048576u);
}

/* ---- Percentage form ------------------------------------------------------ */

static void test_parse_percent(void) {
    uint64_t v = 0u;

    /* Standard uses. */
    ASSERT_TRUE(osh_membudget_parse("50%", 1000u, &v) && v == 500u);
    ASSERT_TRUE(osh_membudget_parse("100%", 4096u, &v) && v == 4096u);
    ASSERT_TRUE(osh_membudget_parse("80%", 10u * GIB, &v) && v == (uint64_t) (0.80 * (double) (10u * GIB)));

    /* Fractional percentage. */
    ASSERT_TRUE(osh_membudget_parse("62.5%", 1024u, &v) && v == 640u);

    /* Whitespace around the whole string. */
    ASSERT_TRUE(osh_membudget_parse("  50%  ", 1000u, &v) && v == 500u);

    /* base == 0: percentage has no meaningful reference → must fail.
     * Without this guard "80%" with unknown RAM would silently produce 0,
     * which would disable the budget rather than reporting bad input. */
    ASSERT_TRUE(!osh_membudget_parse("80%", 0u, &v));
    ASSERT_TRUE(!osh_membudget_parse("0%", 0u, &v));
}

/* ---- Rejection of malformed inputs --------------------------------------- */

static void test_parse_rejects_garbage(void) {
    uint64_t v = 0u;

    /* NULL or empty. */
    ASSERT_TRUE(!osh_membudget_parse(NULL, 0u, &v));
    ASSERT_TRUE(!osh_membudget_parse("", 0u, &v));
    ASSERT_TRUE(!osh_membudget_parse("   ", 0u, &v));

    /* Signed values (budgets are non-negative). */
    ASSERT_TRUE(!osh_membudget_parse("-5", 0u, &v));
    ASSERT_TRUE(!osh_membudget_parse("+5", 0u, &v));
    ASSERT_TRUE(!osh_membudget_parse("-1GB", 0u, &v));

    /* Non-numeric start. */
    ASSERT_TRUE(!osh_membudget_parse("abc", 0u, &v));
    ASSERT_TRUE(!osh_membudget_parse("GB", 0u, &v));

    /* Unrecognised unit suffix. */
    ASSERT_TRUE(!osh_membudget_parse("10XB", 0u, &v));
    ASSERT_TRUE(!osh_membudget_parse("10YiB", 0u, &v));
    ASSERT_TRUE(!osh_membudget_parse("10bytes", 0u, &v));

    /* Trailing junk after a valid unit. */
    ASSERT_TRUE(!osh_membudget_parse("1GB junk", 0u, &v));
    ASSERT_TRUE(!osh_membudget_parse("1GBextra", 0u, &v));

    /* Trailing junk after a percentage. */
    ASSERT_TRUE(!osh_membudget_parse("80% extra", 1000u, &v));
    ASSERT_TRUE(!osh_membudget_parse("80%MB", 1000u, &v));
}

/* ---- Default policy invariants ------------------------------------------- */

static void test_default_policy(void) {
    struct osh_sysinfo info;
    uint64_t budget;

    /* NULL guard. */
    ASSERT_TRUE(osh_membudget_default(NULL) == 0u);

    /* Both total and available known: budget should be between 50 % and 100 %
     * of available (exact value depends on the reserve calculation). */
    info.logical_cores = 8u;
    info.ram_total_bytes = 16u * GIB;
    info.ram_available_bytes = 10u * GIB;
    info.gpu_count = 0u;
    budget = osh_membudget_default(&info);
    ASSERT_TRUE(budget > 0u);
    ASSERT_TRUE(budget < info.ram_available_bytes);       /* always leaves headroom */
    ASSERT_TRUE(budget >= info.ram_available_bytes / 2u); /* at least 50 % */

    /* Available unknown (0): falls back to total RAM. */
    info.ram_available_bytes = 0u;
    budget = osh_membudget_default(&info);
    ASSERT_TRUE(budget > 0u && budget < info.ram_total_bytes);

    /* Both unknown: return 0 so the caller skips enforcement. */
    info.ram_total_bytes = 0u;
    budget = osh_membudget_default(&info);
    ASSERT_TRUE(budget == 0u);

    /* Very small RAM (256 MiB total): reserve floor consumes all of it.
     * Budget must be 0 rather than a nonsensical huge value. */
    info.ram_total_bytes = 256u * MIB;
    info.ram_available_bytes = 200u * MIB;
    budget = osh_membudget_default(&info);
    ASSERT_TRUE(budget == 0u || budget < info.ram_available_bytes);
}

int main(void) {
    test_parse_plain_bytes();
    test_parse_all_units();
    test_parse_spacing();
    test_parse_percent();
    test_parse_rejects_garbage();
    test_default_policy();
    printf("All osh_membudget tests passed.\n");
    return 0;
}
