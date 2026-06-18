/*
 * Unit tests for the memory-budget parser and default policy (osh_membudget).
 * Pure arithmetic, so these run identically on every OS.
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "apps/osh/osh_membudget.h"

#define ASSERT_TRUE(cond)                                                                                              \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            fprintf(stderr, "ASSERT FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__);                                 \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

#define KIB ((uint64_t) 1u << 10)
#define MIB ((uint64_t) 1u << 20)
#define GIB ((uint64_t) 1u << 30)

static void test_parse_plain_bytes(void) {
    uint64_t v = 0u;
    ASSERT_TRUE(osh_membudget_parse("1048576", 0u, &v) && v == 1048576u);
    ASSERT_TRUE(osh_membudget_parse("0", 0u, &v) && v == 0u);
}

static void test_parse_units(void) {
    uint64_t v = 0u;
    ASSERT_TRUE(osh_membudget_parse("1KB", 0u, &v) && v == KIB);
    ASSERT_TRUE(osh_membudget_parse("1KiB", 0u, &v) && v == KIB);
    ASSERT_TRUE(osh_membudget_parse("1k", 0u, &v) && v == KIB);
    ASSERT_TRUE(osh_membudget_parse("2M", 0u, &v) && v == 2u * MIB);
    ASSERT_TRUE(osh_membudget_parse("8GB", 0u, &v) && v == 8u * GIB);
    ASSERT_TRUE(osh_membudget_parse("1.5GiB", 0u, &v) && v == (uint64_t) (1.5 * (double) GIB));
    ASSERT_TRUE(osh_membudget_parse("  4 GiB ", 0u, &v) && v == 4u * GIB);
}

static void test_parse_percent(void) {
    uint64_t v = 0u;
    ASSERT_TRUE(osh_membudget_parse("50%", 1000u, &v) && v == 500u);
    ASSERT_TRUE(osh_membudget_parse("80%", 10u * GIB, &v) && v == (uint64_t) (0.80 * (double) (10u * GIB)));
    ASSERT_TRUE(osh_membudget_parse("100%", 4096u, &v) && v == 4096u);
}

static void test_parse_rejects_garbage(void) {
    uint64_t v = 0u;
    ASSERT_TRUE(!osh_membudget_parse("", 0u, &v));
    ASSERT_TRUE(!osh_membudget_parse("abc", 0u, &v));
    ASSERT_TRUE(!osh_membudget_parse("-5", 0u, &v));
    ASSERT_TRUE(!osh_membudget_parse("10XB", 0u, &v));
    ASSERT_TRUE(!osh_membudget_parse("1GB junk", 0u, &v));
    ASSERT_TRUE(!osh_membudget_parse(NULL, 0u, &v));
}

static void test_default_policy(void) {
    struct osh_sysinfo info;
    uint64_t budget;

    /* Known total + available: budget ≈ 80% of available, below available. */
    info.logical_cores = 8u;
    info.ram_total_bytes = 16u * GIB;
    info.ram_available_bytes = 10u * GIB;
    info.gpu_count = 0u;
    budget = osh_membudget_default(&info);
    ASSERT_TRUE(budget > 0u);
    ASSERT_TRUE(budget < info.ram_available_bytes); /* leaves headroom */
    ASSERT_TRUE(budget >= info.ram_available_bytes / 2u);

    /* Available unknown: falls back to total. */
    info.ram_available_bytes = 0u;
    budget = osh_membudget_default(&info);
    ASSERT_TRUE(budget > 0u && budget < info.ram_total_bytes);

    /* Nothing known: 0 (caller skips enforcement). */
    info.ram_total_bytes = 0u;
    budget = osh_membudget_default(&info);
    ASSERT_TRUE(budget == 0u);

    /* NULL is handled. */
    ASSERT_TRUE(osh_membudget_default(NULL) == 0u);
}

int main(void) {
    test_parse_plain_bytes();
    test_parse_units();
    test_parse_percent();
    test_parse_rejects_garbage();
    test_default_policy();
    return 0;
}
