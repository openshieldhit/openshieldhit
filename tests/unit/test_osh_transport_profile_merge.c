/*
 * Unit tests for osh_transport_profile_merge() — the reduce that folds per-worker
 * transport profiles into a master after a parallel run.  Pure arithmetic, so
 * these run identically on every OS.
 *
 * Coverage:
 *   test_merge_into_zero  — merging one worker into a zeroed dst reproduces it
 *                           exactly (the single-worker / serial case)
 *   test_merge_sums       — counters and per-phase timers add across workers
 *   test_merge_total_max  — total_s combines by maximum, not addition
 *   test_merge_commutative— merge order does not change the result for these
 *                           integer-valued test vectors (see note below)
 *   test_merge_null       — NULL operands are no-ops, not crashes
 *
 * Note on bit-exact assertions: production phase timers are arbitrary doubles
 * read from the monotonic clock, and double addition is not associative, so the
 * real merge is only invariant up to summation order, not byte-for-byte.  These
 * tests deliberately use small integer-valued seconds, which fit exactly in a
 * double and sum without rounding, so `==` comparisons and the order-independence
 * checks below are exact here — a property of the chosen vectors, not a claim
 * about measured timings.
 */

#include <stdio.h>
#include <string.h>

#include "test_assert.h"
#include "transport/osh_transport.h"

/* Bit-exact comparison: these test vectors use integer-valued seconds, which sum
 * in a double without rounding (and max is always exact), so equal results must
 * reproduce exactly here.  This is not true of production timers — see the file
 * header note on summation order. */
static void assert_profile_equal(struct osh_transport_profile const *a, struct osh_transport_profile const *b) {
    ASSERT_TRUE(a->fill_s == b->fill_s);
    ASSERT_TRUE(a->zone_ref_s == b->zone_ref_s);
    ASSERT_TRUE(a->distance_s == b->distance_s);
    ASSERT_TRUE(a->step_s == b->step_s);
    ASSERT_TRUE(a->compact_s == b->compact_s);
    ASSERT_TRUE(a->total_s == b->total_s);
    ASSERT_TRUE(a->steps == b->steps);
    ASSERT_TRUE(a->iterations == b->iterations);
    ASSERT_TRUE(a->nuclear_events == b->nuclear_events);
    ASSERT_TRUE(a->secondaries == b->secondaries);
}

/* ---- Single-worker case: merge into zero reproduces the source exactly --- */

static void test_merge_into_zero(void) {
    struct osh_transport_profile dst;
    struct osh_transport_profile src = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7u, 8u, 9u, 10u};

    memset(&dst, 0, sizeof(dst));
    osh_transport_profile_merge(&dst, &src);
    /* sum-with-zero and max-with-zero both yield src — the serial path is
     * bit-identical whether or not the merge is used. */
    assert_profile_equal(&dst, &src);
}

/* ---- Counters and per-phase timers add across workers -------------------- */

static void test_merge_sums(void) {
    struct osh_transport_profile dst = {1.0, 2.0, 3.0, 4.0, 5.0, 100.0, 7u, 8u, 9u, 10u};
    struct osh_transport_profile src = {10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70u, 80u, 90u, 100u};

    osh_transport_profile_merge(&dst, &src);

    ASSERT_TRUE(dst.fill_s == 11.0);
    ASSERT_TRUE(dst.zone_ref_s == 22.0);
    ASSERT_TRUE(dst.distance_s == 33.0);
    ASSERT_TRUE(dst.step_s == 44.0);
    ASSERT_TRUE(dst.compact_s == 55.0);
    ASSERT_TRUE(dst.steps == 77u);
    ASSERT_TRUE(dst.iterations == 88u);
    ASSERT_TRUE(dst.nuclear_events == 99u);
    ASSERT_TRUE(dst.secondaries == 110u);
    /* total_s is the max (100 vs 60), not the sum. */
    ASSERT_TRUE(dst.total_s == 100.0);
}

/* ---- total_s combines by maximum in both directions ---------------------- */

static void test_merge_total_max(void) {
    struct osh_transport_profile a = {0};
    struct osh_transport_profile b = {0};

    /* dst smaller than src: dst takes src's larger span. */
    a.total_s = 3.0;
    b.total_s = 9.0;
    osh_transport_profile_merge(&a, &b);
    ASSERT_TRUE(a.total_s == 9.0);

    /* dst larger than src: dst keeps its own larger span. */
    a.total_s = 12.0;
    b.total_s = 4.0;
    osh_transport_profile_merge(&a, &b);
    ASSERT_TRUE(a.total_s == 12.0);
}

/* ---- Merge order does not change the result ------------------------------ */

static void test_merge_commutative(void) {
    struct osh_transport_profile base = {1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 1u, 1u, 1u, 1u};
    struct osh_transport_profile b = {2.0, 3.0, 4.0, 5.0, 6.0, 9.0, 2u, 3u, 4u, 5u};
    struct osh_transport_profile c = {20.0, 30.0, 40.0, 50.0, 60.0, 7.0, 20u, 30u, 40u, 50u};
    struct osh_transport_profile dst1 = base;
    struct osh_transport_profile dst2 = base;

    /* dst1: +b then +c ; dst2: +c then +b.  For these integer-valued vectors the
     * sums incur no rounding and max is order-free, so the two must be bit-equal.
     * Arbitrary measured timings would only match up to summation order, not
     * byte-for-byte (see the file header note). */
    osh_transport_profile_merge(&dst1, &b);
    osh_transport_profile_merge(&dst1, &c);
    osh_transport_profile_merge(&dst2, &c);
    osh_transport_profile_merge(&dst2, &b);

    assert_profile_equal(&dst1, &dst2);
    /* total_s is the max of all three spans (9 vs 7 vs 2). */
    ASSERT_TRUE(dst1.total_s == 9.0);
}

/* ---- NULL operands are no-ops, not crashes ------------------------------- */

static void test_merge_null(void) {
    struct osh_transport_profile a = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7u, 8u, 9u, 10u};
    struct osh_transport_profile snapshot = a;

    osh_transport_profile_merge(NULL, &a); /* must not crash */
    osh_transport_profile_merge(&a, NULL); /* must not modify a */
    osh_transport_profile_merge(NULL, NULL);
    assert_profile_equal(&a, &snapshot);
}

int main(void) {
    test_merge_into_zero();
    test_merge_sums();
    test_merge_total_max();
    test_merge_commutative();
    test_merge_null();
    printf("All osh_transport_profile_merge tests passed.\n");
    return 0;
}
