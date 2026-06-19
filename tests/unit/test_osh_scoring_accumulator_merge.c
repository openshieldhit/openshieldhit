/*
 * Unit tests for osh_scoring_accumulator_merge() — the element-wise reduce used
 * to fold per-worker scoring accumulators into a master set.  Pure arithmetic,
 * so these run identically on every OS.
 *
 * Coverage:
 *   test_merge_identity        — adding a zero accumulator leaves dst unchanged
 *   test_merge_correctness     — three partial sums fold to the exact total
 *   test_merge_commutative     — merge order does not change the result
 *   test_merge_len_mismatch    — differing lengths are rejected (OSH_EINVAL)
 *   test_merge_null            — NULL operands are rejected without crashing
 */

#include <stddef.h>

#include "openshieldhit/status.h"
#include "scoring/runtime/osh_scoring_accumulator.h"
#include "test_assert.h"

/* Bit-exact comparison: the merge only adds, so integer-valued doubles that fit
 * exactly in a double must reproduce exactly regardless of order. */
static void assert_data_equals(struct osh_scoring_accumulator const *acc, double const *expect, size_t n) {
    size_t i;
    ASSERT_TRUE(acc->len == n);
    for (i = 0; i < n; ++i) {
        ASSERT_TRUE(acc->data[i] == expect[i]);
    }
}

/* ---- Identity: dst += 0 is a no-op --------------------------------------- */

static void test_merge_identity(void) {
    struct osh_scoring_accumulator dst;
    struct osh_scoring_accumulator zero;
    double const expect[4] = {1.0, 2.0, 3.0, 4.0};
    size_t i;

    ASSERT_TRUE(osh_scoring_accumulator_alloc(&dst, 4u, 0) == OSH_OK);
    ASSERT_TRUE(osh_scoring_accumulator_alloc(&zero, 4u, 0) == OSH_OK);
    for (i = 0; i < 4u; ++i) {
        dst.data[i] = expect[i];
    }

    ASSERT_TRUE(osh_scoring_accumulator_merge(&dst, &zero) == OSH_OK);
    assert_data_equals(&dst, expect, 4u);

    osh_scoring_accumulator_free(&dst);
    osh_scoring_accumulator_free(&zero);
}

/* ---- Correctness: three partial sums fold to the exact total ------------- */

static void test_merge_correctness(void) {
    struct osh_scoring_accumulator total;
    struct osh_scoring_accumulator a;
    struct osh_scoring_accumulator b;
    struct osh_scoring_accumulator c;
    double const expect_data[3] = {1.0 + 10.0 + 100.0, 2.0 + 20.0 + 200.0, 3.0 + 30.0 + 300.0};
    double const expect_data2[3] = {4.0 + 40.0 + 400.0, 5.0 + 50.0 + 500.0, 6.0 + 60.0 + 600.0};
    size_t i;

    /* want_data2 = 1 so the secondary weight array is exercised too. */
    ASSERT_TRUE(osh_scoring_accumulator_alloc(&total, 3u, 1) == OSH_OK);
    ASSERT_TRUE(osh_scoring_accumulator_alloc(&a, 3u, 1) == OSH_OK);
    ASSERT_TRUE(osh_scoring_accumulator_alloc(&b, 3u, 1) == OSH_OK);
    ASSERT_TRUE(osh_scoring_accumulator_alloc(&c, 3u, 1) == OSH_OK);

    for (i = 0; i < 3u; ++i) {
        a.data[i] = (double) (i + 1);          /* 1,2,3 */
        b.data[i] = (double) (10 * (i + 1));   /* 10,20,30 */
        c.data[i] = (double) (100 * (i + 1));  /* 100,200,300 */
        a.data2[i] = (double) (4 + i);         /* 4,5,6 */
        b.data2[i] = (double) (40 + 10 * i);   /* 40,50,60 */
        c.data2[i] = (double) (400 + 100 * i); /* 400,500,600 */
    }

    ASSERT_TRUE(osh_scoring_accumulator_merge(&total, &a) == OSH_OK);
    ASSERT_TRUE(osh_scoring_accumulator_merge(&total, &b) == OSH_OK);
    ASSERT_TRUE(osh_scoring_accumulator_merge(&total, &c) == OSH_OK);

    for (i = 0; i < 3u; ++i) {
        ASSERT_TRUE(total.data[i] == expect_data[i]);
        ASSERT_TRUE(total.data2[i] == expect_data2[i]);
    }

    osh_scoring_accumulator_free(&total);
    osh_scoring_accumulator_free(&a);
    osh_scoring_accumulator_free(&b);
    osh_scoring_accumulator_free(&c);
}

/* ---- Commutativity: merge order does not change the result --------------- */

static void test_merge_commutative(void) {
    struct osh_scoring_accumulator dst1;
    struct osh_scoring_accumulator dst2;
    struct osh_scoring_accumulator b;
    struct osh_scoring_accumulator c;
    double const seed[3] = {7.0, 8.0, 9.0};
    size_t i;

    ASSERT_TRUE(osh_scoring_accumulator_alloc(&dst1, 3u, 0) == OSH_OK);
    ASSERT_TRUE(osh_scoring_accumulator_alloc(&dst2, 3u, 0) == OSH_OK);
    ASSERT_TRUE(osh_scoring_accumulator_alloc(&b, 3u, 0) == OSH_OK);
    ASSERT_TRUE(osh_scoring_accumulator_alloc(&c, 3u, 0) == OSH_OK);
    for (i = 0; i < 3u; ++i) {
        dst1.data[i] = seed[i];
        dst2.data[i] = seed[i];
        b.data[i] = (double) (i + 1);
        c.data[i] = (double) (100 * (i + 1));
    }

    /* dst1: +b then +c ; dst2: +c then +b. Integer-valued => must be bit-equal. */
    ASSERT_TRUE(osh_scoring_accumulator_merge(&dst1, &b) == OSH_OK);
    ASSERT_TRUE(osh_scoring_accumulator_merge(&dst1, &c) == OSH_OK);
    ASSERT_TRUE(osh_scoring_accumulator_merge(&dst2, &c) == OSH_OK);
    ASSERT_TRUE(osh_scoring_accumulator_merge(&dst2, &b) == OSH_OK);

    for (i = 0; i < 3u; ++i) {
        ASSERT_TRUE(dst1.data[i] == dst2.data[i]);
    }

    osh_scoring_accumulator_free(&dst1);
    osh_scoring_accumulator_free(&dst2);
    osh_scoring_accumulator_free(&b);
    osh_scoring_accumulator_free(&c);
}

/* ---- Rejection: mismatched lengths --------------------------------------- */

static void test_merge_len_mismatch(void) {
    struct osh_scoring_accumulator a;
    struct osh_scoring_accumulator b;

    ASSERT_TRUE(osh_scoring_accumulator_alloc(&a, 4u, 0) == OSH_OK);
    ASSERT_TRUE(osh_scoring_accumulator_alloc(&b, 5u, 0) == OSH_OK);

    ASSERT_TRUE(osh_scoring_accumulator_merge(&a, &b) == OSH_EINVAL);

    osh_scoring_accumulator_free(&a);
    osh_scoring_accumulator_free(&b);
}

/* ---- Rejection: NULL operands -------------------------------------------- */

static void test_merge_null(void) {
    struct osh_scoring_accumulator a;

    ASSERT_TRUE(osh_scoring_accumulator_alloc(&a, 2u, 0) == OSH_OK);
    ASSERT_TRUE(osh_scoring_accumulator_merge(NULL, &a) == OSH_EINVAL);
    ASSERT_TRUE(osh_scoring_accumulator_merge(&a, NULL) == OSH_EINVAL);
    ASSERT_TRUE(osh_scoring_accumulator_merge(NULL, NULL) == OSH_EINVAL);
    osh_scoring_accumulator_free(&a);
}

int main(void) {
    test_merge_identity();
    test_merge_correctness();
    test_merge_commutative();
    test_merge_len_mismatch();
    test_merge_null();
    printf("All osh_scoring_accumulator_merge tests passed.\n");
    return 0;
}
