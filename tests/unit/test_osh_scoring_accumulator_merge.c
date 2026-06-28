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
 *   test_merge_shape_mismatch  — mismatched optional-array presence is rejected
 *   test_zero                  — zero() clears every allocated array, keeps len
 *   test_rescale               — rescale() scales data in place, NULL-safe
 *   test_finalize_average      — finalize_average() divides data/data2, guards eps
 *   test_merge_variance_two_batches        — Welford M2 merge: exact dyadic case
 *   test_merge_variance_unequal_single_pass — unequal batches match single-pass
 *                                             weighted M2, order-independent
 *   test_merge_variance_empty_identity     — a weight-0 batch is the merge identity
 *   test_merge_variance_inconsistent_rejected — malformed variance bookkeeping → EINVAL
 */

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include "openshieldhit/status.h"
#include "scoring/runtime/osh_scoring_accumulator.h"
#include "test_assert.h"

/* Relative-tolerance compare: the Welford/Schubert-Gertz merge divides sums by
 * weights, so results are bit-exact only for dyadic test vectors; non-dyadic
 * cases and cross-order comparisons use this. */
static int approx(double a, double b) {
    double const d = fabs(a - b);
    return d <= 1.0e-9 * (1.0 + fabs(b));
}

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

/* ---- Reuse: zero() clears every allocated array but keeps len ------------ */

static void test_zero(void) {
    struct osh_scoring_accumulator acc;
    size_t i;

    ASSERT_TRUE(osh_scoring_accumulator_alloc(&acc, 4u, 1) == OSH_OK);
    for (i = 0; i < 4u; ++i) {
        acc.data[i] = (double) (i + 1);
        acc.data2[i] = (double) (i + 5);
    }

    osh_scoring_accumulator_zero(&acc);
    ASSERT_TRUE(acc.len == 4u); /* zeroing clears values, not the allocation */
    for (i = 0; i < 4u; ++i) {
        ASSERT_TRUE(acc.data[i] == 0.0);
        ASSERT_TRUE(acc.data2[i] == 0.0);
    }

    /* No-op guards: NULL and a zeroed/empty struct must not crash. */
    osh_scoring_accumulator_zero(NULL);
    {
        struct osh_scoring_accumulator empty = {0};
        osh_scoring_accumulator_zero(&empty);
    }

    osh_scoring_accumulator_free(&acc);
}

/* ---- Rejection: mismatched optional-array presence ----------------------- */

static void test_merge_shape_mismatch(void) {
    struct osh_scoring_accumulator with2;
    struct osh_scoring_accumulator without2;

    ASSERT_TRUE(osh_scoring_accumulator_alloc(&with2, 3u, 1) == OSH_OK);    /* has data2 */
    ASSERT_TRUE(osh_scoring_accumulator_alloc(&without2, 3u, 0) == OSH_OK); /* no data2 */

    /* data2 present on exactly one side must be rejected in either direction. */
    ASSERT_TRUE(osh_scoring_accumulator_merge(&with2, &without2) == OSH_EINVAL);
    ASSERT_TRUE(osh_scoring_accumulator_merge(&without2, &with2) == OSH_EINVAL);

    osh_scoring_accumulator_free(&with2);
    osh_scoring_accumulator_free(&without2);
}

/* ---- rescale(): scale the primary array in place ------------------------- */

static void test_rescale(void) {
    struct osh_scoring_accumulator acc;
    size_t i;

    ASSERT_TRUE(osh_scoring_accumulator_alloc(&acc, 3u, 0) == OSH_OK);
    for (i = 0; i < 3u; ++i) {
        acc.data[i] = (double) (i + 1); /* 1,2,3 */
    }

    osh_scoring_accumulator_rescale(&acc, 10.0);
    ASSERT_TRUE(acc.data[0] == 10.0);
    ASSERT_TRUE(acc.data[1] == 20.0);
    ASSERT_TRUE(acc.data[2] == 30.0);

    osh_scoring_accumulator_rescale(NULL, 2.0); /* NULL is a no-op, must not crash */

    osh_scoring_accumulator_free(&acc);
}

/* ---- finalize_average(): data/data2 with an eps guard -------------------- */

static void test_finalize_average(void) {
    struct osh_scoring_accumulator acc;
    struct osh_scoring_accumulator nodata2;

    ASSERT_TRUE(osh_scoring_accumulator_alloc(&acc, 3u, 1) == OSH_OK);
    acc.data[0] = 10.0;
    acc.data2[0] = 2.0; /* -> 5 */
    acc.data[1] = 9.0;
    acc.data2[1] = 3.0; /* -> 3 */
    acc.data[2] = 7.0;
    acc.data2[2] = 0.0; /* weight ~0 -> 0 */

    ASSERT_TRUE(osh_scoring_accumulator_finalize_average(&acc, 1.0e-300) == OSH_OK);
    ASSERT_TRUE(acc.data[0] == 5.0);
    ASSERT_TRUE(acc.data[1] == 3.0);
    ASSERT_TRUE(acc.data[2] == 0.0);

    /* Missing denominator array must be rejected, not dereferenced. */
    ASSERT_TRUE(osh_scoring_accumulator_alloc(&nodata2, 3u, 0) == OSH_OK);
    ASSERT_TRUE(osh_scoring_accumulator_finalize_average(&nodata2, 1.0e-300) == OSH_EINVAL);

    osh_scoring_accumulator_free(&acc);
    osh_scoring_accumulator_free(&nodata2);
}

/* ---- Variance (Welford M2) merge: Schubert-Gertz cross-term ----------------
 *
 * The variance feature is unwired, so osh_scoring_accumulator_alloc() does not
 * create the M2 arrays.  These tests build single-batch accumulators by hand:
 * data + a calloc'd data_var (M2 = 0 for one observation), with weight = the
 * batch's history count and nbatch = 1.  osh_scoring_accumulator_free() releases
 * data_var, so teardown is the usual free().
 */

/* Build a one-observation batch: data = per-bin sum, M2 = 0, weight = w. */
static void make_batch(struct osh_scoring_accumulator *acc, size_t len, double const *data, double weight) {
    size_t i;
    ASSERT_TRUE(osh_scoring_accumulator_alloc(acc, len, 0) == OSH_OK);
    acc->data_var = (double *) calloc(len ? len : 1u, sizeof(double));
    ASSERT_TRUE(acc->data_var != NULL);
    for (i = 0; i < len; ++i) {
        acc->data[i] = data[i];
    }
    acc->weight = weight;
    acc->nbatch = 1u;
}

/* Two dyadic batches merge to the exact hand-computed M2 and additive fields. */
static void test_merge_variance_two_batches(void) {
    struct osh_scoring_accumulator a;
    struct osh_scoring_accumulator b;
    double const da[2] = {8.0, 4.0};  /* means 4, 2 at weight 2 */
    double const db[2] = {16.0, 4.0}; /* means 8, 2 at weight 2 */

    make_batch(&a, 2u, da, 2.0);
    make_batch(&b, 2u, db, 2.0);

    ASSERT_TRUE(osh_scoring_accumulator_merge(&a, &b) == OSH_OK);

    /* bin0: δ = 8 − 4 = 4, M2 = 4² · (2·2/4) = 16.  bin1: δ = 0, M2 = 0. */
    ASSERT_TRUE(a.data_var[0] == 16.0);
    ASSERT_TRUE(a.data_var[1] == 0.0);
    /* Additive fields. */
    ASSERT_TRUE(a.data[0] == 24.0);
    ASSERT_TRUE(a.data[1] == 8.0);
    ASSERT_TRUE(a.weight == 4.0);
    ASSERT_TRUE(a.nbatch == 2u);

    osh_scoring_accumulator_free(&a);
    osh_scoring_accumulator_free(&b);
}

/* Three unequal-size batches: pairwise merge reproduces the single-pass weighted
 * M2, and the result is order-independent (within FP tolerance). */
static void test_merge_variance_unequal_single_pass(void) {
    struct osh_scoring_accumulator a, b, c;    /* order 1: (a⊕b)⊕c */
    struct osh_scoring_accumulator a2, b2, c2; /* order 2: (c⊕b)⊕a */
    /* weights 2,3,5; bin0 means 5,10,1; bin1 means 2,1,1. */
    double const da[2] = {10.0, 4.0};
    double const db[2] = {30.0, 3.0};
    double const dc[2] = {5.0, 5.0};
    /* Single-pass weighted M2 (hand-computed):
     *   bin0: M = 45/10 = 4.5; M2 = 2·.25 + 3·30.25 + 5·12.25 = 152.5
     *   bin1: M = 12/10 = 1.2; M2 = 2·.64 + 3·.04 + 5·.04 = 1.6           */
    double const expect_m2[2] = {152.5, 1.6};

    make_batch(&a, 2u, da, 2.0);
    make_batch(&b, 2u, db, 3.0);
    make_batch(&c, 2u, dc, 5.0);
    ASSERT_TRUE(osh_scoring_accumulator_merge(&a, &b) == OSH_OK);
    ASSERT_TRUE(osh_scoring_accumulator_merge(&a, &c) == OSH_OK);

    ASSERT_TRUE(approx(a.data_var[0], expect_m2[0]));
    ASSERT_TRUE(approx(a.data_var[1], expect_m2[1]));
    ASSERT_TRUE(a.data[0] == 45.0); /* sums are exact */
    ASSERT_TRUE(a.weight == 10.0);
    ASSERT_TRUE(a.nbatch == 3u);

    /* Reverse fold order must agree within tolerance (FP, not bit-exact). */
    make_batch(&a2, 2u, da, 2.0);
    make_batch(&b2, 2u, db, 3.0);
    make_batch(&c2, 2u, dc, 5.0);
    ASSERT_TRUE(osh_scoring_accumulator_merge(&c2, &b2) == OSH_OK);
    ASSERT_TRUE(osh_scoring_accumulator_merge(&c2, &a2) == OSH_OK);

    ASSERT_TRUE(approx(c2.data_var[0], a.data_var[0]));
    ASSERT_TRUE(approx(c2.data_var[1], a.data_var[1]));
    ASSERT_TRUE(c2.weight == 10.0);
    ASSERT_TRUE(c2.nbatch == 3u);

    osh_scoring_accumulator_free(&a);
    osh_scoring_accumulator_free(&b);
    osh_scoring_accumulator_free(&c);
    osh_scoring_accumulator_free(&a2);
    osh_scoring_accumulator_free(&b2);
    osh_scoring_accumulator_free(&c2);
}

/* An empty batch (weight 0) is the merge identity in both directions. */
static void test_merge_variance_empty_identity(void) {
    struct osh_scoring_accumulator empty;
    struct osh_scoring_accumulator b;
    double const zero[2] = {0.0, 0.0};
    double const db[2] = {30.0, 6.0};

    /* empty ⊕ b  ⇒  empty becomes b. */
    make_batch(&empty, 2u, zero, 0.0);
    empty.nbatch = 0u; /* no batches yet */
    make_batch(&b, 2u, db, 3.0);

    ASSERT_TRUE(osh_scoring_accumulator_merge(&empty, &b) == OSH_OK);
    ASSERT_TRUE(empty.data_var[0] == 0.0); /* single observation ⇒ M2 still 0 */
    ASSERT_TRUE(empty.data[0] == 30.0);
    ASSERT_TRUE(empty.weight == 3.0);
    ASSERT_TRUE(empty.nbatch == 1u);

    /* b ⊕ empty(weight 0)  ⇒  b unchanged. */
    {
        struct osh_scoring_accumulator e2;
        make_batch(&e2, 2u, zero, 0.0);
        e2.nbatch = 0u;
        ASSERT_TRUE(osh_scoring_accumulator_merge(&empty, &e2) == OSH_OK);
        ASSERT_TRUE(empty.data_var[0] == 0.0);
        ASSERT_TRUE(empty.data[0] == 30.0);
        ASSERT_TRUE(empty.weight == 3.0);
        ASSERT_TRUE(empty.nbatch == 1u);
        osh_scoring_accumulator_free(&e2);
    }

    osh_scoring_accumulator_free(&empty);
    osh_scoring_accumulator_free(&b);
}

/* Inconsistent variance bookkeeping is rejected rather than silently mis-merged. */
static void test_merge_variance_inconsistent_rejected(void) {
    struct osh_scoring_accumulator a;
    struct osh_scoring_accumulator b;
    double const d[2] = {10.0, 4.0};

    /* weight == 0 but nbatch != 0 is a malformed batch. */
    make_batch(&a, 2u, d, 2.0);
    make_batch(&b, 2u, d, 2.0);
    a.weight = 0.0; /* nbatch stays 1 -> inconsistent */
    ASSERT_TRUE(osh_scoring_accumulator_merge(&a, &b) == OSH_EINVAL);
    ASSERT_TRUE(osh_scoring_accumulator_merge(&b, &a) == OSH_EINVAL); /* src side too */
    osh_scoring_accumulator_free(&a);
    osh_scoring_accumulator_free(&b);

    /* data2_var present without its companion data2 array. */
    {
        struct osh_scoring_accumulator c;
        struct osh_scoring_accumulator e;
        make_batch(&c, 2u, d, 2.0); /* no data2 (alloc want_data2=0) */
        make_batch(&e, 2u, d, 2.0);
        c.data2_var = (double *) calloc(2u, sizeof(double));
        e.data2_var = (double *) calloc(2u, sizeof(double));
        ASSERT_TRUE(c.data2_var != NULL && e.data2_var != NULL);
        ASSERT_TRUE(osh_scoring_accumulator_merge(&c, &e) == OSH_EINVAL);
        osh_scoring_accumulator_free(&c);
        osh_scoring_accumulator_free(&e);
    }
}

int main(void) {
    test_merge_identity();
    test_merge_correctness();
    test_merge_commutative();
    test_merge_len_mismatch();
    test_merge_null();
    test_merge_shape_mismatch();
    test_zero();
    test_rescale();
    test_finalize_average();
    test_merge_variance_two_batches();
    test_merge_variance_unequal_single_pass();
    test_merge_variance_empty_identity();
    test_merge_variance_inconsistent_rejected();
    printf("All osh_scoring_accumulator tests passed.\n");
    return 0;
}
