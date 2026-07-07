/*
 * Unit tests for osh_scoring_finalize_errors() — the batch-means standard-error
 * finalize (issue #209).  It converts each variance page's Welford M2 into a
 * per-bin *relative* standard error, stored back into data_var.  Pure arithmetic
 * over hand-built runtimes, so identical on every OS.
 *
 * Coverage:
 *   test_additive_relative_error  — additive estimator: rel = sqrt(M2·W/(B-1))/|Σx|
 *   test_ratio_relative_error     — AVER estimator: rel = sqrt(rel_num² + rel_den²)
 *   test_one_batch_zero_error     — B < 2 (zero d.o.f.) → error zeroed
 *   test_empty_bin_zero_error     — a zero-sum bin has no defined relative error → 0
 *   test_variance_off_noop        — pages without M2 arrays are left untouched
 *   test_null_rejected            — NULL runtime → OSH_EINVAL
 */

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "openshieldhit/status.h"
#include "scoring/runtime/osh_scoring_accumulator.h"
#include "scoring/runtime/osh_scoring_postprocess.h"
#include "scoring/runtime/osh_scoring_runtime.h"
#include "test_assert.h"

static int approx(double a, double b) {
    double const d = fabs(a - b);
    return d <= 1.0e-9 * (1.0 + fabs(b));
}

/* One page with an accumulator that (optionally) tracks variance. */
static void make_var_page(struct osh_scoring_page_runtime *p,
                          enum osh_scoring_score_kind kind,
                          size_t len,
                          int want_data2,
                          int want_variance) {
    memset(p, 0, sizeof(*p));
    ASSERT_TRUE(osh_scoring_accumulator_alloc_variance(&p->acc, len, want_data2, want_variance) == OSH_OK);
    p->len = len;
    p->score_kind = kind;
    p->has_data2 = (char) (want_data2 ? 1 : 0);
    p->variance = (char) (want_variance ? 1 : 0);
}

static void make_runtime(struct osh_scoring_runtime *rt, struct osh_scoring_page_runtime *pages, size_t npages) {
    memset(rt, 0, sizeof(*rt));
    rt->pages = pages;
    rt->npages = npages;
}

/* Additive estimator (DOSE/FLUENCE/…): the relative SE of the accumulated sum. */
static void test_additive_relative_error(void) {
    struct osh_scoring_page_runtime page;
    struct osh_scoring_runtime rt;

    make_var_page(&page, OSH_SCORING_SCORE_DOSE, 2u, 0, 1);
    /* Two batches (B=2, W=4): bin0 sum=24, M2=16 → rel = sqrt(16·4/1)/24 = 8/24;
     * bin1 sum=8, M2=0 → rel = 0. */
    page.acc.data[0] = 24.0;
    page.acc.data[1] = 8.0;
    page.acc.data_var[0] = 16.0;
    page.acc.data_var[1] = 0.0;
    page.acc.weight = 4.0;
    page.acc.nbatch = 2u;
    make_runtime(&rt, &page, 1u);

    ASSERT_TRUE(osh_scoring_finalize_errors(&rt) == OSH_OK);
    ASSERT_TRUE(approx(page.acc.data_var[0], 8.0 / 24.0));
    ASSERT_TRUE(approx(page.acc.data_var[1], 0.0));

    osh_scoring_accumulator_free(&page.acc);
}

/* AVER ratio estimator (DLET/TLET/Qeff): quadrature of numerator + denominator. */
static void test_ratio_relative_error(void) {
    struct osh_scoring_page_runtime page;
    struct osh_scoring_runtime rt;
    double rel_num;
    double rel_den;
    double expect;

    make_var_page(&page, OSH_SCORING_SCORE_DLET, 1u, 1, 1);
    /* num sum=24, M2=16 → rel_num = 8/24; den sum=8, M2=1 → rel_den = sqrt(4)/8 = 2/8. */
    page.acc.data[0] = 24.0;
    page.acc.data_var[0] = 16.0;
    page.acc.data2[0] = 8.0;
    page.acc.data2_var[0] = 1.0;
    page.acc.weight = 4.0;
    page.acc.nbatch = 2u;
    make_runtime(&rt, &page, 1u);

    rel_num = 8.0 / 24.0;
    rel_den = 2.0 / 8.0;
    expect = sqrt(rel_num * rel_num + rel_den * rel_den);

    ASSERT_TRUE(osh_scoring_finalize_errors(&rt) == OSH_OK);
    ASSERT_TRUE(approx(page.acc.data_var[0], expect));

    osh_scoring_accumulator_free(&page.acc);
}

/* Fewer than two batches has zero degrees of freedom → the error column is zeroed
 * (raw M2 must not leak through). */
static void test_one_batch_zero_error(void) {
    struct osh_scoring_page_runtime page;
    struct osh_scoring_runtime rt;

    make_var_page(&page, OSH_SCORING_SCORE_DOSE, 2u, 0, 1);
    page.acc.data[0] = 10.0;
    page.acc.data[1] = 5.0;
    page.acc.data_var[0] = 7.0; /* stray M2 that must be cleared */
    page.acc.data_var[1] = 3.0;
    page.acc.weight = 2.0;
    page.acc.nbatch = 1u; /* single batch */
    make_runtime(&rt, &page, 1u);

    ASSERT_TRUE(osh_scoring_finalize_errors(&rt) == OSH_OK);
    ASSERT_TRUE(page.acc.data_var[0] == 0.0);
    ASSERT_TRUE(page.acc.data_var[1] == 0.0);

    osh_scoring_accumulator_free(&page.acc);
}

/* A bin whose sum is zero (no deposits) has no relative error defined → 0. */
static void test_empty_bin_zero_error(void) {
    struct osh_scoring_page_runtime page;
    struct osh_scoring_runtime rt;

    make_var_page(&page, OSH_SCORING_SCORE_DOSE, 2u, 0, 1);
    page.acc.data[0] = 0.0; /* empty bin */
    page.acc.data[1] = 12.0;
    page.acc.data_var[0] = 5.0; /* M2 present but sum 0 */
    page.acc.data_var[1] = 9.0;
    page.acc.weight = 3.0;
    page.acc.nbatch = 3u;
    make_runtime(&rt, &page, 1u);

    ASSERT_TRUE(osh_scoring_finalize_errors(&rt) == OSH_OK);
    ASSERT_TRUE(page.acc.data_var[0] == 0.0);
    ASSERT_TRUE(approx(page.acc.data_var[1], sqrt(9.0 * 3.0 / 2.0) / 12.0));

    osh_scoring_accumulator_free(&page.acc);
}

/* Variance off (no M2 arrays): finalize is a no-op and leaves data untouched. */
static void test_variance_off_noop(void) {
    struct osh_scoring_page_runtime page;
    struct osh_scoring_runtime rt;

    make_var_page(&page, OSH_SCORING_SCORE_DOSE, 2u, 0, 0);
    page.acc.data[0] = 4.0;
    page.acc.data[1] = 5.0;
    make_runtime(&rt, &page, 1u);

    ASSERT_TRUE(page.acc.data_var == NULL);
    ASSERT_TRUE(osh_scoring_finalize_errors(&rt) == OSH_OK);
    ASSERT_TRUE(page.acc.data[0] == 4.0 && page.acc.data[1] == 5.0);
    ASSERT_TRUE(page.acc.data_var == NULL);

    osh_scoring_accumulator_free(&page.acc);
}

static void test_null_rejected(void) {
    ASSERT_TRUE(osh_scoring_finalize_errors(NULL) == OSH_EINVAL);
}

int main(void) {
    test_additive_relative_error();
    test_ratio_relative_error();
    test_one_batch_zero_error();
    test_empty_bin_zero_error();
    test_variance_off_noop();
    test_null_rejected();
    printf("All osh_scoring_finalize_errors tests passed.\n");
    return 0;
}
