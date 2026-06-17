/*
 * Unit tests for hot-path scoring filter evaluation
 * (osh_scoring_page_passes_filters).
 *
 * Focus: the NPRIM (primary-index) filter must compare the full uint64_t
 * prim_idx, not a value truncated to 32 bits.  This guards the fix that
 * accompanies widening prim_idx to uint64_t — a history index above UINT_MAX
 * must still filter correctly.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_step.h"
#include "particle/osh_particle.h"
#include "scoring/runtime/osh_scoring_defs.h"
#include "scoring/runtime/osh_scoring_filter_runtime.h"
#include "scoring/runtime/osh_scoring_output_runtime.h"

#define ASSERT_TRUE(cond)                                                                                              \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            fprintf(stderr, "ASSERT FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__);                                 \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

/* Evaluate a single NPRIM rule (field op value) against a history index. */
static int nprim_passes(uint64_t prim_idx, enum osh_scoring_filter_op op, double value) {
    struct osh_scoring_filter_runtime_rule rule;
    struct osh_scoring_page_runtime page;
    struct particle part;
    struct step st;

    memset(&page, 0, sizeof(page));
    memset(&part, 0, sizeof(part));
    memset(&st, 0, sizeof(st));

    rule.field = OSH_SCORING_FILTER_FIELD_NPRIM;
    rule.op = op;
    rule.value = value;

    page.flat_rules = &rule;
    page.nflat_rules = 1u;

    st.prim_idx = prim_idx;

    return osh_scoring_page_passes_filters(&page, &part, &st);
}

static void test_nprim_filter_small_indices(void) {
    /* Baseline behaviour for indices that fit comfortably in 32 bits. */
    ASSERT_TRUE(nprim_passes(5u, OSH_SCORING_FILTER_OP_EQ, 5.0));
    ASSERT_TRUE(!nprim_passes(5u, OSH_SCORING_FILTER_OP_EQ, 6.0));
    ASSERT_TRUE(nprim_passes(5u, OSH_SCORING_FILTER_OP_LT, 6.0));
    ASSERT_TRUE(nprim_passes(7u, OSH_SCORING_FILTER_OP_GE, 7.0));
    ASSERT_TRUE(nprim_passes(7u, OSH_SCORING_FILTER_OP_NE, 8.0));
}

static void test_nprim_filter_above_uint32_max(void) {
    /* 5e9 and 6e9 both exceed UINT_MAX (4294967295).  A 32-bit truncating
     * comparison would wrap these and give wrong answers; the uint64 path
     * must order them correctly. */
    uint64_t const idx = 5000000000ULL;
    double const four_billion = 4000000000.0;
    double const six_billion = 6000000000.0;

    ASSERT_TRUE(nprim_passes(idx, OSH_SCORING_FILTER_OP_GT, four_billion));  /* 5e9 > 4e9 */
    ASSERT_TRUE(!nprim_passes(idx, OSH_SCORING_FILTER_OP_GT, six_billion));  /* 5e9 < 6e9 */
    ASSERT_TRUE(nprim_passes(idx, OSH_SCORING_FILTER_OP_LT, six_billion));   /* 5e9 < 6e9 */
    ASSERT_TRUE(nprim_passes(idx, OSH_SCORING_FILTER_OP_EQ, 5000000000.0));  /* exact */
    ASSERT_TRUE(!nprim_passes(idx, OSH_SCORING_FILTER_OP_EQ, four_billion)); /* distinct */

    /* A truncating comparison would alias idx (5e9 mod 2^32 = 705032704)
     * with a small threshold; assert that does NOT happen. */
    ASSERT_TRUE(!nprim_passes(idx, OSH_SCORING_FILTER_OP_EQ, 705032704.0));
}

int main(void) {
    test_nprim_filter_small_indices();
    test_nprim_filter_above_uint32_max();
    return 0;
}
