#include <string.h>

#include "test_assert.h"
#include "transport/osh_checkpoint_policy.h"

/* ---- init + is_final_only ------------------------------------------------ */

static void test_init_is_final_only_and_exact(void) {
    struct osh_checkpoint_policy p;

    /* A garbage-filled struct must come back fully defaulted. */
    memset(&p, 0xAB, sizeof(p));
    osh_checkpoint_policy_init(&p);

    ASSERT_TRUE(p.mode == OSH_PARTIAL_NONE);
    ASSERT_TRUE(p.completeness == OSH_PARTIAL_EXACT);
    ASSERT_TRUE(p.every_s == 0.0);
    ASSERT_TRUE(p.every_primaries == 0u);
    ASSERT_TRUE(p.batch == 0u);
    ASSERT_TRUE(p.write_files == 0);

    ASSERT_TRUE(osh_checkpoint_policy_is_final_only(&p) == 1);

    /* init(NULL) is a no-op, not a crash. */
    osh_checkpoint_policy_init(NULL);
}

static void test_is_final_only_null_and_modes(void) {
    struct osh_checkpoint_policy p;

    /* A NULL policy is the un-configured baseline: final-only. */
    ASSERT_TRUE(osh_checkpoint_policy_is_final_only(NULL) == 1);

    osh_checkpoint_policy_init(&p);
    ASSERT_TRUE(osh_checkpoint_policy_is_final_only(&p) == 1);

    p.mode = OSH_PARTIAL_LIVE;
    ASSERT_TRUE(osh_checkpoint_policy_is_final_only(&p) == 0);
}

/* ---- next_batch_size ----------------------------------------------------- */

static void test_next_batch_size_final_only_is_whole_run(void) {
    struct osh_checkpoint_policy p;

    osh_checkpoint_policy_init(&p);

    /* Final-only (and NULL) both return the whole remaining range: the K = nstat
     * fast path.  A stray cadence is ignored while mode is NONE. */
    ASSERT_TRUE(osh_checkpoint_next_batch_size(&p, 0.0, 12345u) == 12345u);
    ASSERT_TRUE(osh_checkpoint_next_batch_size(NULL, 0.0, 999u) == 999u);

    p.every_primaries = 10u; /* ignored: mode is still NONE */
    ASSERT_TRUE(osh_checkpoint_next_batch_size(&p, 0.0, 999u) == 999u);
}

static void test_next_batch_size_zero_remaining(void) {
    struct osh_checkpoint_policy p;

    osh_checkpoint_policy_init(&p);
    p.mode = OSH_PARTIAL_LIVE;
    p.every_primaries = 100u;

    ASSERT_TRUE(osh_checkpoint_next_batch_size(&p, 0.0, 0u) == 0u);
    ASSERT_TRUE(osh_checkpoint_next_batch_size(NULL, 0.0, 0u) == 0u);
}

static void test_next_batch_size_count_cadence(void) {
    struct osh_checkpoint_policy p;

    osh_checkpoint_policy_init(&p);
    p.mode = OSH_PARTIAL_LIVE;
    p.every_primaries = 250u;

    /* Full-size batches while plenty remains; clamped to remaining at the tail. */
    ASSERT_TRUE(osh_checkpoint_next_batch_size(&p, 0.0, 1000u) == 250u);
    ASSERT_TRUE(osh_checkpoint_next_batch_size(&p, 0.0, 250u) == 250u);
    ASSERT_TRUE(osh_checkpoint_next_batch_size(&p, 0.0, 100u) == 100u);
    ASSERT_TRUE(osh_checkpoint_next_batch_size(&p, 0.0, 1u) == 1u);
}

static void test_next_batch_size_explicit_batch_overrides_cadence(void) {
    struct osh_checkpoint_policy p;

    osh_checkpoint_policy_init(&p);
    p.mode = OSH_PARTIAL_LIVE;
    p.every_primaries = 250u; /* decoupled from... */
    p.batch = 64u;            /* ...the explicit batch size, which wins. */

    ASSERT_TRUE(osh_checkpoint_next_batch_size(&p, 0.0, 1000u) == 64u);
    ASSERT_TRUE(osh_checkpoint_next_batch_size(&p, 0.0, 40u) == 40u); /* clamp */
}

static void test_next_batch_size_time_cadence_adaptive(void) {
    struct osh_checkpoint_policy p;

    osh_checkpoint_policy_init(&p);
    p.mode = OSH_PARTIAL_LIVE;
    p.every_s = 2.0;

    /* K ≈ round(rate × cadence).  100 pps × 2 s = 200. */
    ASSERT_TRUE(osh_checkpoint_next_batch_size(&p, 100.0, 1000u) == 200u);
    /* Clamp to remaining when the interval outruns the run. */
    ASSERT_TRUE(osh_checkpoint_next_batch_size(&p, 100.0, 150u) == 150u);
    /* A tiny rate still yields at least one primary of progress. */
    ASSERT_TRUE(osh_checkpoint_next_batch_size(&p, 0.1, 1000u) == 1u);
    /* No measured rate yet → span the rest (stay exact, keep progressing). */
    ASSERT_TRUE(osh_checkpoint_next_batch_size(&p, 0.0, 1000u) == 1000u);
}

static void test_next_batch_size_live_without_cadence_is_whole_run(void) {
    struct osh_checkpoint_policy p;

    osh_checkpoint_policy_init(&p);
    p.mode = OSH_PARTIAL_LIVE; /* LIVE but no cadence configured at all */

    ASSERT_TRUE(osh_checkpoint_next_batch_size(&p, 0.0, 777u) == 777u);
}

/* ---- completeness label -------------------------------------------------- */

static void test_completeness_label(void) {
    ASSERT_TRUE(strcmp(osh_checkpoint_completeness_label(OSH_PARTIAL_EXACT), "exact") == 0);
    ASSERT_TRUE(strcmp(osh_checkpoint_completeness_label(OSH_PARTIAL_APPROX), "families_pending") == 0);
}

/* ---- Entry point --------------------------------------------------------- */

int main(void) {
    test_init_is_final_only_and_exact();
    test_is_final_only_null_and_modes();
    test_next_batch_size_final_only_is_whole_run();
    test_next_batch_size_zero_remaining();
    test_next_batch_size_count_cadence();
    test_next_batch_size_explicit_batch_overrides_cadence();
    test_next_batch_size_time_cadence_adaptive();
    test_next_batch_size_live_without_cadence_is_whole_run();
    test_completeness_label();
    return 0;
}
