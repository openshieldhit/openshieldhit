#include <stddef.h>

#include "test_assert.h"
#include "transport/osh_run_control.h"

/* ---- Callback stubs ------------------------------------------------------ */

static int cb_always_stop(void *user) {
    (void) user;
    return 1;
}

/* Reads the stop request out of its caller-owned context (an int *). */
static int cb_from_flag(void *user) {
    return *(int const *) user;
}

/* ---- run_ctl_should_stop: defaults & NULL safety ------------------------- */

static void test_null_control_never_stops(void) {
    /* The un-controlled baseline: no policy means the run is never interrupted,
     * whatever the elapsed time or completed count. */
    ASSERT_TRUE(run_ctl_should_stop(NULL, 0.0, 0u) == 0);
    ASSERT_TRUE(run_ctl_should_stop(NULL, 1.0e9, 1000000u) == 0);
}

static void test_init_disables_everything(void) {
    struct osh_run_control ctl;
    osh_run_control_init(&ctl);

    ASSERT_TRUE(ctl.wall_budget_s == 0.0);
    ASSERT_TRUE(ctl.dump_every_s == 0.0);
    ASSERT_TRUE(ctl.dump_every_primaries == 0u);
    ASSERT_TRUE(ctl.should_stop == NULL);
    ASSERT_TRUE(ctl.should_stop_user == NULL);

    /* Zeroed control: unlimited budget, no callback → never stops. */
    ASSERT_TRUE(run_ctl_should_stop(&ctl, 0.0, 0u) == 0);
    ASSERT_TRUE(run_ctl_should_stop(&ctl, 123456.0, 42u) == 0);
}

/* ---- Wall-time budget boundary ------------------------------------------- */

static void test_wall_budget_boundary(void) {
    struct osh_run_control ctl;
    osh_run_control_init(&ctl);
    ctl.wall_budget_s = 10.0;

    /* Strictly below the budget: keep going. */
    ASSERT_TRUE(run_ctl_should_stop(&ctl, 0.0, 0u) == 0);
    ASSERT_TRUE(run_ctl_should_stop(&ctl, 9.999, 0u) == 0);

    /* At or past the budget: stop (>= boundary is inclusive). */
    ASSERT_TRUE(run_ctl_should_stop(&ctl, 10.0, 0u) == 1);
    ASSERT_TRUE(run_ctl_should_stop(&ctl, 10.0001, 0u) == 1);
    ASSERT_TRUE(run_ctl_should_stop(&ctl, 1.0e6, 0u) == 1);
}

static void test_zero_budget_is_unlimited(void) {
    struct osh_run_control ctl;
    osh_run_control_init(&ctl);
    ctl.wall_budget_s = 0.0; /* 0 means "no limit", not "stop immediately" */

    ASSERT_TRUE(run_ctl_should_stop(&ctl, 0.0, 0u) == 0);
    ASSERT_TRUE(run_ctl_should_stop(&ctl, 1.0e9, 0u) == 0);
}

/* ---- should_stop callback ------------------------------------------------ */

static void test_should_stop_callback(void) {
    struct osh_run_control ctl;
    int flag = 0;

    osh_run_control_init(&ctl);
    ctl.should_stop = cb_from_flag;
    ctl.should_stop_user = &flag;

    /* Callback returns 0, no budget: keep going. */
    ASSERT_TRUE(run_ctl_should_stop(&ctl, 5.0, 7u) == 0);

    /* Callback returns non-zero: stop regardless of elapsed time. */
    flag = 1;
    ASSERT_TRUE(run_ctl_should_stop(&ctl, 0.0, 0u) == 1);

    /* Back to 0 returns to "keep going" — the library re-asks every call. */
    flag = 0;
    ASSERT_TRUE(run_ctl_should_stop(&ctl, 0.0, 0u) == 0);
}

static void test_budget_and_callback_combine_or(void) {
    struct osh_run_control ctl;
    int flag = 0;

    osh_run_control_init(&ctl);
    ctl.wall_budget_s = 100.0;
    ctl.should_stop = cb_from_flag;
    ctl.should_stop_user = &flag;

    /* Neither trigger met. */
    ASSERT_TRUE(run_ctl_should_stop(&ctl, 1.0, 0u) == 0);

    /* Callback alone, well under budget. */
    flag = 1;
    ASSERT_TRUE(run_ctl_should_stop(&ctl, 1.0, 0u) == 1);

    /* Budget alone, callback returns 0. */
    flag = 0;
    ASSERT_TRUE(run_ctl_should_stop(&ctl, 100.0, 0u) == 1);
}

static void test_callback_without_user(void) {
    struct osh_run_control ctl;
    osh_run_control_init(&ctl);
    ctl.should_stop = cb_always_stop; /* ignores user; user stays NULL */

    ASSERT_TRUE(run_ctl_should_stop(&ctl, 0.0, 0u) == 1);
}

/* ---- start: records the monotonic baseline ------------------------------- */

static void test_start_records_baseline(void) {
    struct osh_run_control ctl;
    osh_run_control_init(&ctl);
    osh_run_control_start(&ctl, 1234.5);

    ASSERT_TRUE(ctl.t_start == 1234.5);
    ASSERT_TRUE(ctl.last_dump_s == 1234.5);
    ASSERT_TRUE(ctl.last_dump_primaries == 0u);
}

/* ---- NULL-safe lifecycle calls ------------------------------------------- */

static void test_null_lifecycle_calls_are_noops(void) {
    /* Must not crash. */
    osh_run_control_init(NULL);
    osh_run_control_start(NULL, 1.0);
}

int main(void) {
    test_null_control_never_stops();
    test_init_disables_everything();
    test_wall_budget_boundary();
    test_zero_budget_is_unlimited();
    test_should_stop_callback();
    test_budget_and_callback_combine_or();
    test_callback_without_user();
    test_start_records_baseline();
    test_null_lifecycle_calls_are_noops();
    return 0;
}
