#include "transport/osh_run_control.h"

void osh_run_control_init(struct osh_run_control *ctl) {
    if (!ctl) {
        return;
    }
    ctl->wall_budget_s = 0.0;
    ctl->dump_every_s = 0.0;
    ctl->dump_every_primaries = 0u;
    ctl->should_stop = NULL;
    ctl->should_stop_user = NULL;
    ctl->should_dump = NULL;
    ctl->should_dump_user = NULL;
    ctl->dump_sink = NULL;
    ctl->dump_shadow = NULL;
    ctl->dump_outputs = NULL;
    ctl->dump_noutputs = 0u;
    ctl->t_start = 0.0;
    ctl->last_dump_s = 0.0;
    ctl->last_dump_primaries = 0u;
}

void osh_run_control_start(struct osh_run_control *ctl, double t_now) {
    if (!ctl) {
        return;
    }
    ctl->t_start = t_now;
    /* Dump bookkeeping is run-relative (elapsed seconds), so the baseline is 0,
     * not t_now: the first time cadence fires once dump_every_s has *elapsed*. */
    ctl->last_dump_s = 0.0;
    ctl->last_dump_primaries = 0u;
}

int run_ctl_should_stop(const struct osh_run_control *ctl, double elapsed, size_t completed) {
    (void) completed; /* reserved for future count-based stop policies */

    if (!ctl) {
        return 0;
    }
    if (ctl->wall_budget_s > 0.0 && elapsed >= ctl->wall_budget_s) {
        return 1;
    }
    if (ctl->should_stop && ctl->should_stop(ctl->should_stop_user)) {
        return 1;
    }
    return 0;
}

int run_ctl_should_dump(struct osh_run_control *ctl, double elapsed, size_t completed) {
    int due = 0;

    if (!ctl) {
        return 0;
    }

    /* On-demand trigger first: the callback is edge-triggered (reads and clears
     * its flag), so it must be consulted every call to avoid stranding a pending
     * request — evaluate it up front rather than short-circuiting past it. */
    if (ctl->should_dump && ctl->should_dump(ctl->should_dump_user)) {
        due = 1;
    }
    if (ctl->dump_every_s > 0.0 && (elapsed - ctl->last_dump_s) >= ctl->dump_every_s) {
        due = 1;
    }
    if (ctl->dump_every_primaries > 0u && (completed - ctl->last_dump_primaries) >= ctl->dump_every_primaries) {
        due = 1;
    }

    if (due) {
        /* Rebase both cadences to this dump so the next one is measured from here.
         * Using the observed elapsed/completed (not last + cadence) keeps the
         * bookkeeping self-correcting when a batch overshoots the cadence. */
        ctl->last_dump_s = elapsed;
        ctl->last_dump_primaries = completed;
    }
    return due;
}

int run_ctl_has_scheduled_dump(struct osh_run_control const *ctl) {
    if (!ctl) {
        return 0;
    }
    return (ctl->dump_every_s > 0.0 || ctl->dump_every_primaries > 0u) ? 1 : 0;
}
