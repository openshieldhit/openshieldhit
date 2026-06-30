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
    ctl->t_start = 0.0;
    ctl->last_dump_s = 0.0;
    ctl->last_dump_primaries = 0u;
}

void osh_run_control_start(struct osh_run_control *ctl, double t_now) {
    if (!ctl) {
        return;
    }
    ctl->t_start = t_now;
    ctl->last_dump_s = t_now;
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
