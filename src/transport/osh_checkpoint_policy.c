#include "transport/osh_checkpoint_policy.h"

/* Bootstrap probe for a time cadence before any throughput has been measured
 * (the very first batch): run 1/OSH_CHECKPOINT_TIME_PROBE_FRACTION of what
 * remains rather than the whole run, so we both measure a rate AND reach an
 * early checkpoint at which the first timed dump can fire.  From the second
 * batch on, the measured rate drives K ≈ rate × cadence and this is unused.
 * The fraction scales the probe with the run (small nstat → small probe) and is
 * floored at 1; the exact value only shifts *where* dumps land, never the scored
 * result (batching is RNG-invariant), and a time cadence is non-deterministic by
 * design (issue #195), so no test pins it. */
#define OSH_CHECKPOINT_TIME_PROBE_FRACTION 100u

void osh_checkpoint_policy_init(struct osh_checkpoint_policy *policy) {
    if (!policy) {
        return;
    }
    policy->mode = OSH_PARTIAL_NONE;
    policy->completeness = OSH_PARTIAL_EXACT;
    policy->every_s = 0.0;
    policy->every_primaries = 0u;
    policy->batch = 0u;
    policy->write_files = 0;
}

int osh_checkpoint_policy_is_final_only(struct osh_checkpoint_policy const *policy) {
    return (!policy || policy->mode == OSH_PARTIAL_NONE) ? 1 : 0;
}

size_t
osh_checkpoint_next_batch_size(struct osh_checkpoint_policy const *policy, double measured_rate_pps, size_t remaining) {
    size_t k;

    if (remaining == 0u) {
        return 0u;
    }
    /* FINAL-ONLY: one batch spanning the rest of the run.  This is the K = nstat
     * fast path — a single family pass, byte-for-byte identical to today. */
    if (osh_checkpoint_policy_is_final_only(policy)) {
        return remaining;
    }

    /* LIVE.  Pick the batch size in priority order: an explicit override, then
     * the deterministic count cadence, then the adaptive time cadence. */
    if (policy->batch != 0u) {
        k = policy->batch;
    } else if (policy->every_primaries != 0u) {
        k = policy->every_primaries;
    } else if (policy->every_s > 0.0) {
        if (measured_rate_pps > 0.0) {
            /* K ≈ rate × cadence, rounded to nearest, so a wall-time preview costs
             * the same per wall-second on any machine.  Guard the cast: clamp to
             * remaining before narrowing so a huge rate×time cannot wrap size_t. */
            double target = measured_rate_pps * policy->every_s + 0.5;
            if (target >= (double) remaining) {
                return remaining;
            }
            k = (target < 1.0) ? 1u : (size_t) target;
        } else {
            /* Time cadence, but throughput is not known yet (the first batch):
             * run a small bootstrap probe rather than the whole run, so we both
             * measure a rate and reach an early checkpoint for the first dump. */
            k = remaining / OSH_CHECKPOINT_TIME_PROBE_FRACTION;
            if (k == 0u) {
                k = 1u;
            }
        }
    } else {
        /* LIVE requested but no cadence is usable at all (mode set without a
         * cadence or an explicit batch): span the rest so the run still
         * progresses and the result stays exact. */
        return remaining;
    }

    /* Every branch above set k >= 1 (the cadence fields are non-zero here and the
     * time-cadence path floors at 1), so only the upper clamp to remaining is
     * needed to keep K within [1, remaining]. */
    if (k > remaining) {
        k = remaining;
    }
    return k;
}

char const *osh_checkpoint_completeness_label(enum osh_partial_completeness completeness) {
    return (completeness == OSH_PARTIAL_APPROX) ? "families_pending" : "exact";
}
