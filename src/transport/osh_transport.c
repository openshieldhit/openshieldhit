#include "transport/osh_transport.h"

#include "common/osh_diag.h"
#include "transport/osh_checkpoint_policy.h"
#include "transport/osh_neutron_pool.h"
#include "transport/osh_transport_ion.h"
#include "transport/osh_transport_neutron.h"
#include "transport/osh_transport_photon.h"
#include "transport/osh_transport_scheduler.h"

/* ---- Forward declarations ------------------------------------------------ */

static enum osh_status run_families_over_range(struct osh_transport_context *transport_ctx,
                                               struct osh_beam_runtime *beam_rt,
                                               struct osh_gemca_runtime const *geom_rt,
                                               struct osh_material_runtime const *material_rt,
                                               struct osh_scoring_runtime *score_rt,
                                               size_t hist_lo,
                                               size_t hist_hi,
                                               int neutron_enabled,
                                               size_t *batch_completed_out);

static enum osh_status dispatch_transport_family(enum osh_transport_family family,
                                                 struct osh_transport_context *transport_ctx,
                                                 struct osh_beam_runtime *beam_rt,
                                                 struct osh_gemca_runtime const *geom_rt,
                                                 struct osh_material_runtime const *material_rt,
                                                 struct osh_scoring_runtime *score_rt,
                                                 size_t hist_lo,
                                                 size_t hist_hi,
                                                 size_t *ion_completed_out);

/*
 * Orchestrator for the minimal transport run.
 *
 * Two nested loops:
 *   • an OUTER batch loop (this function) that slices [0, nstat) into checkpoint
 *     batches [done, done+K) and brings the whole simulation to a family-complete,
 *     quiescent point at each boundary (issue #195);
 *   • an INNER family scheduler (run_families_over_range) that, for one batch,
 *     drains ions then the neutrons/fragments they banked, so the batch is
 *     *family-exact* before the checkpoint.
 *
 * The batch size K comes from the checkpoint policy.  FINAL-ONLY (the default —
 * a NULL policy or OSH_PARTIAL_NONE) yields a single batch of K = nstat, so the
 * outer loop runs exactly once and this path is byte-for-byte identical to the
 * un-batched transport: one ion pass over all primaries, then one neutron pass.
 * A LIVE policy yields several family-complete batches at the configured cadence;
 * scored output then matches the final-only result up to floating-point reduction
 * order (the per-history RNG streams are invariant under range splitting).
 *
 * The checkpoint boundary is where the future parallel/preview machinery attaches
 * — per-worker accumulator merge (#161), a variance batch fold (#169), and an
 * optional dump (#193).  Those are no-ops here; this change lands only the
 * batch-aware seam and the quiescence guarantee.
 */
enum osh_status osh_transport_run_minimal(struct osh_transport_context *transport_ctx,
                                          struct osh_beam_runtime *beam_rt,
                                          struct osh_gemca_runtime const *geom_rt,
                                          struct osh_material_runtime const *material_rt,
                                          struct osh_scoring_runtime *score_rt) {
    struct osh_checkpoint_policy const *policy;
    size_t neutron_capacity;
    size_t nstat;
    size_t done;
    int neutron_enabled; /* cache: neutron pool is present for this run */
    enum osh_status rc;

    if (!transport_ctx) {
        return OSH_EINVAL;
    }
    nstat = transport_ctx->params.nstat;
    if (nstat == 0u) {
        return OSH_EINVAL;
    }
    policy = transport_ctx->checkpoint_policy;

    /* Neutron family: initialise the pool once per run (not per batch, so its
     * cumulative n_created diagnostic and its allocation survive across batches;
     * the neutron pass drains it back to empty at every batch anyway). */
    neutron_enabled = (transport_ctx->neutron_pool != NULL);
    if (neutron_enabled) {
        neutron_capacity = (transport_ctx->params.pool_capacity != 0u) ? transport_ctx->params.pool_capacity
                                                                       : (size_t) OSH_TRANSPORT_POOL_CAPACITY;
        if (transport_ctx->neutron_pool->capacity == 0u) {
            rc = osh_neutron_pool_init(transport_ctx->neutron_pool, neutron_capacity);
            if (rc != OSH_OK) {
                return rc;
            }
        } else {
            osh_neutron_pool_reset(transport_ctx->neutron_pool);
        }
    }

    /* ---- Outer checkpoint batch loop ------------------------------------- */
    done = 0u; /* primaries whose histories have finished; the true completed count */
    while (done < nstat) {
        size_t batch_completed = 0u;
        size_t const remaining = nstat - done;
        /* Measured throughput is not plumbed to the scheduler yet, so a bare
         * time cadence falls back to a single batch (see osh_checkpoint_policy);
         * the count cadence and final-only mode do not need it. */
        size_t const k = osh_checkpoint_next_batch_size(policy, 0.0, remaining);

        rc = run_families_over_range(
            transport_ctx, beam_rt, geom_rt, material_rt, score_rt, done, done + k, neutron_enabled, &batch_completed);
        if (rc != OSH_OK) {
            return rc;
        }

        /* Advance by the primaries that actually finished, never the requested
         * batch size: on a clean stop the batch drains fewer than k, and `done`
         * must stay the exact completed count so the checkpoint boundary (and any
         * future dump/merge/variance hook) never over-reports. */
        done += batch_completed;

        /* CHECKPOINT: the simulation is family-quiescent here.  Future work folds
         * a variance batch (#169), merges per-worker accumulators (#161), and
         * fires an optional dump (#193) at this exact point. */

        /* A short batch means the inner loop hit a clean stop and drained early;
         * do not begin another batch. */
        if (batch_completed < k) {
            break;
        }
    }

    transport_ctx->completed_primaries = done;
    return OSH_OK;
}

/**
 * @brief Run the family scheduler for one checkpoint batch [hist_lo, hist_hi).
 *
 * @details
 * The inner scheduler, factored out so the outer batch loop can drive it once per
 * checkpoint.  It transports the ion primaries of this range and then drains every
 * secondary family (neutrons, …) they banked, so the batch reaches family
 * quiescence before returning.
 *
 * has_work discipline (unchanged from the single-pass driver):
 *   • Clear the family's own has_work BEFORE the dispatch call (avoids an
 *     accidental self-reschedule if the kernel ever pushes to its own pool).
 *   • After the call, set has_work for peer families based on what changed.
 *   • scheduler_next() does NOT clear has_work on return.
 *
 * Scheduling policy for one batch:
 *   1. Ion pass        — processes primaries [hist_lo, hist_hi) from beam_rt and
 *                        reports how many finished via @p batch_completed_out.
 *   2. Neutron pass    — drains neutron_pool if non-empty after the ion pass.
 *   3. Ion feedback    — reschedule ions only if nuclear_neutron_ion_feedback is
 *                        set AND the ion pool received (n,x) secondaries.  Not yet
 *                        reachable: neutron transport does not push to the ion pool.
 *
 * @param[out] batch_completed_out  Primaries in the range whose histories finished
 *                                  (== hist_hi - hist_lo unless a clean stop
 *                                  drained the batch early).
 */
static enum osh_status run_families_over_range(struct osh_transport_context *transport_ctx,
                                               struct osh_beam_runtime *beam_rt,
                                               struct osh_gemca_runtime const *geom_rt,
                                               struct osh_material_runtime const *material_rt,
                                               struct osh_scoring_runtime *score_rt,
                                               size_t hist_lo,
                                               size_t hist_hi,
                                               int neutron_enabled,
                                               size_t *batch_completed_out) {
    struct osh_transport_scheduler scheduler;
    enum osh_transport_family family;
    enum osh_status rc;

    if (batch_completed_out) {
        *batch_completed_out = 0u;
    }

    osh_transport_scheduler_reset(&scheduler);

    rc = osh_transport_scheduler_enable(&scheduler, OSH_TRANSPORT_FAMILY_ION);
    if (rc != OSH_OK) {
        return rc;
    }
    osh_transport_scheduler_set_has_work(&scheduler, OSH_TRANSPORT_FAMILY_ION, 1);

    if (neutron_enabled) {
        rc = osh_transport_scheduler_enable(&scheduler, OSH_TRANSPORT_FAMILY_NEUTRON);
        if (rc != OSH_OK) {
            return rc;
        }
        /* has_work stays 0 until the ion pass produces neutrons. */
    }

    while (osh_transport_scheduler_next(&scheduler, &family)) {
        /* Clear before dispatch so the family does not immediately re-schedule
         * itself; it will only run again if a feedback path sets has_work=1. */
        osh_transport_scheduler_set_has_work(&scheduler, family, 0);

        rc = dispatch_transport_family(
            family, transport_ctx, beam_rt, geom_rt, material_rt, score_rt, hist_lo, hist_hi, batch_completed_out);
        if (rc != OSH_OK) {
            return rc;
        }

        if (family == OSH_TRANSPORT_FAMILY_ION && neutron_enabled) {
            /* Neutrons produced during the ion pass are now in the pool. */
            osh_transport_scheduler_set_has_work(
                &scheduler, OSH_TRANSPORT_FAMILY_NEUTRON, (transport_ctx->neutron_pool->n > 0u) ? 1 : 0);
        }

        if (family == OSH_TRANSPORT_FAMILY_NEUTRON && transport_ctx->params.nuclear_neutron_ion_feedback) {
            /* TODO: replace 0 with (ion_pool->n > 0) once osh_transport_neutron.c
             * is wired to push (n,p)/(n,α)/compound ion secondaries to the ion pool. */
            osh_transport_scheduler_set_has_work(&scheduler, OSH_TRANSPORT_FAMILY_ION, 0);
        }
    }

    return OSH_OK;
}

/**
 * @brief Dispatch one transport family to its implementation module.
 *
 * @details
 * This switch is the future join point between the scheduler and the concrete
 * transport kernels.  Keeping the mapping here avoids spreading knowledge of
 * per-family source files through the rest of transport/.
 *
 * The history range [@p hist_lo, @p hist_hi) selects the primaries for the ion
 * family; the secondary families ignore it and simply drain their pools (whatever
 * the ion pass banked for this batch).  @p ion_completed_out receives the ion
 * family's finished-primary count and is left untouched by the other families.
 *
 * @param[in]     family        Scheduled family to transport.
 * @param[in]     beam_rt       Hot beam runtime for primary generation.
 * @param[in]     geom_rt       Compiled geometry runtime.
 * @param[in]     material_rt   Hot material runtime tables.
 * @param[in,out] score_rt      Scoring runtime.
 * @param[in]     hist_lo       Inclusive lower bound of the ion range.
 * @param[in]     hist_hi       Exclusive upper bound of the ion range.
 * @param[out]    ion_completed_out  Ion primaries finished (ion family only).
 *
 * @returns OSH_OK on success, OSH_ENOTSUP for a family without an
 *          implementation, or another OSH_E* from the family kernel.
 */
static enum osh_status dispatch_transport_family(enum osh_transport_family family,
                                                 struct osh_transport_context *transport_ctx,
                                                 struct osh_beam_runtime *beam_rt,
                                                 struct osh_gemca_runtime const *geom_rt,
                                                 struct osh_material_runtime const *material_rt,
                                                 struct osh_scoring_runtime *score_rt,
                                                 size_t hist_lo,
                                                 size_t hist_hi,
                                                 size_t *ion_completed_out) {
    switch (family) {
    case OSH_TRANSPORT_FAMILY_ION:
        return osh_transport_ion_run_range(
            transport_ctx, beam_rt, geom_rt, material_rt, score_rt, hist_lo, hist_hi, ion_completed_out);
    case OSH_TRANSPORT_FAMILY_NEUTRON:
        return osh_transport_neutron_run(transport_ctx, beam_rt, geom_rt, material_rt, score_rt);
    case OSH_TRANSPORT_FAMILY_PHOTON:
        return osh_transport_photon_run_minimal(transport_ctx, beam_rt, geom_rt, material_rt, score_rt);
    case OSH_TRANSPORT_FAMILY_ELECTRON:
    case OSH_TRANSPORT_FAMILY_COUNT:
        break;
    }

    OSH_DIAG_ERRORF(transport_ctx ? transport_ctx->diag : NULL,
                    "transport: %s transport is not implemented",
                    osh_transport_family_name(family));
    return OSH_ENOTSUP;
}
