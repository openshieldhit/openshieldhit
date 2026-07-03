#include "transport/osh_transport.h"

#include <stdlib.h>
#include <string.h>

#include "common/osh_diag.h"
#include "common/osh_particle_pool.h"
#include "common/osh_time.h"
#include "scoring/runtime/osh_scoring_accumulator.h"
#include "scoring/runtime/osh_scoring_compile.h"
#include "scoring/runtime/osh_scoring_runtime.h"
#include "scoring/runtime/osh_scoring_snapshot.h"
#include "transport/osh_checkpoint_policy.h"
#include "transport/osh_neutron_pool.h"
#include "transport/osh_run_control.h"
#include "transport/osh_transport_ion.h"
#include "transport/osh_transport_neutron.h"
#include "transport/osh_transport_photon.h"
#include "transport/osh_transport_scheduler.h"

/* ---- Forward declarations ------------------------------------------------ */

static enum osh_status run_master_batched(struct osh_transport_context *transport_ctx,
                                          struct osh_beam_runtime *beam_rt,
                                          struct osh_gemca_runtime const *geom_rt,
                                          struct osh_material_runtime const *material_rt,
                                          struct osh_scoring_runtime *score_rt,
                                          int neutron_enabled,
                                          size_t *completed_out);

static enum osh_status run_score_replicas(struct osh_transport_context *transport_ctx,
                                          struct osh_beam_runtime *beam_rt,
                                          struct osh_gemca_runtime const *geom_rt,
                                          struct osh_material_runtime const *material_rt,
                                          struct osh_scoring_runtime *score_rt,
                                          size_t nreplicas,
                                          int neutron_enabled,
                                          size_t *completed_out);

static enum osh_status run_families_over_range(struct osh_transport_context *transport_ctx,
                                               struct osh_beam_runtime *beam_rt,
                                               struct osh_gemca_runtime const *geom_rt,
                                               struct osh_material_runtime const *material_rt,
                                               struct osh_scoring_runtime *score_rt,
                                               struct osh_score_target const *target,
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
                                                 struct osh_score_target const *target,
                                                 size_t hist_lo,
                                                 size_t hist_hi,
                                                 size_t *ion_completed_out);

/*
 * Orchestrator for the minimal transport run.
 *
 * After the once-per-run family setup (neutron pool, drop counter), the run takes
 * one of two shapes:
 *   • the DEFAULT path (run_master_batched) — the outer checkpoint batch loop that
 *     slices [0, nstat) into family-complete batches depositing into the shared
 *     master accumulators, with the periodic-dump hook (#193).  A NULL / final-only
 *     checkpoint policy makes this one batch of K = nstat: byte-for-byte identical
 *     to the un-batched transport;
 *   • the REPLICA path (run_score_replicas, params.score_replicas >= 1) — the
 *     sequential private-accumulator harness of issue #230: N contiguous
 *     sub-ranges each transported in turn into their own private set, then merged
 *     into the master.  This is the diagnostic that first exercises the
 *     per-worker accumulator/profile merge every parallel backend will rely on,
 *     with zero concurrency because the replicas run one after another.
 *
 * Both share the INNER family scheduler (run_families_over_range): for one range it
 * drains ions then the neutrons/fragments they banked, so the range is
 * *family-exact* before returning.
 */
enum osh_status osh_transport_run_minimal(struct osh_transport_context *transport_ctx,
                                          struct osh_beam_runtime *beam_rt,
                                          struct osh_gemca_runtime const *geom_rt,
                                          struct osh_material_runtime const *material_rt,
                                          struct osh_scoring_runtime *score_rt) {
    size_t neutron_capacity;
    size_t nstat;
    size_t replicas;
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
    replicas = transport_ctx->params.score_replicas;
    /* A partition of [0, nstat) into contiguous integer sub-ranges only tiles the
     * range with no empty part when N <= nstat; reject the misuse up front with a
     * clear diagnostic rather than launching an empty range (issue #230). */
    if (replicas > nstat) {
        OSH_DIAG_ERRORF(transport_ctx->diag,
                        "transport: --score-replicas %zu exceeds nstat %zu (each replica needs at least one history)",
                        replicas,
                        nstat);
        return OSH_EINVAL;
    }

    /* Neutron family: initialise the pool once per run (not per batch/replica, so
     * its cumulative n_created diagnostic and its allocation survive across ranges;
     * the neutron pass drains it back to empty at every range anyway). */
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

    /* Reset the ion-secondary drop counter for this run (mirrors the neutron
     * pool reset above), so a re-run reports drops for the current run only.
     * It then accumulates across all ranges, since the per-range
     * osh_transport_ion_run_range resets only the pool's live count, never
     * n_dropped. */
    if (transport_ctx->ion_pool != NULL) {
        transport_ctx->ion_pool->n_dropped = 0u;
    }

    done = 0u;
    if (replicas == 0u) {
        /* Default: deposit into the shared master with the checkpoint/dump loop. */
        rc = run_master_batched(transport_ctx, beam_rt, geom_rt, material_rt, score_rt, neutron_enabled, &done);
    } else {
        /* Diagnostic: N sequential private-accumulator replicas, merged into the
         * master before the caller's postprocess + save (issue #230). */
        rc = run_score_replicas(
            transport_ctx, beam_rt, geom_rt, material_rt, score_rt, replicas, neutron_enabled, &done);
    }
    if (rc != OSH_OK) {
        return rc;
    }

    transport_ctx->completed_primaries = done;
    return OSH_OK;
}

/*
 * Default path: the outer checkpoint batch loop over [0, nstat), depositing into
 * the shared master accumulators (target == NULL) and firing the periodic-dump
 * hook at each family-complete boundary.  Factored out of osh_transport_run_minimal
 * so the replica path is a sibling rather than a special case buried in the loop.
 *
 * The batch size K comes from the checkpoint policy.  FINAL-ONLY (a NULL policy or
 * OSH_PARTIAL_NONE) yields a single batch of K = nstat, so the loop runs exactly
 * once and this is byte-for-byte identical to the un-batched transport.  A LIVE
 * policy yields several family-complete batches at the configured cadence; scored
 * output then matches the final-only result up to floating-point reduction order.
 */
static enum osh_status run_master_batched(struct osh_transport_context *transport_ctx,
                                          struct osh_beam_runtime *beam_rt,
                                          struct osh_gemca_runtime const *geom_rt,
                                          struct osh_material_runtime const *material_rt,
                                          struct osh_scoring_runtime *score_rt,
                                          int neutron_enabled,
                                          size_t *completed_out) {
    struct osh_checkpoint_policy const *policy = transport_ctx->checkpoint_policy;
    struct osh_run_control *ctl = transport_ctx->run_control; /* clean-stop / dump policy, or NULL */
    size_t const nstat = transport_ctx->params.nstat;
    size_t done;          /* primaries whose histories have finished; the true completed count */
    double measured_rate; /* primaries/s from the last batch; seeds the adaptive time cadence */
    enum osh_status rc;

    done = 0u;
    measured_rate = 0.0;
    while (done < nstat) {
        size_t batch_completed = 0u;
        size_t const remaining = nstat - done;
        /* Size the batch from the policy.  The measured rate feeds only the
         * adaptive time cadence; count cadence, explicit batch, and final-only
         * ignore it (final-only returns the whole remainder → one pass, unchanged). */
        size_t const k = osh_checkpoint_next_batch_size(policy, measured_rate, remaining);
        double t_batch_start = osh_monotonic_seconds();
        double batch_s;
        int dump_destination_ready; /* set below: the driver wired a sink + shadow for dumps */

        /* target == NULL: deposit straight into the shared master views. */
        rc = run_families_over_range(transport_ctx,
                                     beam_rt,
                                     geom_rt,
                                     material_rt,
                                     score_rt,
                                     NULL,
                                     done,
                                     done + k,
                                     neutron_enabled,
                                     &batch_completed);
        if (rc != OSH_OK) {
            return rc;
        }

        /* Advance by the primaries that actually finished, never the requested
         * batch size: on a clean stop the batch drains fewer than k, and `done`
         * must stay the exact completed count so the checkpoint boundary (and any
         * dump/merge/variance hook) never over-reports. */
        done += batch_completed;

        /* Refresh the throughput estimate for the next batch's time-cadence
         * sizing.  Guard against a zero-length interval (coarse clock / tiny
         * batch): keep the previous estimate rather than dividing by ~0. */
        batch_s = osh_monotonic_seconds() - t_batch_start;
        if (batch_s > 0.0 && batch_completed > 0u) {
            measured_rate = (double) batch_completed / batch_s;
        }

        /* CHECKPOINT: the simulation is family-quiescent here, so a dump taken now
         * is physically EXACT (every secondary family the completed primaries
         * banked has been drained into scoring — issue #195).  Fire a periodic /
         * on-demand snapshot when one is due, but never at the final boundary
         * (done == nstat): the run's own end-of-run save already writes the
         * complete result, so a dump there would be redundant work. */
        /* A dump needs a wired destination (sink + reusable shadow scratch); when
         * the driver configured no dumps these stay NULL and the block is skipped.
         * Named separately from the "not the final boundary" guard so the intent of
         * each half is clear. */
        dump_destination_ready = (ctl != NULL) && (ctl->dump_sink != NULL) && (ctl->dump_shadow != NULL);
        if (done < nstat && dump_destination_ready) {
            /* Measure elapsed from the run-control baseline (ctl->t_start, set by
             * osh_run_control_start before transport), the same clock the wall-budget
             * stop path uses — so the dump time cadence and the stop cadence never
             * drift relative to each other or to a caller reading now - ctl->t_start. */
            double const elapsed = osh_monotonic_seconds() - ctl->t_start;
            if (run_ctl_should_dump(ctl, elapsed, done)) {
                enum osh_status const drc = osh_scoring_snapshot_save(
                    ctl->dump_sink, ctl->dump_shadow, (unsigned long long) done, ctl->dump_outputs, ctl->dump_noutputs);
                if (drc != OSH_OK) {
                    /* Fail-soft: a dump is a preview, never the run's product, so a
                     * write/alloc failure warns and the run keeps going to its
                     * exact final save rather than aborting. */
                    OSH_DIAG_WARNF(transport_ctx->diag,
                                   "transport: partial-result dump failed (rc=%d) after %zu primaries; run continues",
                                   (int) drc,
                                   done);
                }
            }
        }

        /* A short batch means the inner loop hit a clean stop and drained early;
         * do not begin another batch. */
        if (batch_completed < k) {
            break;
        }
    }

    *completed_out = done;
    return OSH_OK;
}

/*
 * Split [0, nstat) into @p nreplicas contiguous integer sub-ranges, transport each
 * one sequentially into its OWN private accumulator set, then merge every set into
 * the shared master before the caller's postprocess + save (issue #230).  This is
 * the first code in the tree that exercises the private-accumulator +
 * osh_scoring_accumulator_merge() reduce that threads/MPI/WASM will all depend on —
 * with zero concurrency risk, because the replicas run one after another.
 *
 * Reproducibility contract (the deterministic core of issue #168):
 *   - N == 1        → one range into a private set merged into an empty master is
 *                     the same summation order as serial, so bit-identical to it.
 *   - N  > 1        → same per-history physics (RNG seeds are a pure function of the
 *                     global index); only cross-partition FP summation order differs,
 *                     so it matches serial within compare_dat.py tolerance.
 *   - fixed N twice → identical partition ⇒ deterministic merge order ⇒ bit-identical
 *                     run to run (the guard against leftover shared-mutable state).
 *   - N vs M (N!=M) → different partition ⇒ different merge grouping; NOT guaranteed
 *                     identical (documented, not asserted).
 *
 * DEVELOPER.md §10: every allocation (private sets + scratch) happens here at setup
 * and every merge at the quiescent boundary after each range fully drains — never on
 * the per-step hot path.
 */
static enum osh_status run_score_replicas(struct osh_transport_context *transport_ctx,
                                          struct osh_beam_runtime *beam_rt,
                                          struct osh_gemca_runtime const *geom_rt,
                                          struct osh_material_runtime const *material_rt,
                                          struct osh_scoring_runtime *score_rt,
                                          size_t nreplicas,
                                          int neutron_enabled,
                                          size_t *completed_out) {
    /* Accumulators for all replicas live in one flat block indexed [r*npages + p]:
     * replica r's deposit target is the slice &private_acc[r*npages].  A flat block
     * (vs. an array of per-replica pointers) is one allocation and one free. */
    struct osh_scoring_accumulator *private_acc;  /* flat nreplicas*npages block, or NULL when npages==0 */
    struct osh_scoring_scratch *private_scratch;  /* one traversal scratch per replica */
    struct osh_scoring_accumulator *master;       /* master accumulator view (merge destination) */
    struct osh_transport_profile *master_profile; /* run master profile, or NULL when profiling off */
    size_t const nstat = transport_ctx->params.nstat;
    size_t const npages = score_rt->npages;
    size_t completed_total;
    size_t r;
    size_t p;
    enum osh_status rc;

    *completed_out = 0u;

    /* ---- Setup: one private accumulator set + scratch per replica ---------- */
    /* Guard the flat-block element count against size_t overflow (npages is small,
     * nreplicas <= nstat, so this never trips in practice — it is a correctness
     * backstop, not an expected path). */
    if (npages > 0u && nreplicas > (size_t) -1 / npages) {
        return OSH_ENOMEM;
    }
    private_acc = (npages > 0u) ? (struct osh_scoring_accumulator *) calloc(nreplicas * npages, sizeof(*private_acc))
                                : NULL;
    private_scratch = (struct osh_scoring_scratch *) calloc(nreplicas, sizeof(*private_scratch));
    if ((npages > 0u && !private_acc) || !private_scratch) {
        free(private_acc);
        free(private_scratch);
        return OSH_ENOMEM;
    }
    rc = OSH_OK;
    for (r = 0u; r < nreplicas; ++r) {
        rc = osh_scoring_runtime_alloc_accumulator_set(score_rt, &private_acc[r * npages]);
        if (rc != OSH_OK) {
            break;
        }
        rc = osh_scoring_runtime_clone_scratch(score_rt, &private_scratch[r]);
        if (rc != OSH_OK) {
            break;
        }
    }
    if (rc != OSH_OK) {
        OSH_DIAG_ERRORF(transport_ctx->diag, "%s", "transport: failed to allocate score-replica private sets");
        goto cleanup;
    }

    /* ---- Run each contiguous sub-range into its own private set, in turn ---- */
    master_profile = transport_ctx->profile;
    completed_total = 0u;
    for (r = 0u; r < nreplicas; ++r) {
        /* Contiguous integer partition: hi(N-1) == nstat exactly, and consecutive
         * ranges share no history, so the union is [0, nstat) with no gap/overlap
         * for any N <= nstat.  64-bit intermediate avoids nstat*r overflow. */
        size_t const lo = (size_t) ((unsigned long long) nstat * r / nreplicas);
        size_t const hi = (size_t) ((unsigned long long) nstat * (r + 1u) / nreplicas);
        size_t completed = 0u;
        struct osh_score_target target;
        struct osh_transport_profile replica_profile;

        target.acc_set = &private_acc[r * npages];
        target.scratch = &private_scratch[r];

        /* Give the replica its own zeroed profile so the per-worker profile reduce
         * (osh_transport_profile_merge) is exercised too — the ion range reads
         * transport_ctx->profile, so repoint it for the duration of this replica.
         * Skipped entirely when profiling is off (master_profile == NULL). */
        if (master_profile) {
            memset(&replica_profile, 0, sizeof(replica_profile));
            transport_ctx->profile = &replica_profile;
        }

        rc = run_families_over_range(
            transport_ctx, beam_rt, geom_rt, material_rt, score_rt, &target, lo, hi, neutron_enabled, &completed);

        if (master_profile) {
            osh_transport_profile_merge(master_profile, &replica_profile);
            transport_ctx->profile = master_profile; /* restore before the next replica / return */
        }
        if (rc != OSH_OK) {
            goto cleanup;
        }

        completed_total += completed;
        /* A short range means a clean stop drained it early; do not launch more
         * replicas (their un-run private sets stay zeroed and merge as identity). */
        if (completed < hi - lo) {
            break;
        }
    }

    /* ---- Merge boundary: fold every private set into the shared master ------ */
    /* Quiescent point — every range has fully drained its families — so this is
     * the exact seam issue #161 reserves.  Merging a private set into the empty
     * master is additive (data += data), so N == 1 reproduces serial bit-for-bit;
     * a zeroed (un-run) set folds in as the identity. */
    master = osh_scoring_runtime_master_accumulators(score_rt);
    for (r = 0u; r < nreplicas; ++r) {
        for (p = 0u; p < npages; ++p) {
            rc = osh_scoring_accumulator_merge(&master[p], &private_acc[r * npages + p]);
            if (rc != OSH_OK) {
                OSH_DIAG_ERRORF(
                    transport_ctx->diag, "transport: score-replica merge failed (rc=%d) on page %zu", (int) rc, p);
                goto cleanup;
            }
        }
    }

    *completed_out = completed_total;
    rc = OSH_OK;

cleanup:
    for (r = 0u; r < nreplicas; ++r) {
        if (private_acc) {
            osh_scoring_runtime_free_accumulator_set(&private_acc[r * npages], npages);
        }
        osh_scoring_runtime_free_scratch(&private_scratch[r]);
    }
    free(private_acc);
    free(private_scratch);
    return rc;
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
 * @param[in]  target               Deposit target threaded to every family's
 *                                  scoring; NULL routes to the shared master views
 *                                  (the default path), a private set routes to a
 *                                  replica's own accumulators (issue #230).
 * @param[out] batch_completed_out  Primaries in the range whose histories finished
 *                                  (== hist_hi - hist_lo unless a clean stop
 *                                  drained the batch early).
 */
static enum osh_status run_families_over_range(struct osh_transport_context *transport_ctx,
                                               struct osh_beam_runtime *beam_rt,
                                               struct osh_gemca_runtime const *geom_rt,
                                               struct osh_material_runtime const *material_rt,
                                               struct osh_scoring_runtime *score_rt,
                                               struct osh_score_target const *target,
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

        rc = dispatch_transport_family(family,
                                       transport_ctx,
                                       beam_rt,
                                       geom_rt,
                                       material_rt,
                                       score_rt,
                                       target,
                                       hist_lo,
                                       hist_hi,
                                       batch_completed_out);
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
 * @param[in]     target        Deposit target for this family's scoring (NULL ⇒
 *                              shared master views).
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
                                                 struct osh_score_target const *target,
                                                 size_t hist_lo,
                                                 size_t hist_hi,
                                                 size_t *ion_completed_out) {
    switch (family) {
    case OSH_TRANSPORT_FAMILY_ION:
        return osh_transport_ion_run_range(
            transport_ctx, beam_rt, geom_rt, material_rt, score_rt, target, hist_lo, hist_hi, ion_completed_out);
    case OSH_TRANSPORT_FAMILY_NEUTRON:
        return osh_transport_neutron_run(transport_ctx, beam_rt, geom_rt, material_rt, score_rt, target);
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
