#include "transport/osh_transport_ion.h"

#include "beam/runtime/osh_beam_runtime.h"
#include "common/osh_diag.h"
#include "common/osh_particle_pool.h"
#include "common/osh_step_segment.h"
#include "common/osh_time.h"
#include "gemca/runtime/osh_gemca_runtime.h"
#include "material/runtime/osh_material_runtime.h"
#include "physics/atomic/osh_physics_scat_moliere.h"
#include "random/osh_rng.h"
#include "scoring/runtime/osh_scoring_runtime.h"
#include "transport/osh_run_control.h"
#include "transport/osh_transport.h"
#include "transport/osh_transport_ion_step.h"
#include "transport/osh_worker_context.h"

/*
 * Transport pool capacity — number of particle histories alive simultaneously.
 *
 * This is the primary tuning knob between cache efficiency and parallelism:
 *   - Small (e.g. 256): pool fits in L1/L2 cache; minimal working-set pressure.
 *   - Medium (e.g. 4096–65536): pool fits in L2/L3; good CPU SIMD throughput.
 *   - capacity == NSTAT: all primaries live at once; natural for GPU offload.
 *   - capacity == 1: scalar reference — each primary fully transported alone.
 *
 * Changing this constant does not affect physics results, only performance.
 *
 * OSH_TRANSPORT_MAX_STEPS_PER_PRIMARY is a safety guard against particles
 * stuck in degenerate geometry (near-zero steps, persistent boundary nudges,
 * or pathological material configurations).  The wavefront loop multiplies
 * this by nstat to get a total step budget for the run; if exceeded,
 * OSH_ESTATE is returned.  1 000 000 steps is generous for CSDA transport
 * through realistic geometry; a proton Bragg peak typically takes O(1 000)
 * steps with DELTAE = 0.02.
 */
#define OSH_TRANSPORT_MAX_STEPS_PER_PRIMARY 1000000u

#define OSH_TRANSPORT_PROGRESS_MIN_INTERVAL_S 1.0
#define OSH_TRANSPORT_PROGRESS_MAX_INTERVAL_S 10.0
#define OSH_TRANSPORT_PROGRESS_MIN_CHUNK 1000u
#define OSH_TRANSPORT_PROGRESS_TARGET_CHUNKS 20u

/* ---- Forward declarations ------------------------------------------------ */

static size_t transport_progress_chunk_size(size_t total);
static void
report_transport_progress(struct osh_diag_sink const *diag, size_t completed, size_t total, double elapsed_s);
static enum osh_status validate_transport_modes(struct osh_transport_context const *transport_ctx);

/* ---- Wavefront loop ------------------------------------------------------ */

/**
 * @brief Transport every history in a worker's assigned range to termination.
 *
 * @details
 * The wavefront (BFS) CSDA loop, driven by osh_transport_ion_run_range() over an
 * explicit history range rather than always the whole run.  The worker context isolates the *transport-local* mutable
 * state — the particle pool and per-step batch scratch — for one slice [hist_lo, hist_hi). Primary generation is
 * already partition-safe: this loop fills the pool through osh_beam_runtime_fill_pool_at() with an explicit global base
 * derived from the worker's own range, never the shared beam cursor.  Profile counters are per-worker too: this loop
 * accumulates into @c wctx->profile, never a shared profile, so workers do not race on them; the driver folds them into
 * the run's master with osh_transport_profile_merge().  One piece of shared state still stands between this and
 * concurrent execution: @p score_rt's accumulators (deposits land in the shared master); per-worker copies merged at
 * the end are the plan, see the @c accumulators field on @ref osh_worker_context.  Today a single worker covers [0,
 * nstat) and points its profile straight at the master, so none of that sharing is yet exercised.
 *
 *   1. When the pool is empty and primaries remain, fill it from beam_runtime
 *      (up to the worker's pool capacity).
 *   2. Batch-query zone refs and current-medium boundary distances for all live slots.
 *   3. Build a one-segment current-medium step and call osh_transport_ion_step()
 *      for every live slot.  Particles that die (energy cutoff, geometry exit,
 *      blackhole) are marked by zeroing e[slot].
 *   4. Compact the pool, removing dead entries.
 *   5. Repeat until every history in the range is done and the pool is empty.
 *
 * Each primary is seeded from its global history index (rndoffset + hist_lo +
 * the worker-local primary index), assembled by passing hist_lo + primaries_done
 * as the explicit global base to osh_beam_runtime_fill_pool_at().  Splitting the
 * run into arbitrary disjoint sub-ranges and replaying them in any order
 * therefore reproduces the canonical per-history streams exactly — the seed is a
 * pure function of the index, with no shared, mutable beam cursor in the path.
 * Deposits currently land in the shared master accumulators in @p score_rt;
 * per-worker private accumulators (@c wctx->accumulators) will route here once
 * parallel scoring lands.
 */
static enum osh_status run_history_range(struct osh_worker_context *wctx,
                                         struct osh_transport_context *transport_ctx,
                                         struct osh_beam_runtime *beam_rt,
                                         struct osh_gemca_runtime const *geom_rt,
                                         struct osh_material_runtime const *material_rt,
                                         struct osh_scoring_runtime *score_rt,
                                         size_t *completed_out) {
    struct osh_rng_seeding seeding;
    struct osh_transport_params const *params = &transport_ctx->params;
    struct osh_particle_pool *pool = wctx->pool;
    struct osh_zone_ref *zone_refs = wctx->zone_refs;
    double *dist_batch = wctx->dist_batch;
    size_t const nstat = wctx->hist_hi - wctx->hist_lo; /* histories in this range */
    size_t primaries_done;
    size_t primaries_completed;
    size_t n_fill;
    size_t n_wavefront;
    size_t width; /* primaries injected this fill: wavefront cap, headroom left for secondaries */
    size_t i;
    size_t step_budget;
    size_t steps_taken;
    size_t last_report_completed;
    size_t progress_chunk;
    size_t next_report_completed;
    struct osh_step_segment step_seg;
    /* Deposit target for this worker: its private accumulator set + scratch, or
     * NULL fields that osh_score_target_* resolve to the shared master views —
     * the single-worker baseline, bit-for-bit unchanged (issue #230). */
    struct osh_score_target target;
    struct osh_transport_profile *prof = wctx->profile; /* per-worker; never the shared master on the hot path */
    struct osh_run_control const *ctl = transport_ctx->run_control; /* clean-stop / wall-budget policy, or NULL */
    double ctl_t_start;
    double t_start;
    double t_last_report;
    double t_phase = 0.0;
    int stop = 0; /* set once a clean stop is requested: halt injection, keep draining */
    enum osh_status rc = OSH_OK;

    /* Per-history seeding context.  Every primary derives independent BEAM and
     * PHYSICS streams from its global history index, so a history sees the same
     * random draws on any pool capacity, thread, or rank.  (Scored output is
     * therefore invariant up to floating-point summation order in the shared
     * scoring accumulators, not byte-for-byte.)  Keeping BEAM and PHYSICS on
     * separate purposes makes NUCRE/MSCAT/STRAGG comparisons launch the same
     * primary histories.
     *
     * The global history index is rndoffset + (hist_lo + worker-local index).
     * We keep hist_base at rndoffset alone and carry the worker's hist_lo in the
     * explicit global base passed to osh_beam_runtime_fill_pool_at() below, so
     * the full global index is assembled in exactly one place.  Splitting the
     * run into disjoint sub-ranges and running them in any order then reproduces
     * the canonical per-history streams, because each primary's seed depends
     * only on its index — never on a shared, mutable beam cursor. */
    seeding.type = OSH_RNG_TYPE_PCG32;
    seeding.seed = (uint64_t) params->rndseed;
    seeding.hist_base = (uint64_t) params->rndoffset;

    /* Deposit target for the score-step calls below.  The worker owns the pair;
     * NULL fields resolve to the shared master views inside osh_scoring_score_step
     * (osh_score_target_*), so a worker with no private set is bit-for-bit serial. */
    target.acc_set = wctx->accumulators;
    target.scratch = wctx->scratch;

    /* Total step budget: nstat × max_steps, capped at SIZE_MAX to avoid overflow. */
    if (nstat > (size_t) -1 / OSH_TRANSPORT_MAX_STEPS_PER_PRIMARY) {
        step_budget = (size_t) -1;
    } else {
        step_budget = nstat * (size_t) OSH_TRANSPORT_MAX_STEPS_PER_PRIMARY;
    }
    steps_taken = 0u;
    primaries_done = 0u;
    last_report_completed = 0u;
    progress_chunk = transport_progress_chunk_size(nstat);
    next_report_completed = progress_chunk;
    t_start = osh_monotonic_seconds();
    t_last_report = t_start;
    /* Measure the wall budget from the run-level start.  The driver always arms
     * the control with osh_run_control_start() before transport begins, so a
     * non-NULL ctl always carries a valid baseline — including a legitimate 0.0
     * (the monotonic epoch is unspecified).  Fall back to this loop's start only
     * when the run is uncontrolled. */
    ctl_t_start = ctl ? ctl->t_start : t_start;

    report_transport_progress(transport_ctx->diag, 0u, nstat, 0.0);

    /* Loop until every requested primary has finished OR a clean stop was
     * requested — but in both cases keep going while the pool is non-empty, so
     * in-flight histories drain to termination rather than being abandoned.
     * After a stop the pool empties to 0 and the loop exits with an *exact*
     * completed count (primaries_done == primaries actually injected). */
    while ((primaries_done < nstat && !stop) || pool->n > 0u) {
        if (prof) {
            prof->iterations++;
        }

        /* Fill the pool when it is empty and primaries remain — never after a
         * stop request, so no new primary is injected once we begin draining.
         *
         * Inject at most `wavefront_width` primaries, NOT the full pool capacity:
         * the extra capacity is reserved as headroom for the nuclear secondaries
         * this wavefront produces, so a full primary wavefront no longer leaves
         * zero room for its recoils/fragments (issue #213).  A width of 0 (a
         * context built without a simulation) falls back to the full capacity —
         * the pre-#213 behaviour. */
        if (!stop && pool->n == 0u && primaries_done < nstat) {
            width = transport_ctx->ion_wavefront_width;
            if (width == 0u || width > pool->capacity) {
                width = pool->capacity;
            }
            n_fill = nstat - primaries_done;
            if (n_fill > width) {
                n_fill = width;
            }
            if (prof) {
                t_phase = osh_monotonic_seconds();
            }
            /* Cursor-free fill: the global index of the first primary in this
             * chunk is this worker's range base plus the primaries it has
             * already emitted.  No shared beam cursor is read or mutated, so
             * disjoint workers can fill concurrently and remain reproducible. */
            rc = osh_beam_runtime_fill_pool_at(
                beam_rt, &seeding, pool, n_fill, (uint64_t) wctx->hist_lo + (uint64_t) primaries_done);
            if (prof) {
                prof->fill_s += osh_monotonic_seconds() - t_phase;
            }
            if (rc != OSH_OK) {
                return rc;
            }
            primaries_done += n_fill;
        }

        /* Batch geometry: zone-ref lookup (zone + current HU/material) and
         * current-medium boundary distance for all live particles. */
        if (prof) {
            t_phase = osh_monotonic_seconds();
        }
        osh_gemca_runtime_get_zone_ref_batch(
            geom_rt, pool->x, pool->y, pool->z, pool->ux, pool->uy, pool->uz, pool->n, zone_refs);
        if (prof) {
            double const t_now = osh_monotonic_seconds();
            prof->zone_ref_s += t_now - t_phase;
            t_phase = t_now;
        }
        osh_gemca_runtime_get_distance_batch(
            geom_rt, pool->x, pool->y, pool->z, pool->ux, pool->uy, pool->uz, zone_refs, pool->n, dist_batch);
        if (prof) {
            prof->distance_s += osh_monotonic_seconds() - t_phase;
        }

        /* Snapshot wavefront size — secondaries injected this pass are processed
         * on the next iteration, not in the current one. */
        n_wavefront = pool->n;

        /* Advance every live particle by one step */
        if (prof) {
            t_phase = osh_monotonic_seconds();
        }
        for (i = 0u; i < n_wavefront; ++i) {
            if (steps_taken >= step_budget) {
                OSH_DIAG_ERRORF(transport_ctx->diag,
                                "transport: step budget exceeded after %zu steps (pool slot %zu, primaries_done=%zu, "
                                "zone=%zu, e=%.17g, pos=(%.17g, %.17g, %.17g), dir=(%.17g, %.17g, %.17g))",
                                steps_taken,
                                i,
                                primaries_done,
                                zone_refs[i].zone_idx,
                                pool->e[i],
                                pool->x[i],
                                pool->y[i],
                                pool->z[i],
                                pool->ux[i],
                                pool->uy[i],
                                pool->uz[i]);
                rc = OSH_ESTATE;
                return rc;
            }

            step_seg.ds = dist_batch[i];

            /* Each slot draws from its own carried stream, so the random
             * sequence a history sees does not depend on its position in the
             * wavefront or on how many secondaries its siblings spawned. */
            rc = osh_transport_ion_step(pool,
                                        i,
                                        &zone_refs[i],
                                        &step_seg,
                                        1u,
                                        geom_rt,
                                        transport_ctx,
                                        prof,
                                        material_rt,
                                        score_rt,
                                        &target,
                                        &pool->rng[i]);
            if (rc != OSH_OK) {
                OSH_DIAG_ERRORF(transport_ctx->diag,
                                "transport: slot %zu failed with rc=%d zone=%zu boundary_ds=%.17g e=%.17g pos=(%.17g, "
                                "%.17g, %.17g) dir=(%.17g, %.17g, %.17g)",
                                i,
                                (int) rc,
                                zone_refs[i].zone_idx,
                                dist_batch[i],
                                pool->e[i],
                                pool->x[i],
                                pool->y[i],
                                pool->z[i],
                                pool->ux[i],
                                pool->uy[i],
                                pool->uz[i]);
                return rc;
            }
            ++steps_taken;
        }
        if (prof) {
            double const t_now = osh_monotonic_seconds();
            prof->step_s += t_now - t_phase;
            t_phase = t_now;
        }

        /* Compact dead entries (e[i] <= 0) */
        osh_particle_pool_compact(pool);
        if (prof) {
            prof->compact_s += osh_monotonic_seconds() - t_phase;
        }
        primaries_completed = primaries_done - pool->n;
        if (primaries_completed > last_report_completed) {
            double const t_now = osh_monotonic_seconds();
            int const chunk_reached = (primaries_completed >= next_report_completed);
            int const min_interval_elapsed = ((t_now - t_last_report) >= OSH_TRANSPORT_PROGRESS_MIN_INTERVAL_S);
            int const max_interval_elapsed = ((t_now - t_last_report) >= OSH_TRANSPORT_PROGRESS_MAX_INTERVAL_S);
            if (primaries_completed == nstat || (chunk_reached && min_interval_elapsed) || max_interval_elapsed) {
                report_transport_progress(transport_ctx->diag, primaries_completed, nstat, t_now - t_start);
                last_report_completed = primaries_completed;
                t_last_report = t_now;
                while (next_report_completed <= primaries_completed && next_report_completed < nstat) {
                    next_report_completed += progress_chunk;
                }
                if (next_report_completed > nstat) {
                    next_report_completed = nstat;
                }
            }
        }

        /* Safe point: re-evaluate the stop request once per wavefront pass,
         * after compaction so the count reflects only finished histories.  A
         * stop halts *injection* only (see the fill guard above); the pool keeps
         * draining and, on return, the family scheduler still runs the remaining
         * passes over the secondaries these primaries banked — so the partial
         * result is family-exact for the completed primaries (issue #195). */
        if (!stop && ctl) {
            stop = run_ctl_should_stop(ctl, osh_monotonic_seconds() - ctl_t_start, primaries_completed);
        }
    }

    /* The pool is fully drained here, so this is the exact number of primaries
     * whose histories finished — equal to nstat unless a clean stop fired. */
    if (completed_out) {
        *completed_out = primaries_done - pool->n;
    }

    {
        double const t_end = osh_monotonic_seconds();
        double const total_s = t_end - t_start;
        size_t const completed = primaries_done - pool->n;
        double const avg_pps = (total_s > 0.0) ? ((double) completed / total_s) : 0.0;
        unsigned int const th = (unsigned int) (total_s / 3600.0);
        unsigned int const tm = (unsigned int) ((total_s - th * 3600.0) / 60.0);
        unsigned int const ts = (unsigned int) (total_s - th * 3600.0 - tm * 60.0);

        if (prof) {
            /* Accumulate, never assign: the outer batch loop (osh_transport.c) calls
             * this range once per checkpoint batch against the same run master, so
             * wall time and step count must sum across sequential batches exactly as
             * the per-phase timers and iteration counters above already do.  A bare
             * assignment would leave total_s/steps reflecting only the last batch —
             * internally contradictory (e.g. step_s > total_s) and wrong for the
             * benchmark accounting.  For the default single batch the master is zeroed
             * once at profiling-enable, so 0 + x is byte-identical to the old path.
             * This is a distinct axis from osh_transport_profile_merge(), which folds
             * *concurrent* workers and so combines total_s by max, not sum. */
            prof->total_s += total_s;
            prof->steps += (unsigned long long) steps_taken;
        }

        if (last_report_completed < completed) {
            report_transport_progress(transport_ctx->diag, completed, nstat, total_s);
        }

        /* A clean stop is normal operation, not an error — say so explicitly so
         * the true (partial) primary count is visible in the run log. */
        if (stop && completed < nstat) {
            OSH_DIAG_INFOF(transport_ctx->diag,
                           "transport: clean stop after %zu of %zu requested primaries "
                           "(wall budget or stop request); pool drained, output normalised by the true count",
                           completed,
                           nstat);
        }

        if (th > 0u) {
            OSH_DIAG_INFOF(transport_ctx->diag,
                           "Transport complete: %zu primaries in %u:%02u:%02u  (avg %.0f primaries/s)",
                           completed,
                           th,
                           tm,
                           ts,
                           avg_pps);
        } else if (tm > 0u) {
            OSH_DIAG_INFOF(transport_ctx->diag,
                           "Transport complete: %zu primaries in %02u:%02u  (avg %.0f primaries/s)",
                           completed,
                           tm,
                           ts,
                           avg_pps);
        } else {
            OSH_DIAG_INFOF(transport_ctx->diag,
                           "Transport complete: %zu primaries in %u s  (avg %.0f primaries/s)",
                           completed,
                           ts,
                           avg_pps);
        }
    }

    return rc;
}

/* ---- Per-worker profile reduction ---------------------------------------- */

void osh_transport_profile_merge(struct osh_transport_profile *dst, struct osh_transport_profile const *src) {
    if (!dst || !src) {
        return;
    }
    /* Per-phase timers and event counters are additive: their sum across workers
     * is aggregate work-in-phase (which legitimately exceeds elapsed wall time). */
    dst->fill_s += src->fill_s;
    dst->zone_ref_s += src->zone_ref_s;
    dst->distance_s += src->distance_s;
    dst->step_s += src->step_s;
    dst->compact_s += src->compact_s;
    dst->steps += src->steps;
    dst->iterations += src->iterations;
    dst->nuclear_events += src->nuclear_events;
    dst->secondaries += src->secondaries;
    /* Wall time is not additive across concurrent workers — elapsed time is the
     * longest worker's span, so combine by maximum.  Merging into a zeroed dst
     * (the single-worker case) reproduces src->total_s exactly. */
    if (src->total_s > dst->total_s) {
        dst->total_s = src->total_s;
    }
}

/* ---- Public entry: single-worker run over an explicit history range ------ */

/**
 * @brief Transport the ion primaries of one history range [hist_lo, hist_hi).
 *
 * @details
 * Validates the run parameters, borrows the simulation's pre-allocated ion pool
 * and geometry scratch (transport_ctx->ion_pool / zone_refs / dist_batch, sized
 * at simulation-create time and reused across families and batches) through
 * osh_worker_context_attach(), and transports exactly [hist_lo, hist_hi).  The
 * borrowed pool is reset to empty on entry; the wavefront loop drains it back to
 * empty before returning, so consecutive batches never carry live histories
 * across a checkpoint.  Partitioning across several workers (threads/ranks) is a
 * matter of giving each a context over a disjoint sub-range with its own private
 * pool and merging their scoring accumulators; run_history_range() does not
 * change.
 */
enum osh_status osh_transport_ion_run_range(struct osh_transport_context *transport_ctx,
                                            struct osh_beam_runtime *beam_rt,
                                            struct osh_gemca_runtime const *geom_rt,
                                            struct osh_material_runtime const *material_rt,
                                            struct osh_scoring_runtime *score_rt,
                                            struct osh_score_target const *target,
                                            size_t hist_lo,
                                            size_t hist_hi,
                                            size_t *completed_out) {
    struct osh_worker_context wctx;
    struct osh_transport_params const *params;
    enum osh_status rc;

    if (completed_out) {
        *completed_out = 0u;
    }
    if (!transport_ctx || !beam_rt || !geom_rt || !material_rt || !score_rt) {
        return OSH_EINVAL;
    }
    if (!transport_ctx->ion_pool || !transport_ctx->zone_refs || !transport_ctx->dist_batch) {
        return OSH_EINVAL;
    }
    if (hist_hi <= hist_lo) {
        return OSH_EINVAL;
    }
    params = &transport_ctx->params;
    if (params->nstat == 0u) {
        return OSH_EINVAL;
    }
    if (params->deltae <= 0.0f || params->deltae >= 1.0f) {
        return OSH_EINVAL;
    }
    rc = validate_transport_modes(transport_ctx);
    if (rc != OSH_OK) {
        return rc;
    }

    /* Borrow the simulation's pre-allocated ion pool and geometry scratch instead
     * of allocating a private set, and reset the reused pool to empty before
     * transporting this range. */
    transport_ctx->ion_pool->n = 0u;
    osh_worker_context_attach(
        &wctx, hist_lo, hist_hi, transport_ctx->ion_pool, transport_ctx->zone_refs, transport_ctx->dist_batch);

    /* Wire this worker's deposit target: a caller-supplied private set (the replica
     * driver, issue #230) or NULL for the shared-master baseline.  run_history_range
     * builds the score target from these fields; NULL fields resolve to the master
     * so the un-partitioned path is bit-for-bit unchanged. */
    if (target) {
        wctx.accumulators = target->acc_set;
        wctx.naccumulators = target->acc_set ? score_rt->npages : 0u;
        wctx.scratch = target->scratch;
    }

    /* Single worker: point its profile straight at the run's master so counters
     * land there directly — bit-identical to the un-parallelised path, no merge.
     * Batches accumulate into the same master, so per-run totals are unaffected
     * by how the range is sliced.  A parallel driver would instead give each
     * worker a private profile and fold them in with osh_transport_profile_merge(). */
    wctx.profile = transport_ctx->profile;

    return run_history_range(&wctx, transport_ctx, beam_rt, geom_rt, material_rt, score_rt, completed_out);
}

/* ---- Progress helpers ---------------------------------------------------- */

static size_t transport_progress_chunk_size(size_t total) {
    size_t chunk;

    if (total == 0u) {
        return 1u;
    }
    chunk = (total + OSH_TRANSPORT_PROGRESS_TARGET_CHUNKS - 1u) / OSH_TRANSPORT_PROGRESS_TARGET_CHUNKS;
    if (chunk < OSH_TRANSPORT_PROGRESS_MIN_CHUNK && total > OSH_TRANSPORT_PROGRESS_MIN_CHUNK) {
        chunk = OSH_TRANSPORT_PROGRESS_MIN_CHUNK;
    }
    if (chunk == 0u) {
        chunk = 1u;
    }
    return chunk;
}

static void
report_transport_progress(struct osh_diag_sink const *diag, size_t completed, size_t total, double elapsed_s) {
    double primaries_per_second;
    double eta_s;
    size_t remaining;
    int pct;
    unsigned int eta_h;
    unsigned int eta_m;
    unsigned int eta_sec;

    remaining = (completed < total) ? (total - completed) : 0u;
    primaries_per_second = (elapsed_s > 0.0) ? ((double) completed / elapsed_s) : 0.0;
    pct = (total > 0u) ? (int) (completed * 100u / total) : 0;

    if (primaries_per_second > 0.0 && remaining > 0u) {
        eta_s = (double) remaining / primaries_per_second;
        eta_h = (unsigned int) (eta_s / 3600.0);
        eta_m = (unsigned int) ((eta_s - eta_h * 3600.0) / 60.0);
        eta_sec = (unsigned int) (eta_s - eta_h * 3600.0 - eta_m * 60.0);
        if (eta_h > 0u) {
            OSH_DIAG_INFOF(diag,
                           "Progress: %zu/%zu  %3d%%  %.0f primaries/s  ETA %u:%02u:%02u",
                           completed,
                           total,
                           pct,
                           primaries_per_second,
                           eta_h,
                           eta_m,
                           eta_sec);
        } else {
            OSH_DIAG_INFOF(diag,
                           "Progress: %zu/%zu  %3d%%  %.0f primaries/s  ETA %02u:%02u",
                           completed,
                           total,
                           pct,
                           primaries_per_second,
                           eta_m,
                           eta_sec);
        }
    } else {
        OSH_DIAG_INFOF(diag, "Progress: %zu/%zu  %3d%%  %.0f primaries/s", completed, total, pct, primaries_per_second);
    }
}

static enum osh_status validate_transport_modes(struct osh_transport_context const *transport_ctx) {
    switch ((enum osh_transport_mcs_mode) transport_ctx->params.mcs_mode) {
    case OSH_TRANSPORT_MCS_OFF:
    case OSH_TRANSPORT_MCS_GAUSSIAN:
        break;
    case OSH_TRANSPORT_MCS_MOLIERE:
        /* Build the Molière reduced-angle tables once, here on the single
         * setup thread, before the parallel transport loop starts. */
        osh_physics_moliere_init();
        break;
    case OSH_TRANSPORT_MCS_WENTZEL:
        OSH_DIAG_ERRORF(transport_ctx->diag, "%s", "transport: Wentzel-VI/Urban MCS is not implemented");
        return OSH_ENOTSUP;
    default:
        OSH_DIAG_ERRORF(
            transport_ctx->diag, "transport: unsupported MCS mode value %d", (int) transport_ctx->params.mcs_mode);
        return OSH_EINVAL;
    }

    switch ((enum osh_transport_straggling_mode) transport_ctx->params.straggling_mode) {
    case OSH_TRANSPORT_STRAGGLING_OFF:
    case OSH_TRANSPORT_STRAGGLING_GAUSSIAN:
    case OSH_TRANSPORT_STRAGGLING_VAVILOV:
        break;
    case OSH_TRANSPORT_STRAGGLING_URBAN:
        OSH_DIAG_ERRORF(transport_ctx->diag, "%s", "transport: Urban straggling (STRAGG 3) is not yet implemented");
        return OSH_ENOTSUP;
    default:
        OSH_DIAG_ERRORF(transport_ctx->diag,
                        "transport: unsupported straggling mode value %d",
                        (int) transport_ctx->params.straggling_mode);
        return OSH_EINVAL;
    }

    return OSH_OK;
}
