#include "transport/osh_transport_ion.h"

#include <stdlib.h>

#include "beam/runtime/osh_beam_runtime.h"
#include "common/osh_diag.h"
#include "common/osh_particle_pool.h"
#include "common/osh_step_segment.h"
#include "common/osh_time.h"
#include "gemca/runtime/osh_gemca_runtime.h"
#include "material/runtime/osh_material_runtime.h"
#include "random/osh_rng.h"
#include "transport/osh_transport.h"
#include "transport/osh_transport_ion_step.h"

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
/* Overridable from the build line (e.g. -DOSH_TRANSPORT_POOL_CAPACITY=256)
 * so the benchmark harness can sweep capacities without editing sources. */
#ifndef OSH_TRANSPORT_POOL_CAPACITY
#define OSH_TRANSPORT_POOL_CAPACITY 4096u
#endif
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
 * @brief Wavefront (BFS) CSDA transport loop for ion primaries.
 *
 * @details
 * The outer loop runs until all params->nstat primaries have been generated and
 * transported to termination:
 *
 *   1. When the pool is empty and primaries remain, fill it from beam_runtime
 *      (up to OSH_TRANSPORT_POOL_CAPACITY primaries).
 *   2. Batch-query zone refs and current-medium boundary distances for all live slots.
 *   3. Build a one-segment current-medium step and call osh_transport_ion_step()
 *      for every live slot.  Particles that
 *      die (energy cutoff, geometry exit, blackhole) are marked by
 *      zeroing e[slot].
 *   4. Compact the pool, removing dead entries.
 *   5. Repeat until primaries_done == params->nstat and pool is empty.
 *
 * The pool capacity controls the trade-off between cache pressure (small) and
 * parallelism (large).  Physics results are independent of capacity.
 */
enum osh_status osh_transport_ion_run_minimal(struct osh_transport_context *transport_ctx,
                                              struct osh_beam_runtime *beam_rt,
                                              struct osh_gemca_runtime const *geom_rt,
                                              struct osh_material_runtime const *material_rt,
                                              struct osh_scoring_runtime *score_rt) {
    struct osh_rng_seeding seeding;
    struct osh_transport_params const *params;
    struct osh_particle_pool *pool = NULL;
    struct osh_zone_ref *zone_refs = NULL;
    double *dist_batch = NULL;
    size_t capacity;
    size_t primaries_done;
    size_t primaries_completed;
    size_t n_fill;
    size_t n_wavefront;
    size_t i;
    size_t step_budget;
    size_t steps_taken;
    size_t last_report_completed;
    size_t progress_chunk;
    size_t next_report_completed;
    struct osh_step_segment step_seg;
    struct osh_transport_profile *prof;
    double t_start;
    double t_last_report;
    double t_phase = 0.0;
    enum osh_status rc = OSH_OK;

    if (!transport_ctx || !beam_rt || !geom_rt || !material_rt || !score_rt) {
        return OSH_EINVAL;
    }
    prof = transport_ctx->profile;
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

    /* Pool capacity: runtime override (params->pool_capacity) when set, else the
     * compiled default.  Never exceeds nstat (no point allocating more slots
     * than primaries).  Capacity is a pure performance knob: per-history RNG
     * streams make each history consume the same random draws at any capacity,
     * so scored results match up to floating-point summation order in the
     * shared scoring accumulators. */
    capacity = (params->pool_capacity != 0u) ? params->pool_capacity : (size_t) OSH_TRANSPORT_POOL_CAPACITY;
    if (capacity > params->nstat) {
        capacity = params->nstat;
    }
    rc = osh_particle_pool_alloc(capacity, &pool);
    if (rc != OSH_OK) {
        return rc;
    }

    /* Per-run batch scratch for the geometry queries, sized to the pool.
     * Allocated once here (never on the per-step hot path) and freed at
     * cleanup; replaces fixed-size stack arrays so capacity can be chosen at
     * runtime without a recompile.  Guard the size product against size_t
     * overflow now that capacity is caller-controlled. */
    if (capacity > ((size_t) -1) / sizeof(*zone_refs)) {
        rc = OSH_ENOMEM;
        goto cleanup;
    }
    zone_refs = (struct osh_zone_ref *) malloc(capacity * sizeof(*zone_refs));
    dist_batch = (double *) malloc(capacity * sizeof(*dist_batch));
    if (!zone_refs || !dist_batch) {
        rc = OSH_ENOMEM;
        goto cleanup;
    }

    /* Per-history seeding context.  Every primary derives independent BEAM and
     * PHYSICS streams from its global history index (rndoffset + prim_idx), so a
     * history sees the same random draws on any pool capacity, thread, or rank.
     * (Scored output is therefore invariant up to floating-point summation order
     * in the shared scoring accumulators, not byte-for-byte.)  rndoffset is the
     * global history-index base: disjoint ranges (e.g. one per MPI rank) yield
     * disjoint, non-overlapping streams.  Keeping BEAM and PHYSICS on separate
     * purposes makes NUCRE/MSCAT/STRAGG comparisons launch the same primary
     * histories. */
    seeding.type = OSH_RNG_TYPE_PCG32;
    seeding.seed = (uint64_t) params->rndseed;
    seeding.hist_base = (uint64_t) params->rndoffset;

    /* Total step budget: nstat × max_steps, capped at SIZE_MAX to avoid overflow. */
    if (params->nstat > (size_t) -1 / OSH_TRANSPORT_MAX_STEPS_PER_PRIMARY) {
        step_budget = (size_t) -1;
    } else {
        step_budget = params->nstat * (size_t) OSH_TRANSPORT_MAX_STEPS_PER_PRIMARY;
    }
    steps_taken = 0u;
    primaries_done = 0u;
    last_report_completed = 0u;
    progress_chunk = transport_progress_chunk_size(params->nstat);
    next_report_completed = progress_chunk;
    t_start = osh_monotonic_seconds();
    t_last_report = t_start;

    report_transport_progress(transport_ctx->diag, 0u, params->nstat, 0.0);

    while (primaries_done < params->nstat || pool->n > 0u) {
        if (prof) {
            prof->iterations++;
        }

        /* Fill the pool when it is empty and primaries remain */
        if (pool->n == 0u && primaries_done < params->nstat) {
            n_fill = params->nstat - primaries_done;
            if (n_fill > pool->capacity) {
                n_fill = pool->capacity;
            }
            if (prof) {
                t_phase = osh_monotonic_seconds();
            }
            rc = osh_beam_runtime_fill_pool(beam_rt, &seeding, pool, n_fill);
            if (prof) {
                prof->fill_s += osh_monotonic_seconds() - t_phase;
            }
            if (rc != OSH_OK) {
                goto cleanup;
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
                goto cleanup;
            }

            step_seg.ds = dist_batch[i];

            /* Each slot draws from its own carried stream, so the random
             * sequence a history sees does not depend on its position in the
             * wavefront or on how many secondaries its siblings spawned. */
            rc = osh_transport_ion_step(
                pool, i, &zone_refs[i], &step_seg, 1u, geom_rt, transport_ctx, material_rt, score_rt, &pool->rng[i]);
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
                goto cleanup;
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
            if (primaries_completed == params->nstat || (chunk_reached && min_interval_elapsed)
                || max_interval_elapsed) {
                report_transport_progress(transport_ctx->diag, primaries_completed, params->nstat, t_now - t_start);
                last_report_completed = primaries_completed;
                t_last_report = t_now;
                while (next_report_completed <= primaries_completed && next_report_completed < params->nstat) {
                    next_report_completed += progress_chunk;
                }
                if (next_report_completed > params->nstat) {
                    next_report_completed = params->nstat;
                }
            }
        }
    }

    {
        double const t_end = osh_monotonic_seconds();
        double const total_s = t_end - t_start;
        double const avg_pps = (total_s > 0.0) ? ((double) params->nstat / total_s) : 0.0;
        unsigned int const th = (unsigned int) (total_s / 3600.0);
        unsigned int const tm = (unsigned int) ((total_s - th * 3600.0) / 60.0);
        unsigned int const ts = (unsigned int) (total_s - th * 3600.0 - tm * 60.0);

        if (prof) {
            prof->total_s = total_s;
            prof->steps = (unsigned long long) steps_taken;
        }

        if (last_report_completed < params->nstat) {
            report_transport_progress(transport_ctx->diag, params->nstat, params->nstat, total_s);
        }

        if (th > 0u) {
            OSH_DIAG_INFOF(transport_ctx->diag,
                           "Transport complete: %zu primaries in %u:%02u:%02u  (avg %.0f primaries/s)",
                           params->nstat,
                           th,
                           tm,
                           ts,
                           avg_pps);
        } else if (tm > 0u) {
            OSH_DIAG_INFOF(transport_ctx->diag,
                           "Transport complete: %zu primaries in %02u:%02u  (avg %.0f primaries/s)",
                           params->nstat,
                           tm,
                           ts,
                           avg_pps);
        } else {
            OSH_DIAG_INFOF(transport_ctx->diag,
                           "Transport complete: %zu primaries in %u s  (avg %.0f primaries/s)",
                           params->nstat,
                           ts,
                           avg_pps);
        }
    }

cleanup:
    free(zone_refs);
    free(dist_batch);
    osh_particle_pool_free(pool);
    return rc;
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
    case OSH_TRANSPORT_MCS_MOLIERE:
        break;
    case OSH_TRANSPORT_MCS_GAUSSIAN:
        OSH_DIAG_ERRORF(transport_ctx->diag, "%s", "transport: Gaussian MCS is not implemented");
        return OSH_ENOTSUP;
    default:
        OSH_DIAG_ERRORF(
            transport_ctx->diag, "transport: unsupported MCS mode value %d", (int) transport_ctx->params.mcs_mode);
        return OSH_EINVAL;
    }

    switch ((enum osh_transport_straggling_mode) transport_ctx->params.straggling_mode) {
    case OSH_TRANSPORT_STRAGGLING_OFF:
    case OSH_TRANSPORT_STRAGGLING_GAUSSIAN:
        break;
    case OSH_TRANSPORT_STRAGGLING_VAVILOV:
        OSH_DIAG_ERRORF(transport_ctx->diag, "%s", "transport: Vavilov straggling is not implemented");
        return OSH_ENOTSUP;
    default:
        OSH_DIAG_ERRORF(transport_ctx->diag,
                        "transport: unsupported straggling mode value %d",
                        (int) transport_ctx->params.straggling_mode);
        return OSH_EINVAL;
    }

    return OSH_OK;
}
