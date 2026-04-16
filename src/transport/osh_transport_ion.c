#include "transport/osh_transport_ion.h"

#include <time.h>

#include "beam/osh_beam.h"
#include "beam/runtime/osh_beam_runtime.h"
#include "common/osh_logger.h"
#include "common/osh_particle_pool.h"
#include "gemca/runtime/osh_gemca_runtime.h"
#include "random/osh_rng.h"
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
#define OSH_TRANSPORT_POOL_CAPACITY 4096u
#define OSH_TRANSPORT_MAX_STEPS_PER_PRIMARY 1000000u

#define OSH_TRANSPORT_PROGRESS_MIN_INTERVAL_S 1.0
#define OSH_TRANSPORT_PROGRESS_MAX_INTERVAL_S 10.0
#define OSH_TRANSPORT_PROGRESS_MIN_CHUNK 1000u
#define OSH_TRANSPORT_PROGRESS_TARGET_CHUNKS 20u

/* ---- Forward declarations ------------------------------------------------ */

static double monotonic_seconds(void);
static size_t transport_progress_chunk_size(size_t total);
static void report_transport_progress(size_t completed, size_t total, double elapsed_s);

/* ---- Wavefront loop ------------------------------------------------------ */

/**
 * @brief Wavefront (BFS) CSDA transport loop for ion primaries.
 *
 * @details
 * The outer loop runs until all beam->nstat primaries have been generated and
 * transported to termination:
 *
 *   1. When the pool is empty and primaries remain, fill it from beam_runtime
 *      (up to OSH_TRANSPORT_POOL_CAPACITY primaries).
 *   2. Batch-query zone indices and boundary distances for all live slots.
 *   3. Call osh_transport_ion_step_one() for every live slot.  Particles that
 *      die (energy cutoff, geometry exit, blackhole) are marked by
 *      zeroing e[slot].
 *   4. Compact the pool, removing dead entries.
 *   5. Repeat until primaries_done == beam->nstat and pool is empty.
 *
 * The pool capacity controls the trade-off between cache pressure (small) and
 * parallelism (large).  Physics results are independent of capacity.
 */
enum osh_status osh_transport_ion_run_minimal(struct beam_workspace const *beam,
                                              struct gemca_runtime const *geom_rt,
                                              struct osh_material_runtime const *tables,
                                              struct osh_scoring_runtime *scoring) {
    struct osh_rng rng;
    struct osh_beam_runtime *beam_rt = NULL;
    struct osh_particle_pool *pool = NULL;
    size_t zone_batch[OSH_TRANSPORT_POOL_CAPACITY]; /* stack-safe: compile-time constant */
    double dist_batch[OSH_TRANSPORT_POOL_CAPACITY]; /* ~32 KiB at 4096 doubles */
    size_t capacity;
    size_t primaries_done;
    size_t primaries_completed;
    size_t n_fill;
    size_t i;
    size_t step_budget;
    size_t steps_taken;
    size_t last_report_completed;
    size_t progress_chunk;
    size_t next_report_completed;
    double t_start;
    double t_last_report;
    enum osh_status rc = OSH_OK;

    if (!beam || !geom_rt || !tables || !scoring) {
        return OSH_EINVAL;
    }
    if (beam->nstat == 0u) {
        return OSH_EINVAL;
    }
    if (beam->deltae <= 0.0f || beam->deltae >= 1.0f) {
        return OSH_EINVAL;
    }

    rc = osh_beam_runtime_setup(beam, &beam_rt);
    if (rc != OSH_OK) {
        return rc;
    }

    capacity =
        (beam->nstat < (size_t) OSH_TRANSPORT_POOL_CAPACITY) ? beam->nstat : (size_t) OSH_TRANSPORT_POOL_CAPACITY;
    rc = osh_particle_pool_alloc(capacity, &pool);
    if (rc != OSH_OK) {
        osh_beam_runtime_free(beam_rt);
        return rc;
    }

    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, (uint64_t) beam->rndseed, (uint64_t) beam->rndoffset);

    /* Total step budget: nstat × max_steps, capped at SIZE_MAX to avoid overflow. */
    if (beam->nstat > (size_t) -1 / OSH_TRANSPORT_MAX_STEPS_PER_PRIMARY) {
        step_budget = (size_t) -1;
    } else {
        step_budget = beam->nstat * (size_t) OSH_TRANSPORT_MAX_STEPS_PER_PRIMARY;
    }
    steps_taken = 0u;
    primaries_done = 0u;
    primaries_completed = 0u;
    last_report_completed = 0u;
    progress_chunk = transport_progress_chunk_size(beam->nstat);
    next_report_completed = progress_chunk;
    t_start = monotonic_seconds();
    t_last_report = t_start;

    report_transport_progress(0u, beam->nstat, 0.0);

    while (primaries_done < beam->nstat || pool->n > 0u) {
        /* Fill the pool when it is empty and primaries remain */
        if (pool->n == 0u && primaries_done < beam->nstat) {
            n_fill = beam->nstat - primaries_done;
            if (n_fill > pool->capacity) {
                n_fill = pool->capacity;
            }
            rc = osh_beam_runtime_fill_pool(beam_rt, &rng, pool, n_fill);
            if (rc != OSH_OK) {
                goto cleanup;
            }
            primaries_done += n_fill;
        }

        /* Batch geometry: zone lookup and boundary distance for all live particles */
        osh_gemca_runtime_get_zone_batch(
            geom_rt, pool->x, pool->y, pool->z, pool->ux, pool->uy, pool->uz, pool->n, zone_batch);
        osh_gemca_runtime_get_distance_batch(
            geom_rt, pool->x, pool->y, pool->z, pool->ux, pool->uy, pool->uz, zone_batch, pool->n, dist_batch);

        /* Advance every live particle by one step */
        for (i = 0u; i < pool->n; ++i) {
            if (steps_taken >= step_budget) {
                osh_error("transport: step budget exceeded after %zu steps (pool slot %zu, primaries_done=%zu, "
                          "zone=%zu, e=%.17g, pos=(%.17g, %.17g, %.17g), dir=(%.17g, %.17g, %.17g))",
                          steps_taken,
                          i,
                          primaries_done,
                          zone_batch[i],
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
            rc = osh_transport_ion_step_one(
                pool, i, zone_batch[i], dist_batch[i], geom_rt, beam, tables, scoring, (double) beam->deltae, &rng);
            if (rc != OSH_OK) {
                osh_error("transport: slot %zu failed with rc=%d zone=%zu boundary_ds=%.17g e=%.17g pos=(%.17g, %.17g, "
                          "%.17g) dir=(%.17g, %.17g, %.17g)",
                          i,
                          (int) rc,
                          zone_batch[i],
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

        /* Compact dead entries (e[i] <= 0) */
        osh_particle_pool_compact(pool);
        primaries_completed = primaries_done - pool->n;
        if (primaries_completed > last_report_completed) {
            double const t_now = monotonic_seconds();
            int const chunk_reached = (primaries_completed >= next_report_completed);
            int const min_interval_elapsed = ((t_now - t_last_report) >= OSH_TRANSPORT_PROGRESS_MIN_INTERVAL_S);
            int const max_interval_elapsed = ((t_now - t_last_report) >= OSH_TRANSPORT_PROGRESS_MAX_INTERVAL_S);
            if (primaries_completed == beam->nstat || (chunk_reached && min_interval_elapsed) || max_interval_elapsed) {
                report_transport_progress(primaries_completed, beam->nstat, t_now - t_start);
                last_report_completed = primaries_completed;
                t_last_report = t_now;
                while (next_report_completed <= primaries_completed && next_report_completed < beam->nstat) {
                    next_report_completed += progress_chunk;
                }
                if (next_report_completed > beam->nstat) {
                    next_report_completed = beam->nstat;
                }
            }
        }
    }

    if (last_report_completed < beam->nstat) {
        double const t_now = monotonic_seconds();
        report_transport_progress(beam->nstat, beam->nstat, t_now - t_start);
    }

cleanup:
    osh_particle_pool_free(pool);
    osh_beam_runtime_free(beam_rt);
    return rc;
}

/* ---- Progress helpers ---------------------------------------------------- */

static double monotonic_seconds(void) {
#if defined(CLOCK_MONOTONIC)
    {
        struct timespec ts;
        if (clock_gettime(CLOCK_MONOTONIC, &ts) == 0) {
            return (double) ts.tv_sec + 1.0e-9 * (double) ts.tv_nsec;
        }
    }
#endif

#if defined(TIME_UTC)
    {
        struct timespec ts;
        if (timespec_get(&ts, TIME_UTC) == TIME_UTC) {
            return (double) ts.tv_sec + 1.0e-9 * (double) ts.tv_nsec;
        }
    }
#endif

    return (double) time(NULL);
}

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

static void report_transport_progress(size_t completed, size_t total, double elapsed_s) {
    double primaries_per_second;
    size_t remaining;

    remaining = (completed < total) ? (total - completed) : 0u;
    primaries_per_second = (elapsed_s > 0.0) ? ((double) completed / elapsed_s) : 0.0;
    osh_info("Transport progress: %zu/%zu completed, %zu left, %.1f primaries/s",
             completed,
             total,
             remaining,
             primaries_per_second);
}
