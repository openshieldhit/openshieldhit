#include "transport/osh_transport.h"

#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "beam/osh_beam.h"
#include "beam/osh_beamdef.h"
#include "beam/runtime/osh_beam_runtime.h"
#include "common/osh_coord.h"
#include "common/osh_logger.h"
#include "common/osh_particle_pool.h"
#include "common/osh_ray.h"
#include "gemca/runtime/osh_gemca_runtime.h" /* replaces gemca/osh_gemca2.h — all gemca access goes through the runtime */
#include "material/osh_material.h"
#include "material/runtime/osh_material_runtime.h"
#include "physics/osh_physics_bethe.h"
#include "physics/osh_physics_moliere.h"
#include "physics/osh_physics_straggling.h"
#include "random/osh_rng.h"
#include "scoring/runtime/osh_scoring_step.h"

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
#define OSH_TRANSPORT_BOUNDARY_EPS 1e-8
#define OSH_TRANSPORT_POOL_CAPACITY 4096u
#define OSH_TRANSPORT_MAX_STEPS_PER_PRIMARY 1000000u
#define OSH_TRANSPORT_PROGRESS_MIN_INTERVAL_S 1.0
#define OSH_TRANSPORT_PROGRESS_MAX_INTERVAL_S 10.0
#define OSH_TRANSPORT_PROGRESS_MIN_CHUNK 1000u
#define OSH_TRANSPORT_PROGRESS_TARGET_CHUNKS 20u

/*
 * Maximum projected RMS scattering angle allowed in a single substep [rad].
 *
 * The Highland / Molière Gaussian model is reliable while θ₀ is small.
 * This threshold drives the s_theta substep criterion: steps are split so
 * that θ₀ ≤ OSH_TRANSPORT_THETA_MAX before MCS is sampled.  The value 0.1 rad
 * is a conservative limit well inside the Gaussian regime; for clinical beams
 * the DELTAE criterion typically dominates and s_theta is never reached.
 */
#define OSH_TRANSPORT_THETA_MAX_RAD 0.1

/* ---- Forward declarations ------------------------------------------------ */

static enum osh_status transport_step_one(struct osh_particle_pool *pool,
                                          size_t slot,
                                          size_t zone_idx,
                                          double boundary_ds,
                                          struct gemca_runtime const *geom_rt,
                                          struct beam_workspace const *beam,
                                          struct material_workspace const *materials,
                                          struct osh_material_runtime const *tables,
                                          struct osh_scoring_runtime *scoring,
                                          double deltae,
                                          struct osh_rng *rng);
static double cutoff_total_energy(struct beam_workspace const *beam,
                                  struct osh_material_runtime const *tables,
                                  struct particle const *part);
static double energy_from_residual_range(struct osh_material_runtime const *tables,
                                         size_t material_idx,
                                         size_t projectile_idx,
                                         double residual_range);
static double energy_grid_value(struct osh_material_runtime const *tables, size_t energy_idx);
static enum osh_status find_projectile_index(struct osh_material_runtime const *tables,
                                             struct particle const *part,
                                             size_t *projectile_idx_out);
static int is_blackhole_material(size_t material_idx);
static int is_vacuum_material(size_t material_idx);
static double monotonic_seconds(void);
static size_t transport_progress_chunk_size(size_t total);
static void report_transport_progress(size_t completed, size_t total, double elapsed_s);

static void step_from_pool(struct step *st,
                           struct osh_particle_pool const *pool,
                           size_t slot,
                           double step_len,
                           double exit_energy,
                           double rho,
                           int medium,
                           int zone);

/* ---- Public API ---------------------------------------------------------- */

/**
 * @brief Run the minimal straight-line CSDA transport loop.
 *
 * @details
 * Wavefront (BFS) transport loop: all live primaries advance one step per
 * round.  The outer loop runs until all beam->nstat primaries have been
 * generated and transported to death:
 *
 *   1. When the pool is empty and primaries remain, fill it from beam_runtime
 *      (up to OSH_TRANSPORT_POOL_CAPACITY primaries).
 *   2. Call transport_step_one() for every live pool slot.  Particles that die
 *      (energy cutoff, geometry exit, blackhole) are marked by zeroing e[slot].
 *   3. Compact the pool, removing dead entries.
 *   4. Repeat until primaries_done == beam->nstat and pool is empty.
 *
 * The pool capacity controls the trade-off between cache pressure (small pool)
 * and parallelism (large pool).  Physics results are independent of capacity.
 *
 * Current limitations:
 *   - no nuclear interactions, no secondaries
 *   - SOBP and single-spot beam modes only (PHSP returns OSH_ENOTSUP)
 */
enum osh_status osh_transport_run_minimal(struct beam_workspace const *beam,
                                          struct gemca_runtime const *geom_rt,
                                          struct material_workspace const *materials,
                                          struct osh_material_runtime const *tables,
                                          struct osh_scoring_runtime *scoring) {
    struct osh_rng rng;
    struct osh_beam_runtime *beam_rt = NULL;
    struct osh_particle_pool *pool = NULL;
    size_t zone_batch[OSH_TRANSPORT_POOL_CAPACITY]; /* stack-safe: capacity is compile-time constant */
    double dist_batch[OSH_TRANSPORT_POOL_CAPACITY]; /* stack-safe: ~32 KiB at 4096 doubles */
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

    if (!beam || !geom_rt || !materials || !tables || !scoring) {
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

        /* Batch geometry: zone lookup and boundary distance for all live particles.
         * The scratch arrays are stack-allocated — safe because pool->n is always
         * <= pool->capacity <= OSH_TRANSPORT_POOL_CAPACITY at this point. */
        osh_gemca_runtime_get_zone_batch(
            geom_rt, pool->x, pool->y, pool->z, pool->ux, pool->uy, pool->uz, pool->n, zone_batch);
        osh_gemca_runtime_get_distance_batch(
            geom_rt, pool->x, pool->y, pool->z, pool->ux, pool->uy, pool->uz, zone_batch, pool->n, dist_batch);

        /* Advance every live particle by one step using the precomputed geometry */
        for (i = 0u; i < pool->n; ++i) {
            if (steps_taken >= step_budget) {
                osh_error("transport: step budget exceeded after %zu steps (pool slot %zu, primaries_done=%zu, zone=%zu, e=%.17g, pos=(%.17g, %.17g, %.17g), dir=(%.17g, %.17g, %.17g))",
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
                rc = OSH_ESTATE; /* step budget exceeded — likely stuck particle */
                goto cleanup;
            }
            rc = transport_step_one(pool,
                                    i,
                                    zone_batch[i],
                                    dist_batch[i],
                                    geom_rt,
                                    beam,
                                    materials,
                                    tables,
                                    scoring,
                                    (double) beam->deltae,
                                    &rng);
            if (rc != OSH_OK) {
                osh_error("transport: slot %zu failed with rc=%d zone=%zu boundary_ds=%.17g e=%.17g pos=(%.17g, %.17g, %.17g) dir=(%.17g, %.17g, %.17g)",
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

        /* Compact dead entries (e[i] <= 0) so the next fill appends cleanly */
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

/* ---- Per-slot transport step --------------------------------------------- */

/**
 * @brief Advance particle at @p slot by one straight-line CSDA step.
 *
 * @details
 * Reads the particle's phase-space state from the pool SoA arrays, computes
 * one CSDA step bounded by either a DELTAE energy-loss criterion or the next
 * geometry boundary, scores the step, and writes the updated position and
 * energy back into the pool.
 *
 * Termination: when the particle reaches the energy cutoff, exits the geometry,
 * or enters a blackhole material, its energy slot is zeroed so that the next
 * osh_particle_pool_compact() call removes it.  These physics-termination
 * conditions are handled silently (OSH_OK return) so that one bad particle
 * does not abort the entire run.  Only allocation or I/O failures from scoring
 * propagate as non-OK status.
 *
 * Coordinate system: beam primaries are always placed in OSH_COORD_UNIVERSE by
 * osh_beam_runtime_fill_pool().  Secondaries (future) inherit their parent's
 * system at birth.  The pool does not store the system field; this function
 * always reconstructs the geometry ray with system == OSH_COORD_UNIVERSE,
 * which is correct for beam-primary histories in the current single-species
 * implementation.
 *
 * Boundary nudge: when a particle sits exactly on a boundary
 * (distance <= OSH_TRANSPORT_BOUNDARY_EPS), it is nudged one epsilon along its
 * direction and the function returns without scoring a zero-length step.  The
 * particle remains alive; transport_step_one() is called again for it in the
 * next wavefront round.
 *
 * Random-hinge MCS (physics detail):
 * Multiple Coulomb scattering is applied via the random hinge method rather
 * than at the end of the step.  A hinge point h is sampled uniformly in
 * [0, step_len].  The particle travels straight from p to the hinge along the
 * incident direction u0, then is deflected by Gaussian angles (tx, ty) drawn
 * from the Highland distribution, and continues in the scattered direction u1
 * for the remaining (step_len - h).  The actual exit position is therefore:
 *
 *   q = p + h·u0 + tail·u1
 *
 * This avoids the parallel-boundary pathology of the 2-pass approach (where
 * re-querying geometry in the scattered direction can give infinite boundary
 * distances for backscatter) and correctly distributes the deflection
 * throughout the step rather than applying it entirely at the endpoint.
 *
 * To keep every scored step inside a single GEMCA zone, the post-hinge leg is
 * clipped by a second distance query from the hinge point along u1:
 *
 *   tail = min(step_len - h, get_dist(ph, u1))
 *
 * where ph = p + h·u0.  The actual transported track length is then
 *
 *   ds = h + tail
 *
 * Three substep criteria are applied before sampling the hinge:
 *   1. boundary distance (geometric)
 *   2. ds_csda  — DELTAE fraction of CSDA range
 *   3. ds_theta — θ₀ ≤ OSH_TRANSPORT_THETA_MAX_RAD (keeps Highland valid)
 */
static enum osh_status transport_step_one(struct osh_particle_pool *pool,
                                          size_t slot,
                                          size_t zone_idx,
                                          double boundary_ds,
                                          struct gemca_runtime const *geom_rt,
                                          struct beam_workspace const *beam,
                                          struct material_workspace const *materials,
                                          struct osh_material_runtime const *tables,
                                          struct osh_scoring_runtime *scoring,
                                          double deltae,
                                          struct osh_rng *rng) {
    struct particle const *part;
    size_t projectile_idx;
    double cutoff_total;
    struct gemca_rt_zone const *zone;
    struct material const *material;
    double rho;
    double step_len;
    double exit_energy;
    double requested_step_len;
    double a_proj;
    double e0_total;
    double e1_target;
    double r0;
    double r1_csda;   /* CSDA range at the DELTAE target energy */
    double ds_csda;   /* DELTAE substep length [cm] */
    double ds_theta;  /* MCS substep length [cm] */
    double residual_range;
    double h;         /* hinge point [cm] along incident direction */
    double tail_len;  /* post-hinge leg [cm] along scattered direction */
    double qx;        /* exit position — bent-path endpoint */
    double qy;
    double qz;
    struct step st;
    int hit_boundary;
    int requested_hit_boundary;
    int requested_is_csda;
    enum osh_status rc;
    /* Per-material atomic scalars (loaded in non-vacuum branch) */
    double mat_z_mean;
    double mat_z_over_a;
    double mat_x0_gcm2;
    double proj_mass_mev;
    /* Mid-step physics quantities */
    double z_eff_0;   /* z_eff at entry energy — used for s_theta */
    double e_mid;
    double z_eff;
    double ds_gcm2;
    double proposed_ds_gcm2;
    double proposed_exit_energy;
    double sigma_strag;
    double theta0;
    double boundary_tail_ds;
    double w_scat[3]; /* exit direction after MCS (= u0 if no scatter) */
    double nudge_dir[3];
    struct ray hinge_ray;
    int enable_mcs;
    int enable_straggling;

    part = pool->species[slot];
    a_proj = (part->a > 0u) ? (double) part->a : 1.0;
    e0_total = pool->e[slot];

    rc = find_projectile_index(tables, part, &projectile_idx);
    if (rc != OSH_OK) {
        pool->e[slot] = 0.0; /* unknown species — kill silently */
        return OSH_OK;
    }

    cutoff_total = cutoff_total_energy(beam, tables, part);
    if (pool->e[slot] <= cutoff_total) {
        pool->e[slot] = 0.0;
        return OSH_OK;
    }

    if (zone_idx == OSH_GEMCA_ZONE_INDEX_INVALID) {
        pool->e[slot] = 0.0; /* escaped geometry */
        return OSH_OK;
    }

    zone = &geom_rt->zones[zone_idx];
    if (is_blackhole_material(zone->material_idx)) {
        pool->e[slot] = 0.0;
        return OSH_OK;
    }

    material = osh_material_by_index(materials, zone->material_idx);
    if (!material) {
        return OSH_ESTATE;
    }

    if (boundary_ds < 0.0) {
        return OSH_ESTATE;
    }
    if (boundary_ds <= OSH_TRANSPORT_BOUNDARY_EPS) {
        /* Nudge past the boundary; step again in the next wavefront round */
        pool->x[slot] += pool->ux[slot] * OSH_TRANSPORT_BOUNDARY_EPS;
        pool->y[slot] += pool->uy[slot] * OSH_TRANSPORT_BOUNDARY_EPS;
        pool->z[slot] += pool->uz[slot] * OSH_TRANSPORT_BOUNDARY_EPS;
        return OSH_OK;
    }

    rho = (material->rho > 0.0) ? material->rho : 0.0;

    /* Exit direction default: incident direction (no scatter in vacuum) */
    w_scat[0] = pool->ux[slot];
    w_scat[1] = pool->uy[slot];
    w_scat[2] = pool->uz[slot];
    nudge_dir[0] = pool->ux[slot];
    nudge_dir[1] = pool->uy[slot];
    nudge_dir[2] = pool->uz[slot];
    h = 0.0;
    tail_len = 0.0;
    requested_hit_boundary = 0;
    requested_step_len = 0.0;
    requested_is_csda = 0;
    hit_boundary = 0;
    enable_mcs = (beam && beam->scatter != OSH_BEAM_MSCAT_OFF);
    enable_straggling = (beam && beam->straggl != OSH_BEAM_STRAGG_OFF);

    if (is_vacuum_material(zone->material_idx) || rho <= 0.0) {
        /* Vacuum: straight shot to next boundary, no energy loss or scatter */
        step_len = boundary_ds;
        requested_step_len = step_len;
        requested_hit_boundary = 1;
        exit_energy = e0_total;
        h = step_len;
        tail_len = 0.0;
        hit_boundary = 1;
    } else {
        /*
         * ---- Step-length determination: three substep criteria -------------
         *
         * The substep is the minimum of:
         *   1. boundary_ds   — next geometry boundary in incident direction
         *   2. ds_csda       — DELTAE fraction of the CSDA range
         *   3. ds_theta      — maximum step keeping θ₀ ≤ OSH_TRANSPORT_THETA_MAX_RAD
         *
         * After the minimum is chosen, exit energy is recovered from the
         * CSDA residual range for the actual step_len.  Using residual range
         * for all three cases is slightly less accurate than using e1_target
         * when step_len == ds_csda (round-trip through inverse range table),
         * but the difference is below the energy grid spacing and avoids a
         * branch in the exit-energy computation.
         */
        mat_z_mean    = (double) tables->z_mean[zone->material_idx];
        mat_z_over_a  = (double) tables->z_over_a[zone->material_idx];
        mat_x0_gcm2   = (double) tables->rad_length[zone->material_idx];
        proj_mass_mev = tables->projectile_mass_mev[projectile_idx];

        e1_target = e0_total * (1.0 - deltae);
        if (e1_target < cutoff_total) {
            e1_target = cutoff_total;
        }
        if (e1_target >= e0_total) {
            pool->e[slot] = 0.0;
            return OSH_OK;
        }

        r0      = osh_material_runtime_range_lookup(tables, zone->material_idx, projectile_idx, e0_total / a_proj);
        r1_csda = osh_material_runtime_range_lookup(tables, zone->material_idx, projectile_idx, e1_target / a_proj);
        ds_csda = (r0 - r1_csda) / rho;
        if (ds_csda <= 0.0) {
            pool->e[slot] = 0.0;
            return OSH_OK;
        }

        /* MCS substep criterion evaluated at entry energy */
        if (enable_mcs) {
            z_eff_0  = osh_physics_bethe_z_eff(e0_total / a_proj, (double) part->z, a_proj, mat_z_mean);
            ds_theta = osh_physics_moliere_s_theta(e0_total, proj_mass_mev, z_eff_0,
                                                   rho, mat_x0_gcm2, OSH_TRANSPORT_THETA_MAX_RAD);
        } else {
            ds_theta = 0.0;
        }

        requested_step_len = boundary_ds;
        if (ds_csda < requested_step_len) {
            requested_step_len = ds_csda;
        }
        if (ds_theta > 0.0 && ds_theta < requested_step_len) {
            requested_step_len = ds_theta;
        }
        requested_hit_boundary = (requested_step_len >= boundary_ds - OSH_TRANSPORT_BOUNDARY_EPS);
        requested_is_csda = (fabs(requested_step_len - ds_csda) <= 1.0e-12 * fmax(1.0, requested_step_len));

        /*
         * Geometry-limited step: stay on the entry chord to the boundary and
         * apply the direction change at the endpoint.
         *
         * Physics-limited step: sample a random hinge, scatter there, then
         * clip the post-hinge leg with a second get_distance() so the final
         * endpoint remains inside the current zone.
         */
        if (requested_hit_boundary) {
            h = requested_step_len;
            tail_len = 0.0;
            step_len = requested_step_len;
            hit_boundary = 1;
        } else {
            h = requested_step_len * osh_rng_double(rng);
            tail_len = requested_step_len - h;

            residual_range = r0 - rho * requested_step_len;
            if (residual_range <= 0.0) {
                pool->e[slot] = 0.0;
                return OSH_OK;
            }
            if (requested_is_csda) {
                proposed_exit_energy = e1_target;
            } else {
                proposed_exit_energy =
                    energy_from_residual_range(tables, zone->material_idx, projectile_idx, residual_range) * a_proj;
                if (proposed_exit_energy > e0_total) {
                    proposed_exit_energy = e0_total;
                }
                if (proposed_exit_energy < cutoff_total) {
                    proposed_exit_energy = cutoff_total;
                }
            }

            e_mid = 0.5 * (e0_total + proposed_exit_energy);
            z_eff = osh_physics_bethe_z_eff(e_mid / a_proj, (double) part->z, a_proj, mat_z_mean);
            proposed_ds_gcm2 = rho * requested_step_len;

            if (enable_mcs && mat_x0_gcm2 > 0.0) {
                theta0 = osh_physics_moliere_theta0(e_mid, proj_mass_mev, z_eff, proposed_ds_gcm2, mat_x0_gcm2);
                if (theta0 > 0.0) {
                    double const v_in[3] = {pool->ux[slot], pool->uy[slot], pool->uz[slot]};
                    osh_physics_moliere_scatter(v_in, w_scat, theta0, rng);
                }
            }

            hinge_ray.p[0] = pool->x[slot] + pool->ux[slot] * h;
            hinge_ray.p[1] = pool->y[slot] + pool->uy[slot] * h;
            hinge_ray.p[2] = pool->z[slot] + pool->uz[slot] * h;
            hinge_ray.cp[0] = w_scat[0];
            hinge_ray.cp[1] = w_scat[1];
            hinge_ray.cp[2] = w_scat[2];
            hinge_ray.system = OSH_COORD_UNIVERSE;

            boundary_tail_ds = osh_gemca_runtime_get_distance(geom_rt, zone_idx, &hinge_ray);
            if (boundary_tail_ds < 0.0) {
                osh_error("transport: negative hinge boundary distance zone=%zu h=%.17g tail_req=%.17g hinge=(%.17g, %.17g, %.17g) w=(%.17g, %.17g, %.17g)",
                          zone_idx,
                          h,
                          tail_len,
                          hinge_ray.p[0],
                          hinge_ray.p[1],
                          hinge_ray.p[2],
                          hinge_ray.cp[0],
                          hinge_ray.cp[1],
                          hinge_ray.cp[2]);
                return OSH_ESTATE;
            }
            if (boundary_tail_ds < tail_len) {
                tail_len = boundary_tail_ds;
                hit_boundary = 1;
                nudge_dir[0] = w_scat[0];
                nudge_dir[1] = w_scat[1];
                nudge_dir[2] = w_scat[2];
            }
            step_len = h + tail_len;
        }

        /* Actual exit energy and straggling use the clipped in-zone step. */
        residual_range = r0 - rho * step_len;
        if (residual_range <= 0.0) {
            pool->e[slot] = 0.0;
            return OSH_OK;
        }
        if (requested_is_csda && fabs(step_len - requested_step_len) <= 1.0e-12 * fmax(1.0, requested_step_len)) {
            exit_energy = e1_target;
        } else {
            exit_energy = energy_from_residual_range(tables, zone->material_idx, projectile_idx, residual_range) * a_proj;
            if (exit_energy > e0_total) {
                exit_energy = e0_total;
            }
            if (exit_energy < cutoff_total) {
                exit_energy = cutoff_total;
            }
        }

        e_mid = 0.5 * (e0_total + exit_energy);
        z_eff = osh_physics_bethe_z_eff(e_mid / a_proj, (double) part->z, a_proj, mat_z_mean);
        ds_gcm2 = rho * step_len;

        /* Gaussian energy straggling (Bohr variance) */
        sigma_strag = enable_straggling ? osh_physics_straggling_sigma(z_eff, mat_z_over_a, ds_gcm2) : 0.0;
        if (sigma_strag > 0.0) {
            exit_energy += osh_rng_gauss(rng, 0.0, sigma_strag);
            if (exit_energy > e0_total) {
                exit_energy = e0_total;
            }
            if (exit_energy < cutoff_total) {
                exit_energy = cutoff_total;
            }
        }

        if (requested_hit_boundary && enable_mcs && mat_x0_gcm2 > 0.0) {
            theta0 = osh_physics_moliere_theta0(e_mid, proj_mass_mev, z_eff, ds_gcm2, mat_x0_gcm2);
            if (theta0 > 0.0) {
                double const v_in[3] = {pool->ux[slot], pool->uy[slot], pool->uz[slot]};
                osh_physics_moliere_scatter(v_in, w_scat, theta0, rng);
            }
        }
    }

    if (step_len <= 0.0) {
        pool->ux[slot] = w_scat[0];
        pool->uy[slot] = w_scat[1];
        pool->uz[slot] = w_scat[2];
        if (hit_boundary) {
            pool->x[slot] += nudge_dir[0] * OSH_TRANSPORT_BOUNDARY_EPS;
            pool->y[slot] += nudge_dir[1] * OSH_TRANSPORT_BOUNDARY_EPS;
            pool->z[slot] += nudge_dir[2] * OSH_TRANSPORT_BOUNDARY_EPS;
        }
        return OSH_OK;
    }

    /*
     * Exit position: two-segment bent path.
     *
     *   q = p + h·u0 + tail·u1
 *
     * For vacuum (h == 0, w_scat == u0) this reduces to q = p + step_len·u0.
     * For material steps the hinge bends the path; q is no longer collinear
     * with p and u0 but ds = step_len is still the correct track length.
     *
     * Scoring uses p, q, v=u0, w=u1, ds, de — no hinge position stored.
     * Zone-based scorers deposit de in the zone containing p (single zone per
     * step, guaranteed by the boundary criterion).  Voxel DDA scorers trace
     * the chord p→q; the approximation error is O(θ₀² × step_len) ≈ 10⁻⁵ cm
     * for clinical beams, well below voxel resolution.
     */
    qx = pool->x[slot] + pool->ux[slot] * h + w_scat[0] * tail_len;
    qy = pool->y[slot] + pool->uy[slot] * h + w_scat[1] * tail_len;
    qz = pool->z[slot] + pool->uz[slot] * h + w_scat[2] * tail_len;

    step_from_pool(&st, pool, slot, step_len, exit_energy, rho, (int) zone->material_idx, (int) zone_idx);
    st.q[0] = qx;
    st.q[1] = qy;
    st.q[2] = qz;
    st.w[0] = w_scat[0];
    st.w[1] = w_scat[1];
    st.w[2] = w_scat[2];

    rc = osh_scoring_score_step(scoring, part, &st);
    if (rc != OSH_OK) {
        osh_error("transport: scoring rejected step rc=%d ds=%.17g p=(%.17g, %.17g, %.17g) q=(%.17g, %.17g, %.17g) v=(%.17g, %.17g, %.17g) w=(%.17g, %.17g, %.17g)",
                  (int) rc,
                  st.ds,
                  st.p[0],
                  st.p[1],
                  st.p[2],
                  st.q[0],
                  st.q[1],
                  st.q[2],
                  st.v[0],
                  st.v[1],
                  st.v[2],
                  st.w[0],
                  st.w[1],
                  st.w[2]);
        return rc;
    }

    /* Update pool: position to bent endpoint, energy, heading to u1 */
    pool->x[slot] = qx;
    pool->y[slot] = qy;
    pool->z[slot] = qz;
    pool->e[slot] = exit_energy;
    pool->ux[slot] = w_scat[0];
    pool->uy[slot] = w_scat[1];
    pool->uz[slot] = w_scat[2];

    if (hit_boundary) {
        /*
         * Nudge past the boundary in the direction of the leg that actually
         * reached it:
         *   - straight / geometry-limited boundary steps use u0
         *   - hinge-clipped steps use u1
         *
         * This avoids trapping a boundary-limited step in the old zone when an
         * end-of-step scatter happens to backscatter into the medium.
         */
        pool->x[slot] += nudge_dir[0] * OSH_TRANSPORT_BOUNDARY_EPS;
        pool->y[slot] += nudge_dir[1] * OSH_TRANSPORT_BOUNDARY_EPS;
        pool->z[slot] += nudge_dir[2] * OSH_TRANSPORT_BOUNDARY_EPS;
    }

    return OSH_OK;
}

/* ---- Physics helpers ----------------------------------------------------- */

static double cutoff_total_energy(struct beam_workspace const *beam,
                                  struct osh_material_runtime const *tables,
                                  struct particle const *part) {
    double a_proj;
    double cutoff_total;
    double cutoff_from_transport;
    double cutoff_from_beam;

    a_proj = (part->a > 0u) ? (double) part->a : 1.0;
    cutoff_total = OSH_BEAM_TMIN;
    cutoff_from_beam = 0.0;
    if (beam) {
        cutoff_from_beam = (double) beam->tcut * a_proj;
        if (cutoff_from_beam > cutoff_total) {
            cutoff_total = cutoff_from_beam;
        }
    }
    cutoff_from_transport = tables->emin * a_proj;
    if (cutoff_from_transport > cutoff_total) {
        cutoff_total = cutoff_from_transport;
    }
    return cutoff_total;
}

static double energy_from_residual_range(struct osh_material_runtime const *tables,
                                         size_t material_idx,
                                         size_t projectile_idx,
                                         double residual_range) {
    float const *range_col;
    size_t lo;
    size_t hi;
    size_t mid;
    double r_lo;
    double r_hi;
    double e_lo;
    double e_hi;
    double frac;

    range_col = tables->range_csda + (material_idx * tables->nprojectiles + projectile_idx) * tables->nenergy;
    if (residual_range <= (double) range_col[0]) {
        return tables->emin;
    }
    if (residual_range >= (double) range_col[tables->nenergy - 1u]) {
        return tables->emax;
    }

    lo = 0u;
    hi = tables->nenergy - 1u;
    while (hi - lo > 1u) {
        mid = lo + (hi - lo) / 2u;
        if ((double) range_col[mid] <= residual_range) {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    r_lo = (double) range_col[lo];
    r_hi = (double) range_col[hi];
    e_lo = energy_grid_value(tables, lo);
    e_hi = energy_grid_value(tables, hi);
    if (r_hi <= r_lo) {
        return e_lo;
    }
    frac = (residual_range - r_lo) / (r_hi - r_lo);
    return e_lo * (1.0 - frac) + e_hi * frac;
}

static double energy_grid_value(struct osh_material_runtime const *tables, size_t energy_idx) {
    return exp(tables->log_emin + (double) energy_idx / tables->inv_dlog);
}

static enum osh_status find_projectile_index(struct osh_material_runtime const *tables,
                                             struct particle const *part,
                                             size_t *projectile_idx_out) {
    unsigned int z_match;
    unsigned int a_match;
    size_t i;

    if (!tables || !part || !projectile_idx_out) {
        return OSH_EINVAL;
    }
    z_match = part->z;
    a_match = part->a;
    if (part->pdg == 2212 && z_match == 0u && a_match == 0u) {
        z_match = 1u;
        a_match = 1u;
    }
    for (i = 0; i < tables->nprojectiles; ++i) {
        if (tables->projectile_z[i] == z_match && tables->projectile_a[i] == a_match) {
            *projectile_idx_out = i;
            return OSH_OK;
        }
    }
    for (i = 0; i < tables->nprojectiles; ++i) {
        if (tables->projectile_z[i] == z_match) {
            *projectile_idx_out = i;
            return OSH_OK;
        }
    }
    return OSH_ENOTSUP;
}

static int is_blackhole_material(size_t material_idx) {
    return material_idx == OSH_MATERIAL_INDEX_BLACKHOLE;
}

static int is_vacuum_material(size_t material_idx) {
    return material_idx == OSH_MATERIAL_INDEX_VACUUM;
}

static double monotonic_seconds(void) {
    struct timespec ts;

    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double) ts.tv_sec + 1.0e-9 * (double) ts.tv_nsec;
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

/* ---- Pool ↔ geometry/scoring adapters ------------------------------------ */

/**
 * @brief Build a scoring step from pool slot @p slot and step kinematics.
 *
 * @details
 * Reads per-history metadata (wt, prim_idx, gen) directly from the pool SoA
 * arrays.  This properly wires up the statistical weight and history context
 * that scoring filters and tallies need.  With the old DFS design these fields
 * were placeholder constants (wt=1.0, prim_idx=0, gen=0); they are now
 * populated from the pool that was filled by osh_beam_runtime_fill_pool().
 */
static void step_from_pool(struct step *st,
                           struct osh_particle_pool const *pool,
                           size_t slot,
                           double step_len,
                           double exit_energy,
                           double rho,
                           int medium,
                           int zone) {
    st->p[0] = pool->x[slot];
    st->p[1] = pool->y[slot];
    st->p[2] = pool->z[slot];
    st->p[3] = pool->e[slot];
    st->q[0] = pool->x[slot] + pool->ux[slot] * step_len;
    st->q[1] = pool->y[slot] + pool->uy[slot] * step_len;
    st->q[2] = pool->z[slot] + pool->uz[slot] * step_len;
    st->q[3] = exit_energy;
    st->v[0] = pool->ux[slot];
    st->v[1] = pool->uy[slot];
    st->v[2] = pool->uz[slot];
    st->w[0] = pool->ux[slot];
    st->w[1] = pool->uy[slot];
    st->w[2] = pool->uz[slot];
    st->ds = step_len;
    st->de = pool->e[slot] - exit_energy;
    if (st->de < 0.0) {
        st->de = 0.0;
    }
    st->rho = rho;
    st->wt = pool->wt[slot];
    st->medium = medium;
    st->zone = zone;
    st->system = OSH_COORD_UNIVERSE;
    st->prim_idx = pool->prim_idx[slot];
    st->gen = pool->gen[slot];
}
