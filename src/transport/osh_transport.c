#include "transport/osh_transport.h"

#include <stdlib.h>

#include "beam/osh_beam.h"
#include "beam/osh_beamdef.h"
#include "beam/runtime/osh_beam_runtime.h"
#include "common/osh_coord.h"
#include "common/osh_particle_pool.h"
#include "common/osh_ray.h"
#include "gemca/runtime/osh_gemca_runtime.h" /* replaces gemca/osh_gemca2.h — all gemca access goes through the runtime */
#include "material/osh_material.h"
#include "material/runtime/osh_material_runtime.h"
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

/* ---- Forward declarations ------------------------------------------------ */

static enum osh_status transport_step_one(struct osh_particle_pool *pool,
                                          size_t slot,
                                          struct gemca_runtime *geom_rt,
                                          struct beam_workspace const *beam,
                                          struct material_workspace const *materials,
                                          struct osh_material_runtime const *tables,
                                          struct osh_scoring_runtime *scoring,
                                          double deltae);
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
static void ray_from_pool(struct ray *ray, struct osh_particle_pool const *pool, size_t slot);
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
 * Current limitations (straight-line CSDA only, no secondaries):
 *   - no multiple scattering, no energy straggling
 *   - no nuclear interactions, no secondaries
 *   - SOBP and single-spot beam modes only (PHSP returns OSH_ENOTSUP)
 */
enum osh_status osh_transport_run_minimal(struct beam_workspace const *beam,
                                          struct gemca_runtime *geom_rt,
                                          struct material_workspace const *materials,
                                          struct osh_material_runtime const *tables,
                                          struct osh_scoring_runtime *scoring) {
    struct osh_rng rng;
    struct osh_beam_runtime *beam_rt = NULL;
    struct osh_particle_pool *pool = NULL;
    size_t capacity;
    size_t primaries_done;
    size_t n_fill;
    size_t i;
    size_t step_budget;
    size_t steps_taken;
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

        /* Advance every live particle by one step */
        for (i = 0u; i < pool->n; ++i) {
            if (steps_taken >= step_budget) {
                rc = OSH_ESTATE; /* step budget exceeded — likely stuck particle */
                goto cleanup;
            }
            rc = transport_step_one(pool, i, geom_rt, beam, materials, tables, scoring, (double) beam->deltae);
            if (rc != OSH_OK) {
                goto cleanup;
            }
            ++steps_taken;
        }

        /* Compact dead entries (e[i] <= 0) so the next fill appends cleanly */
        osh_particle_pool_compact(pool);
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
 */
static enum osh_status transport_step_one(struct osh_particle_pool *pool,
                                          size_t slot,
                                          struct gemca_runtime *geom_rt,
                                          struct beam_workspace const *beam,
                                          struct material_workspace const *materials,
                                          struct osh_material_runtime const *tables,
                                          struct osh_scoring_runtime *scoring,
                                          double deltae) {
    struct particle const *part;
    struct ray ray;
    size_t projectile_idx;
    double cutoff_total;
    size_t zone_idx;
    struct gemca_rt_zone const *zone;
    struct material const *material;
    double boundary_ds;
    double rho;
    double step_len;
    double exit_energy;
    double a_proj;
    double e0_total;
    double e1_target;
    double r0;
    double r1_target;
    double ds_csda;
    double residual_range;
    struct step st;
    int hit_boundary;
    enum osh_status rc;

    part = pool->species[slot];

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

    ray_from_pool(&ray, pool, slot);
    zone_idx = osh_gemca_runtime_get_zone(geom_rt, &ray);
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

    boundary_ds = osh_gemca_runtime_get_distance(geom_rt, zone_idx, &ray);
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
    if (is_vacuum_material(zone->material_idx) || rho <= 0.0) {
        step_len = boundary_ds;
        exit_energy = pool->e[slot];
    } else {
        a_proj = (part->a > 0u) ? (double) part->a : 1.0;
        e0_total = pool->e[slot];
        e1_target = e0_total * (1.0 - deltae);
        if (e1_target < cutoff_total) {
            e1_target = cutoff_total;
        }
        if (e1_target >= e0_total) {
            pool->e[slot] = 0.0;
            return OSH_OK;
        }

        r0 = osh_material_runtime_range_lookup(tables, zone->material_idx, projectile_idx, e0_total / a_proj);
        r1_target = osh_material_runtime_range_lookup(tables, zone->material_idx, projectile_idx, e1_target / a_proj);
        ds_csda = (r0 - r1_target) / rho;
        if (ds_csda <= 0.0) {
            pool->e[slot] = 0.0;
            return OSH_OK;
        }

        if (ds_csda <= boundary_ds) {
            step_len = ds_csda;
            exit_energy = e1_target;
        } else {
            step_len = boundary_ds;
            residual_range = r0 - rho * step_len;
            if (residual_range < 0.0) {
                residual_range = 0.0;
            }
            exit_energy =
                energy_from_residual_range(tables, zone->material_idx, projectile_idx, residual_range) * a_proj;
            if (exit_energy < cutoff_total) {
                exit_energy = cutoff_total;
            }
        }
    }

    if (step_len <= 0.0) {
        pool->e[slot] = 0.0;
        return OSH_OK;
    }
    if (exit_energy > pool->e[slot]) {
        exit_energy = pool->e[slot];
    }

    step_from_pool(&st, pool, slot, step_len, exit_energy, rho, (int) zone->material_idx, (int) zone_idx);
    rc = osh_scoring_score_step(scoring, part, &st);
    if (rc != OSH_OK) {
        return rc;
    }

    /* Advance position and update energy in the pool SoA arrays */
    hit_boundary = (step_len >= boundary_ds - OSH_TRANSPORT_BOUNDARY_EPS);
    pool->x[slot] += pool->ux[slot] * step_len;
    pool->y[slot] += pool->uy[slot] * step_len;
    pool->z[slot] += pool->uz[slot] * step_len;
    pool->e[slot] = exit_energy;

    if (hit_boundary) {
        pool->x[slot] += pool->ux[slot] * OSH_TRANSPORT_BOUNDARY_EPS;
        pool->y[slot] += pool->uy[slot] * OSH_TRANSPORT_BOUNDARY_EPS;
        pool->z[slot] += pool->uz[slot] * OSH_TRANSPORT_BOUNDARY_EPS;
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

/* ---- Pool ↔ geometry/scoring adapters ------------------------------------ */

/**
 * @brief Build a geometry ray from pool slot @p slot.
 *
 * @details
 * All beam primaries are generated in OSH_COORD_UNIVERSE by
 * osh_beam_runtime_fill_pool().  The pool does not store the coordinate system
 * field; it is hardcoded here as OSH_COORD_UNIVERSE.  When secondary particle
 * support is added, secondaries born in other coordinate systems will require
 * either an explicit system field in the pool or a per-entry flag.
 */
static void ray_from_pool(struct ray *ray, struct osh_particle_pool const *pool, size_t slot) {
    ray->p[0] = pool->x[slot];
    ray->p[1] = pool->y[slot];
    ray->p[2] = pool->z[slot];
    ray->cp[0] = pool->ux[slot];
    ray->cp[1] = pool->uy[slot];
    ray->cp[2] = pool->uz[slot];
    ray->system = OSH_COORD_UNIVERSE;
}

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
