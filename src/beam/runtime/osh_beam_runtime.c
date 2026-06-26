#include "beam/runtime/osh_beam_runtime.h"

#include <stdint.h>
#include <stdlib.h>

#include "beam/osh_beam_model.h"
#include "beam/osh_beamdef.h"
#include "common/osh_ray.h"
#include "particle/osh_particle.h"

/* ---- Forward declarations ------------------------------------------------ */

static enum osh_status setup_spots(struct osh_beam_runtime *rt);
static enum osh_status setup_phsp(struct osh_beam_runtime *rt);
static enum osh_status fill_from_spots(struct osh_beam_runtime *rt,
                                       struct osh_rng_seeding const *seeding,
                                       struct osh_particle_pool *pool,
                                       size_t n,
                                       uint64_t global_prim_base);
static enum osh_status fill_from_phsp(struct osh_beam_runtime *rt,
                                      struct osh_rng_seeding const *seeding,
                                      struct osh_particle_pool *pool,
                                      size_t n,
                                      uint64_t global_prim_base);

/* ---- Lifecycle ----------------------------------------------------------- */

enum osh_status osh_beam_compile(struct beam_workspace const *workspace, struct osh_beam_runtime **rt_out) {
    struct osh_beam_runtime *rt;
    enum osh_status rc;

    if (!workspace || !rt_out) {
        return OSH_EINVAL;
    }

    rt = (struct osh_beam_runtime *) malloc(sizeof(*rt));
    if (!rt) {
        return OSH_ENOMEM;
    }

    rt->workspace = workspace;
    rt->primary = (struct particle) {0};
    rt->primaries_generated = 0u;

    switch (workspace->beam_mode) {
    case OSH_BEAM_MODE_SPOTS:
    case OSH_BEAM_MODE_SOBP:
        rc = setup_spots(rt);
        break;
    case OSH_BEAM_MODE_PHSP:
        rc = setup_phsp(rt);
        break;
    default:
        rc = OSH_ENOTSUP;
        break;
    }

    if (rc != OSH_OK) {
        free(rt);
        return rc;
    }

    *rt_out = rt;
    return OSH_OK;
}

void osh_beam_runtime_free(struct osh_beam_runtime **rt) {
    if (!rt || !*rt) {
        return;
    }

    if ((*rt)->workspace && (*rt)->workspace->beam_mode == OSH_BEAM_MODE_PHSP) {
        if ((*rt)->source.phsp.mcpl_handle) {
            /* TODO: mcpl_close((*rt)->source.phsp.mcpl_handle); */
            (*rt)->source.phsp.mcpl_handle = NULL;
        }
    }

    free(*rt);
    *rt = NULL;
}

/* ---- Primary generation -------------------------------------------------- */

enum osh_status osh_beam_runtime_fill_pool_at(struct osh_beam_runtime *rt,
                                              struct osh_rng_seeding const *seeding,
                                              struct osh_particle_pool *pool,
                                              size_t n,
                                              uint64_t global_prim_base) {
    if (!rt || !seeding || !pool) {
        return OSH_EINVAL;
    }
    if (n == 0u) {
        return OSH_OK;
    }
    if (pool->n > pool->capacity || n > pool->capacity - pool->n) {
        return OSH_EINVAL;
    }
    /* The last primary in this fill has prim_idx = global_prim_base + (n - 1)
     * and hist_index = seeding->hist_base + prim_idx.  Reject any range that
     * would wrap either value rather than silently reusing indices/RNG streams.
     * (Unreachable in practice -- it needs ~2^64 histories -- but this is the
     * contract for parallel fills, so fail loudly instead of corrupting tallies.) */
    {
        uint64_t const last_offset = (uint64_t) n - 1u;
        uint64_t first_hist_index;

        if (global_prim_base > UINT64_MAX - last_offset) {
            return OSH_EINVAL;
        }
        if (seeding->hist_base > UINT64_MAX - global_prim_base) {
            return OSH_EINVAL;
        }
        first_hist_index = seeding->hist_base + global_prim_base;
        if (first_hist_index > UINT64_MAX - last_offset) {
            return OSH_EINVAL;
        }
    }

    switch (rt->workspace->beam_mode) {
    case OSH_BEAM_MODE_SPOTS:
    case OSH_BEAM_MODE_SOBP:
        return fill_from_spots(rt, seeding, pool, n, global_prim_base);
    case OSH_BEAM_MODE_PHSP:
        return fill_from_phsp(rt, seeding, pool, n, global_prim_base);
    default:
        return OSH_ENOTSUP;
    }
}

enum osh_status osh_beam_runtime_fill_pool(struct osh_beam_runtime *rt,
                                           struct osh_rng_seeding const *seeding,
                                           struct osh_particle_pool *pool,
                                           size_t n) {
    enum osh_status rc;

    if (!rt) {
        return OSH_EINVAL;
    }

    /* Serial convenience wrapper: the global history base is the running cursor
     * of primaries this runtime has already emitted.  All the work — and all
     * validation — lives in _at(); only this wrapper reads or advances the
     * shared cursor, so the cursor-free _at() path stays safe to call from many
     * workers over disjoint, explicitly-based ranges.  The cursor is advanced
     * only on success, and _at() rejects a base that would wrap, so the cursor
     * can never silently overflow into reused indices. */
    rc = osh_beam_runtime_fill_pool_at(rt, seeding, pool, n, rt->primaries_generated);
    if (rc == OSH_OK) {
        rt->primaries_generated += (uint64_t) n;
    }
    return rc;
}

/* ---- Spots/SOBP implementation ------------------------------------------- */

static enum osh_status setup_spots(struct osh_beam_runtime *rt) {
    /*
     * No extra setup needed for the scalar spots path — sampling reads
     * directly from workspace->spots[] via osh_beam_new_primary().
     *
     * Future optimisation: pre-flatten per-spot scalars (t0, tsigma, size,
     * div, cor, tm) into SoA arrays in rt->source.spots so that fill_from_spots
     * can drive SIMD loops without chasing pointers into the AoS beam_spot
     * array.  See the comment in osh_beam_runtime.h for the planned layout.
     */
    rt->source.spots._reserved = 0;
    if (!rt->workspace->prepared) {
        return OSH_EINVAL;
    }
    if (!rt->workspace->has_primary || !osh_particle_from_pdg(&rt->primary, rt->workspace->primary.pdg)) {
        return OSH_EINVAL;
    }
    return OSH_OK;
}

static enum osh_status fill_from_spots(struct osh_beam_runtime *rt,
                                       struct osh_rng_seeding const *seeding,
                                       struct osh_particle_pool *pool,
                                       size_t n,
                                       uint64_t global_prim_base) {
    size_t i;
    size_t slot;
    struct ray_v ray;
    enum osh_status rc;

    /*
     * Scalar fill loop: one primary per iteration.
     *
     * Each primary owns two independent streams keyed by its global history
     * index: a transient BEAM stream that samples the source phase space, and
     * the slot's persistent PHYSICS stream (pool->rng[slot]) used later by the
     * transport loop.  Seeding by index keeps the source sampling identical
     * regardless of fill chunking, pool capacity, or physics options, and
     * keeps it independent of the transport stream.
     *
     * The global history index of primary i is global_prim_base + i, supplied
     * by the caller rather than read from an internal cursor.  That is what lets
     * two workers fill disjoint slices [a, b) and [b, c) of the same run — in
     * any order, on any thread — and still reproduce the single-pass per-history
     * streams exactly: the seed depends only on the index, never on call order.
     *
     * osh_beam_new_primary() returns one ray_v (AoS) at a time. We copy the
     * seven phase-space scalars into the pool's SoA arrays and assign
     * per-history metadata, reusing one resolved species descriptor per run.
     *
     * Future: replace with a vectorized kernel that samples directly into
     * the SoA arrays, eliminating the intermediate ray_v and enabling
     * auto-vectorization of the energy/position sampling arithmetic.
     */
    for (i = 0u; i < n; ++i) {
        uint64_t const prim_idx = global_prim_base + i;
        uint64_t const hist_index = seeding->hist_base + prim_idx;
        struct osh_rng beam_rng;

        slot = pool->n + i;

        osh_rng_seed_history(&beam_rng, seeding->type, seeding->seed, hist_index, OSH_RNG_PURPOSE_BEAM);
        rc = osh_beam_new_primary(rt->workspace, &beam_rng, &ray);
        if (rc != OSH_OK) {
            return rc;
        }

        pool->x[slot] = ray.p[0];
        pool->y[slot] = ray.p[1];
        pool->z[slot] = ray.p[2];
        pool->ux[slot] = ray.v[0];
        pool->uy[slot] = ray.v[1];
        pool->uz[slot] = ray.v[2];
        pool->e[slot] = ray.p[3];
        pool->wt[slot] = 1.0; /* SOBP spot weight is absorbed into the
                               * inverse-CDF spot selection; per-history
                               * weight remains 1.0 for unweighted transport */
        pool->prim_idx[slot] = prim_idx;
        pool->gen[slot] = 0u;
        pool->species[slot] = &rt->primary;

        /* Persistent transport stream for this slot, independent of BEAM. */
        osh_rng_seed_history(&pool->rng[slot], seeding->type, seeding->seed, hist_index, OSH_RNG_PURPOSE_PHYSICS);
    }

    pool->n += n;
    /* No rt->primaries_generated mutation here: the global base is caller-driven
     * (see osh_beam_runtime_fill_pool_at).  The serial wrapper owns the cursor. */
    return OSH_OK;
}

/* ---- PHSP stub ----------------------------------------------------------- */

static enum osh_status setup_phsp(struct osh_beam_runtime *rt) {
    /*
     * PHSP (MCPL) source — not yet implemented.
     *
     * When implemented, this function should:
     *   1. Open the MCPL file at workspace->phsp->fname.
     *   2. Read and validate the file header.
     *   3. Store the opaque handle in rt->source.phsp.mcpl_handle.
     *   4. Set rt->source.phsp.total to the number of entries.
     *   5. Set rt->source.phsp.file_pos to 0.
     *
     * TODO: link against libmcpl and implement.
     * Reference: https://mctools.github.io/mcpl/
     */
    rt->source.phsp.mcpl_handle = NULL;
    rt->source.phsp.file_pos = 0u;
    rt->source.phsp.total = 0u;
    return OSH_ENOTSUP;
}

static enum osh_status fill_from_phsp(struct osh_beam_runtime *rt,
                                      struct osh_rng_seeding const *seeding,
                                      struct osh_particle_pool *pool,
                                      size_t n,
                                      uint64_t global_prim_base) {
    /*
     * PHSP (MCPL) fill — not yet implemented.
     *
     * When implemented, this function should:
     *   1. Read min(n, total - file_pos) entries from the MCPL file.
     *   2. Convert each MCPL particle record to pool SoA layout.
     *   3. Advance rt->source.phsp.file_pos; rewind when exhausted.
     *   4. Assign prim_idx from global_prim_base + i and gen = 0.
     *   5. Seed pool->rng[slot] via osh_rng_seed_history(..., hist_index,
     *      OSH_RNG_PURPOSE_PHYSICS) exactly as the spots path does, so PHSP
     *      transport is reproducible and capacity-independent too.
     *
     * Key design constraint: MCPL files can exceed several GB, so entries
     * must be read in chunks — never all at once.  The pool capacity acts
     * as the chunk size; a single fill_pool call reads at most
     * pool->capacity entries from disk.
     *
     * TODO: implement once libmcpl integration is in place.
     */
    (void) rt;
    (void) seeding;
    (void) pool;
    (void) n;
    (void) global_prim_base;
    return OSH_ENOTSUP;
}
