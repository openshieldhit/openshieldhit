#include "beam/runtime/osh_beam_runtime.h"

#include <stdlib.h>

#include "beam/osh_beam_model.h"
#include "beam/osh_beamdef.h"
#include "common/osh_ray.h"
#include "particle/osh_particle.h"

/* ---- Forward declarations ------------------------------------------------ */

static enum osh_status setup_spots(struct osh_beam_runtime *rt);
static enum osh_status setup_phsp(struct osh_beam_runtime *rt);
static enum osh_status
fill_from_spots(struct osh_beam_runtime *rt, struct osh_rng *rng, struct osh_particle_pool *pool, size_t n);
static enum osh_status
fill_from_phsp(struct osh_beam_runtime *rt, struct osh_rng *rng, struct osh_particle_pool *pool, size_t n);

/* ---- Lifecycle ----------------------------------------------------------- */

enum osh_status osh_beam_runtime_setup(struct beam_workspace const *workspace, struct osh_beam_runtime **rt_out) {
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
    rt->primary = (struct particle){0};
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

void osh_beam_runtime_free(struct osh_beam_runtime *rt) {
    if (!rt) {
        return;
    }

    /*
     * For PHSP mode: close the MCPL file handle when implemented.
     * Currently a no-op since setup_phsp returns OSH_ENOTSUP and no handle
     * is ever opened.
     */
    if (rt->workspace && rt->workspace->beam_mode == OSH_BEAM_MODE_PHSP) {
        if (rt->source.phsp.mcpl_handle) {
            /* TODO: mcpl_close(rt->source.phsp.mcpl_handle); */
            rt->source.phsp.mcpl_handle = NULL;
        }
    }

    free(rt);
}

/* ---- Primary generation -------------------------------------------------- */

enum osh_status
osh_beam_runtime_fill_pool(struct osh_beam_runtime *rt, struct osh_rng *rng, struct osh_particle_pool *pool, size_t n) {
    if (!rt || !rng || !pool) {
        return OSH_EINVAL;
    }
    if (n == 0u) {
        return OSH_OK;
    }
    if (pool->n + n > pool->capacity) {
        return OSH_EINVAL;
    }

    switch (rt->workspace->beam_mode) {
    case OSH_BEAM_MODE_SPOTS:
    case OSH_BEAM_MODE_SOBP:
        return fill_from_spots(rt, rng, pool, n);
    case OSH_BEAM_MODE_PHSP:
        return fill_from_phsp(rt, rng, pool, n);
    default:
        return OSH_ENOTSUP;
    }
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

static enum osh_status
fill_from_spots(struct osh_beam_runtime *rt, struct osh_rng *rng, struct osh_particle_pool *pool, size_t n) {
    size_t i;
    size_t slot;
    struct ray_v ray;
    enum osh_status rc;

    /*
     * Scalar fill loop: one primary per iteration.
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
        rc = osh_beam_new_primary(rt->workspace, rng, &ray);
        if (rc != OSH_OK) {
            return rc;
        }

        slot = pool->n + i;
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
        pool->prim_idx[slot] = (uint32_t) (rt->primaries_generated + i);
        pool->gen[slot] = 0u;
        pool->species[slot] = &rt->primary;
    }

    pool->n += n;
    rt->primaries_generated += n;
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

static enum osh_status
fill_from_phsp(struct osh_beam_runtime *rt, struct osh_rng *rng, struct osh_particle_pool *pool, size_t n) {
    /*
     * PHSP (MCPL) fill — not yet implemented.
     *
     * When implemented, this function should:
     *   1. Read min(n, total - file_pos) entries from the MCPL file.
     *   2. Convert each MCPL particle record to pool SoA layout.
     *   3. Advance rt->source.phsp.file_pos; rewind when exhausted.
     *   4. Assign prim_idx from rt->primaries_generated and gen = 0.
     *
     * Key design constraint: MCPL files can exceed several GB, so entries
     * must be read in chunks — never all at once.  The pool capacity acts
     * as the chunk size; a single fill_pool call reads at most
     * pool->capacity entries from disk.
     *
     * TODO: implement once libmcpl integration is in place.
     */
    (void) rt;
    (void) rng;
    (void) pool;
    (void) n;
    return OSH_ENOTSUP;
}
