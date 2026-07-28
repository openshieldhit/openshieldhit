#ifndef OSH_BEAM_RUNTIME_H
#define OSH_BEAM_RUNTIME_H

#include <stddef.h>

#include "beam/osh_beam.h"
#include "common/osh_particle_pool.h"
#include "openshieldhit/status.h"
#include "particle/osh_particle.h"
#include "random/osh_rng.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Hot runtime state for beam primary generation.
 *
 * @details
 * Built once from a cold @ref beam_workspace and then used to fill
 * @ref osh_particle_pool instances during the transport loop.  The runtime
 * does not own the workspace; the workspace must outlive the runtime.
 *
 * Two source types are supported, distinguished by workspace->beam_mode:
 *
 * OSH_BEAM_MODE_SPOTS / OSH_BEAM_MODE_SOBP — spot-list sampling.
 *   Both single-spot and SOBP beams use the same code path: spot selection
 *   is weighted inverse-CDF sampling over workspace->spots[].  The spots
 *   union member is currently a placeholder; a future optimisation will
 *   pre-flatten per-spot scalars into SoA arrays here so that the fill loop
 *   can be SIMD-vectorized without accessing the AoS beam_spot structs.
 *
 * OSH_BEAM_MODE_PHSP — phase-space file streaming (MCPL).
 *   MCPL files can be 300 MB – several GB, so the full file is never loaded
 *   into memory.  Instead the phsp union member holds the open file handle
 *   and a read cursor; osh_beam_runtime_fill_pool() reads one chunk at a
 *   time.  This mode is not yet implemented — fill_pool returns
 *   OSH_ENOTSUP until MCPL support is added.
 *
 * primaries_generated is serial-path bookkeeping: it counts the primaries
 * emitted through the convenience wrapper osh_beam_runtime_fill_pool() so that
 * successive chunked fills get contiguous, globally-unique prim_idx values.
 * The cursor-free osh_beam_runtime_fill_pool_at() neither reads nor writes it —
 * a parallel driver supplies each worker's global base explicitly instead, so
 * the shared cursor is never a point of contention.
 */
struct osh_beam_runtime {
    struct beam_workspace const *workspace; /* cold storage — not owned */
    struct particle primary;                /* resolved species for beam primaries */
    uint64_t primaries_generated;           /* serial-path cursor; see header doc */

    union {
        /**
         * Spot-list source state (SPOTS and SOBP modes).
         *
         * Currently holds no extra data beyond the workspace reference.
         * Reserved for future SoA flattening of per-spot scalars:
         *   double *t0, *tsigma          — energies per spot
         *   double *size, *div, *cor     — phase-space parameters per spot
         *   double *tm                   — affine matrices (16 doubles × nspots)
         * Once flattened, fill_pool can drive SIMD loops over spots without
         * chasing pointers into the AoS beam_spot array.
         */
        struct {
            int _reserved; /* placeholder — no extra state yet */
        } spots;

        /**
         * Phase-space file streaming state (PHSP mode, MCPL format).
         *
         * mcpl_handle is an opaque pointer to the open MCPL file object.
         * NULL means the file has not been opened yet; setup opens it.
         * file_pos is the index of the next entry to read (0-based).
         * total is the number of entries in the file.
         *
         * fill_pool reads min(n, total - file_pos) entries per call.
         * When file_pos reaches total the file is rewound and reading
         * resumes from the start (recycling — consistent with how
         * SHIELD-HIT12A handled PHSP sources shorter than NSTAT).
         *
         * Not yet implemented: fill_pool returns OSH_ENOTSUP.
         * TODO: link against libmcpl and implement chunk reading.
         */
        struct {
            void *mcpl_handle; /* opaque libmcpl handle; NULL = not open */
            size_t file_pos;   /* next entry index to read */
            size_t total;      /* total entries in the file */
        } phsp;
    } source;
};

/* ---- Lifecycle ----------------------------------------------------------- */

/**
 * @brief Allocate and initialise a beam runtime from cold storage.
 *
 * @details
 * For SPOTS/SOBP mode: validates the workspace and sets up the spots state.
 * For PHSP mode:       opens the MCPL file and reads the entry count.
 *                      Currently returns OSH_ENOTSUP (not yet implemented).
 *
 * The workspace must remain valid and unmodified for the lifetime of the
 * returned runtime.
 *
 * @param[in]  workspace  Fully initialised cold beam workspace.
 * @param[out] rt_out     Receives the allocated runtime on success.
 *
 * @returns OSH_OK on success, OSH_E* on failure.
 */
enum osh_status osh_beam_compile(struct beam_workspace const *workspace, struct osh_beam_runtime **rt_out);

/**
 * @brief Free a beam runtime.
 *
 * @details
 * Releases the runtime struct.  For PHSP mode the MCPL file handle is
 * closed.  The workspace is not freed (it is not owned by the runtime).
 * Safe to call with NULL (or *rt NULL).  Sets *rt to NULL after freeing.
 *
 * @param[in,out] rt  Pointer to the runtime pointer; *rt is set to NULL.
 */
void osh_beam_runtime_free(struct osh_beam_runtime **rt);

/* ---- Primary generation -------------------------------------------------- */

/**
 * @brief Fill @p n new primary histories into the pool at an explicit global base.
 *
 * @details
 * The cursor-free core of primary generation.  Appends @p n entries to @p pool
 * starting at pool->n, then increments pool->n by @p n.  The caller must ensure
 * pool->n + n <= pool->capacity before calling; OSH_EINVAL is returned otherwise.
 *
 * The global history index of the i-th primary in this fill is
 * @p global_prim_base + i — supplied by the caller, not read from any internal
 * cursor.  This is the property that makes primary generation safe to partition
 * across workers: each worker fills its own slice of the run by passing its own
 * base (e.g. wctx->hist_lo + primaries_done), and because each primary is seeded
 * purely from its global index — gen is set to 0 (beam primary), prim_idx to the
 * global index, a transient BEAM stream samples the source phase space and the
 * slot's persistent PHYSICS stream (pool->rng[slot]) is initialised — the same
 * history draws the same numbers no matter which worker, thread, or rank emits
 * it, or in what order.  Concatenating the primaries of any contiguous partition
 * of [0, N) reproduces the single-pass [0, N) sequence exactly.
 *
 * @ref osh_beam_runtime::primaries_generated is neither read nor written here;
 * the caller owns the global index space.
 *
 * For PHSP mode: not yet implemented; returns OSH_ENOTSUP.
 *
 * @param[in,out] rt               Beam runtime (read-only here; cursor untouched).
 * @param[in]     seeding          Run-wide RNG seeding context (engine, seed).
 * @param[in,out] pool             Pool to fill; must have capacity >= pool->n + n.
 * @param[in]     n                Number of primaries to generate.
 * @param[in]     global_prim_base Global history index of the first primary.
 *
 * @returns OSH_OK on success, OSH_E* on failure.
 */
enum osh_status osh_beam_runtime_fill_pool_at(struct osh_beam_runtime *rt,
                                              struct osh_rng_seeding const *seeding,
                                              struct osh_particle_pool *pool,
                                              size_t n,
                                              uint64_t global_prim_base);

/**
 * @brief Fill @p n new primary histories into the pool (serial convenience).
 *
 * @details
 * Thin wrapper over @ref osh_beam_runtime_fill_pool_at that uses the runtime's
 * own @ref osh_beam_runtime::primaries_generated cursor as the global base and
 * advances it by @p n on success.  This is the right call for a single-stream
 * serial driver that fills one runtime in order; for partitioned/parallel work
 * call @ref osh_beam_runtime_fill_pool_at with an explicit base instead, so the
 * shared cursor is never read or mutated concurrently.
 *
 * Behaviour is otherwise identical to @ref osh_beam_runtime_fill_pool_at.
 *
 * @param[in,out] rt       Beam runtime (primaries_generated is advanced).
 * @param[in]     seeding  Run-wide RNG seeding context (engine, seed).
 * @param[in,out] pool     Pool to fill; must have capacity >= pool->n + n.
 * @param[in]     n        Number of primaries to generate.
 *
 * @returns OSH_OK on success, OSH_E* on failure.
 */
enum osh_status osh_beam_runtime_fill_pool(struct osh_beam_runtime *rt,
                                           struct osh_rng_seeding const *seeding,
                                           struct osh_particle_pool *pool,
                                           size_t n);

#ifdef __cplusplus
}
#endif

#endif /* OSH_BEAM_RUNTIME_H */
