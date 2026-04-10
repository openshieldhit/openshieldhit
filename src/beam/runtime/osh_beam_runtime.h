#ifndef OSH_BEAM_RUNTIME_H
#define OSH_BEAM_RUNTIME_H

#include <stddef.h>

#include "beam/osh_beam.h"
#include "common/osh_particle_pool.h"
#include "common/osh_rc.h"
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
 * primaries_generated tracks the total number of primary histories emitted
 * across all fill_pool calls on this runtime.  It is used to assign globally
 * unique prim_idx values when the transport loop chunks a large beam into
 * multiple pool fills.
 */
struct osh_beam_runtime {
    struct beam_workspace const *workspace; /* cold storage — not owned */
    size_t primaries_generated;             /* total primaries emitted so far */

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
enum osh_status osh_beam_runtime_setup(struct beam_workspace const *workspace, struct osh_beam_runtime **rt_out);

/**
 * @brief Free a beam runtime.
 *
 * @details
 * Releases the runtime struct.  For PHSP mode the MCPL file handle is
 * closed.  The workspace is not freed (it is not owned by the runtime).
 * Safe to call with NULL.
 *
 * @param[in] rt  Runtime to free; may be NULL.
 */
void osh_beam_runtime_free(struct osh_beam_runtime *rt);

/* ---- Primary generation -------------------------------------------------- */

/**
 * @brief Fill @p n new primary histories into the pool.
 *
 * @details
 * Appends @p n entries to @p pool starting at pool->n, then increments
 * pool->n by @p n.  The caller must ensure pool->n + n <= pool->capacity
 * before calling; OSH_EINVAL is returned otherwise.
 *
 * For SPOTS/SOBP mode: samples spot, energy, and phase space for each
 * primary using the beam model in beam/osh_beam_model.c.  Each entry is
 * written directly into the pool SoA arrays, bypassing the intermediate
 * AoS ray_v representation.  gen is set to 0 (beam primary) and prim_idx
 * is assigned from rt->primaries_generated so that indices are unique and
 * contiguous across multiple fill calls on the same runtime.
 *
 * For PHSP mode: not yet implemented; returns OSH_ENOTSUP.
 *
 * @param[in,out] rt    Beam runtime (primaries_generated is updated).
 * @param[in,out] rng   Random-number generator state.
 * @param[in,out] pool  Pool to fill; must have capacity >= pool->n + n.
 * @param[in]     n     Number of primaries to generate.
 *
 * @returns OSH_OK on success, OSH_E* on failure.
 */
enum osh_status
osh_beam_runtime_fill_pool(struct osh_beam_runtime *rt, struct osh_rng *rng, struct osh_particle_pool *pool, size_t n);

#ifdef __cplusplus
}
#endif

#endif /* OSH_BEAM_RUNTIME_H */
