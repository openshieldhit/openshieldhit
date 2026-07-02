#ifndef OSH_PARTICLE_POOL_H
#define OSH_PARTICLE_POOL_H

#include <stddef.h>
#include <stdint.h>

#include "openshieldhit/status.h"
#include "particle/osh_particle.h"

#ifdef __cplusplus
extern "C" {
#endif

struct osh_rng; /* per-slot RNG state; full definition in random/osh_rng.h */

/*
 * Particle pool — the shared interface between beam generation and transport.
 *
 * Design note — current placement:
 * The pool struct lives in common/ so that beam/runtime/ (which fills it) and
 * transport/ (which drains it) can both include it without creating a
 * dependency cycle.  This mirrors how struct step lives in common/osh_step.h
 * and is used by both transport and scoring.
 *
 * Design note — future redesign:
 * Once multiple pool types exist (ion pool, neutron pool, photon pool), it may
 * be worth promoting the pool into its own top-level module src/pool/ alongside
 * beam/ and transport/.  That would make the pool a first-class simulation
 * primitive with its own parse/runtime split, thread-local variants, and
 * GPU-friendly memory management.  For now, common/ is the right home.
 */

/**
 * @brief A pool of live particle histories in SoA (structure-of-arrays) layout.
 *
 * @details
 * The pool is the central data structure for batched particle transport.
 * It holds N particle histories simultaneously, with each field stored as a
 * contiguous flat array of length @p capacity.  Only the first @p n entries
 * are valid at any time; entries beyond n are uninitialized scratch space.
 *
 * SoA layout is chosen over AoS (array of struct ray_v) so that arithmetic
 * kernels operating on one field across all N particles (e.g. advancing all
 * x-coordinates by one step) access contiguous memory and can be
 * auto-vectorized by the compiler without gather/scatter instructions.
 *
 * Ownership and lifetime:
 *   - For embedded use (pool struct owned by the caller): initialise with
 *     osh_particle_pool_init() and destroy with osh_particle_pool_free().
 *   - For heap-allocated use: osh_particle_pool_alloc() allocates both the
 *     struct and the slab.  Release with osh_particle_pool_free(pool) then
 *     free(pool) (or keep _alloc/_free paired in the same scope).
 *   - The species[] pointers point into the particle registry (beam_workspace
 *     or a future shared registry); the pool does not own the species structs.
 *   - The pool does not own any geometry or material state.
 *
 * Filling and draining:
 *   - osh_beam_runtime_fill_pool() (beam/runtime/) appends fresh primaries.
 *   - The transport engine processes entries one step at a time and marks
 *     dead particles by zeroing their energy (e[i] <= 0).
 *   - osh_particle_pool_compact() removes dead entries, closing gaps so the
 *     next fill always appends to a dense prefix.
 *
 * Thread safety:
 *   - The pool is not thread-safe.  Multi-threaded transport uses one pool
 *     per thread (thread-local pools), with atomic scoring tallies.
 *
 * Batch size guidance:
 *   - capacity == NSTAT: all primaries live simultaneously; maximum
 *     parallelism, natural for GPU.
 *   - capacity == 4096..65536: pool fits in L2/L3 cache; optimal for CPU
 *     SIMD.  The outer loop runs ceil(NSTAT / capacity) times.
 *   - capacity == 1: degenerates to scalar serial transport; useful as a
 *     reference implementation for validation.
 */
struct osh_particle_pool {
    /* ---- Phase space (SoA) ---- */
    double *x;  /* position x [cm]             — capacity doubles */
    double *y;  /* position y [cm]             — capacity doubles */
    double *z;  /* position z [cm]             — capacity doubles */
    double *ux; /* direction unit vector x     — capacity doubles */
    double *uy; /* direction unit vector y     — capacity doubles */
    double *uz; /* direction unit vector z     — capacity doubles */
    double *e;  /* total kinetic energy [MeV]  — capacity doubles */

    /* ---- Per-history metadata (SoA) ---- */
    double *wt;         /* statistical weight; 1.0 = unweighted          — capacity doubles   */
    uint64_t *prim_idx; /* 0-based index of beam-primary ancestor        — capacity uint64_t  */
    uint8_t *gen;       /* generation: 0 = beam primary, 1 = secondary … — capacity uint8_t   */

    /* ---- Per-slot RNG state (SoA) ---- */
    /* One independent physics stream per live history, carried with the slot
     * so that random draws follow the particle rather than the wavefront
     * schedule.  A history thus sees the same draws on any pool capacity,
     * thread, or rank (scored output is invariant up to floating-point
     * summation order in the scorer, not byte-for-byte).  Primaries are seeded
     * by global history index at fill time; secondaries are split from their
     * parent. */
    struct osh_rng *rng; /* one RNG state per entry — capacity osh_rng */

    /* ---- Species (pointer array, NOT owned) ---- */
    struct particle const **species; /* one pointer per entry into the particle registry */

    /* ---- Bookkeeping ---- */
    size_t n;         /* number of live entries in [0, n)   */
    size_t capacity;  /* total allocated slots               */
    size_t n_dropped; /* secondaries lost to overflow (diag); mirrors the
                       * neutron pool's n_dropped.  Never reset by compaction,
                       * so it accumulates across wavefront passes and
                       * checkpoint batches for the run-level diagnostic. */
};

/* ---- Lifecycle ----------------------------------------------------------- */

/**
 * @brief Initialise a caller-allocated particle pool.
 *
 * @details
 * Allocates the SoA slab and wires all internal pointers.  Use this when the
 * pool struct itself is embedded in another object (e.g. osh_simulation).
 * Destroy with osh_particle_pool_free().
 *
 * @param[in,out] pool      Caller-allocated struct to initialise (must not be NULL).
 * @param[in]     capacity  Number of particle slots to allocate.
 *
 * @returns OSH_OK on success, OSH_ENOMEM on allocation failure,
 *          OSH_EINVAL if pool is NULL or capacity is zero.
 */
enum osh_status osh_particle_pool_init(struct osh_particle_pool *pool, size_t capacity);

/**
 * @brief Allocate a particle pool and its slab on the heap.
 *
 * @details
 * Allocates both the pool struct and the SoA slab.  The pool is
 * initialised with n = 0.  Release with osh_particle_pool_free(*pool_out)
 * followed by free(*pool_out).
 *
 * @param[in]  capacity  Number of particle slots to allocate.
 * @param[out] pool_out  Receives the allocated pool pointer on success.
 *
 * @returns OSH_OK on success, OSH_ENOMEM on allocation failure,
 *          OSH_EINVAL if capacity is zero or pool_out is NULL.
 */
enum osh_status osh_particle_pool_alloc(size_t capacity, struct osh_particle_pool **pool_out);

/**
 * @brief Free the slab owned by a particle pool (in-place destructor).
 *
 * @details
 * Frees the SoA slab but does NOT free the pool struct itself.  Use this
 * for pools initialised with osh_particle_pool_init() (embedded lifetime).
 * Safe to call with NULL.
 *
 * @param[in] pool  Pool whose slab should be freed; may be NULL.
 */
void osh_particle_pool_free(struct osh_particle_pool *pool);

/* ---- Operations ---------------------------------------------------------- */

/**
 * @brief Remove dead entries from the pool, compacting live entries to [0, n).
 *
 * @details
 * An entry is considered dead when e[i] <= 0.  Live entries are shifted
 * down to fill gaps; order within live entries is preserved.  On return,
 * n reflects the number of live entries remaining.
 *
 * This is an O(capacity) in-place scan — call once per transport round,
 * not once per step.
 *
 * @param[in,out] pool  Pool to compact.
 */
void osh_particle_pool_compact(struct osh_particle_pool *pool);

#ifdef __cplusplus
}
#endif

#endif /* OSH_PARTICLE_POOL_H */
