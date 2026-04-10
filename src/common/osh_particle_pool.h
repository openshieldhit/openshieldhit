#ifndef OSH_PARTICLE_POOL_H
#define OSH_PARTICLE_POOL_H

#include <stddef.h>
#include <stdint.h>

#include "common/osh_rc.h"
#include "particle/osh_particle.h"

#ifdef __cplusplus
extern "C" {
#endif

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
 *   - All arrays are allocated in a single contiguous slab by
 *     osh_particle_pool_alloc() and freed by osh_particle_pool_free().
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
    uint32_t *prim_idx; /* 0-based index of beam-primary ancestor        — capacity uint32_t  */
    uint8_t *gen;       /* generation: 0 = beam primary, 1 = secondary … — capacity uint8_t   */

    /* ---- Species (pointer array, NOT owned) ---- */
    struct particle const **species; /* one pointer per entry into the particle registry */

    /* ---- Bookkeeping ---- */
    size_t n;        /* number of live entries in [0, n)   */
    size_t capacity; /* total allocated slots               */
};

/* ---- Lifecycle ----------------------------------------------------------- */

/**
 * @brief Allocate a particle pool with the given capacity.
 *
 * @details
 * All SoA arrays are allocated in a single contiguous slab and partitioned
 * internally, so a single free() releases all memory.  The pool is
 * initialised with n = 0.
 *
 * @param[in]  capacity  Number of particle slots to allocate.
 * @param[out] pool_out  Receives the allocated pool pointer on success.
 *
 * @returns OSH_OK on success, OSH_ENOMEM on allocation failure,
 *          OSH_EINVAL if capacity is zero or pool_out is NULL.
 */
enum osh_status osh_particle_pool_alloc(size_t capacity, struct osh_particle_pool **pool_out);

/**
 * @brief Free a particle pool and all its arrays.
 *
 * @details
 * Safe to call with NULL.  The species pointers are not freed (they are not
 * owned by the pool).
 *
 * @param[in] pool  Pool to free; may be NULL.
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
