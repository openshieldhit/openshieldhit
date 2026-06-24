#ifndef OSH_NEUTRON_POOL_H
#define OSH_NEUTRON_POOL_H

#include <stddef.h>
#include <stdint.h>

#include "openshieldhit/status.h"

#ifdef __cplusplus
extern "C" {
#endif

struct osh_rng; /* full definition in random/osh_rng.h */

/**
 * @file osh_neutron_pool.h
 * @brief Neutron secondary pool in SoA layout.
 *
 * @details
 * Neutrons produced by nuclear reactions (abrasion, Fermi break-up, (n,2n))
 * are not transported through the CSDA charged-particle loop.  They are
 * collected here and drained by the neutron transport kernel.
 *
 * Memory layout follows osh_particle_pool: a single contiguous slab holds all
 * arrays, so osh_neutron_pool_free() is a single free().  No species array is
 * needed — all entries are neutrons.
 *
 * The pool is used as a stack: transport pops from slot n-1, new secondaries
 * are appended.  Call osh_neutron_pool_reset() between beam primaries to clear
 * the pool for the next history.
 */
struct osh_neutron_pool {
    /* SoA phase space */
    double *x;  /**< position x [cm]           */
    double *y;  /**< position y [cm]           */
    double *z;  /**< position z [cm]           */
    double *ux; /**< direction unit vector x   */
    double *uy; /**< direction unit vector y   */
    double *uz; /**< direction unit vector z   */
    double *e;  /**< kinetic energy [MeV]      */
    double *wt; /**< statistical weight        */

    /* Per-history metadata */
    uint64_t *prim_idx;  /**< beam-primary ancestor index              */
    uint8_t *gen;        /**< generation; saturates at 255             */
    struct osh_rng *rng; /**< per-slot RNG state (one stream/neutron)  */

    /* Bookkeeping */
    size_t n;         /**< live entries in [0, n)              */
    size_t capacity;  /**< allocated slots                     */
    size_t n_created; /**< cumulative neutrons accepted (diag) */
    size_t n_dropped; /**< neutrons lost to overflow   (diag) */
};

/**
 * @brief Allocate a neutron pool with the given capacity.
 *
 * All SoA arrays are carved from one malloc() call.  The pool is initialised
 * with n = n_created = n_dropped = 0.
 *
 * @param[in]  capacity  Number of neutron slots to allocate (must be > 0).
 * @param[out] out       Receives the allocated pool pointer on success.
 * @returns OSH_OK, OSH_ENOMEM, or OSH_EINVAL.
 */
enum osh_status osh_neutron_pool_alloc(size_t capacity, struct osh_neutron_pool **out);

/**
 * @brief Free a neutron pool and all its arrays.  Safe to call with NULL.
 */
void osh_neutron_pool_free(struct osh_neutron_pool *pool);

/**
 * @brief Reset the pool between beam primaries.
 *
 * Zeroes n, n_created, and n_dropped; does not free memory.
 */
void osh_neutron_pool_reset(struct osh_neutron_pool *pool);

/**
 * @brief Remove dead entries (e == 0) by in-place squeezing.
 *
 * Mirrors osh_particle_pool_compact().  Called at the end of each neutron
 * wavefront pass to remove killed neutrons before the next pass begins.
 */
void osh_neutron_pool_compact(struct osh_neutron_pool *pool);

#ifdef __cplusplus
}
#endif

#endif /* OSH_NEUTRON_POOL_H */
