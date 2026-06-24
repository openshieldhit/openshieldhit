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
 * arrays.  No species array is needed — all entries are neutrons.
 *
 * The pool is an append-only bag drained by the neutron transport kernel: new
 * secondaries are appended and dead entries (e == 0) are removed by
 * osh_neutron_pool_compact().  Call osh_neutron_pool_reset() before starting a
 * new transport run to clear any residual entries.
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
 * @brief Initialise an embedded neutron pool with the given capacity.
 *
 * Allocates the internal SoA slab but not the struct itself.  Destroy
 * with osh_neutron_pool_free().
 */
enum osh_status osh_neutron_pool_init(struct osh_neutron_pool *pool, size_t capacity);

/**
 * @brief Free the slab owned by an embedded neutron pool (in-place destructor).
 *
 * Frees the SoA slab but does NOT free the pool struct itself.  Safe to call
 * with NULL.  Zeroes all pointer fields after freeing.
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
