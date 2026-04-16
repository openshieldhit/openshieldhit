#ifndef OSH_TRANSPORT_ION_STEP_H
#define OSH_TRANSPORT_ION_STEP_H

/*
 * Internal header: declares the per-slot ion step function used by the ion
 * transport loop.  Not part of the public API.
 */

#include <stddef.h> /* size_t */

#include "common/osh_rc.h"

struct osh_particle_pool;
struct gemca_runtime;
struct beam_workspace;
struct osh_material_runtime;
struct osh_scoring_runtime;
struct osh_rng;

/**
 * @brief Advance one ion pool slot by a single CSDA transport substep.
 *
 * @details
 * Reads the particle state from @p pool[slot], executes one substep (CSDA
 * energy loss, optional Highland/Molière MCS via random hinge, optional Bohr
 * straggling), scores it, and writes the updated state back.
 *
 * The random-hinge treatment follows the approach described by Fippel and
 * Soukup for fast proton transport.
 *
 * @par References
 * Fippel M, Soukup M. A Monte Carlo dose calculation algorithm for proton
 * therapy. Med Phys. 2004;31(8):2263-2273. doi:10.1118/1.1769631.
 *
 * On physics termination (energy cutoff, geometry exit, blackhole)
 * pool->e[slot] is zeroed and OSH_OK is returned; the slot is collected by
 * the next osh_particle_pool_compact() call.  Only allocation or I/O errors
 * from scoring propagate as non-OK status.
 */
enum osh_status osh_transport_ion_step_one(struct osh_particle_pool *pool,
                                           size_t slot,
                                           size_t zone_idx,
                                           double boundary_ds,
                                           struct gemca_runtime const *geom_rt,
                                           struct beam_workspace const *beam,
                                           struct osh_material_runtime const *tables,
                                           struct osh_scoring_runtime *scoring,
                                           double deltae,
                                           struct osh_rng *rng);

#endif /* OSH_TRANSPORT_ION_STEP_H */
