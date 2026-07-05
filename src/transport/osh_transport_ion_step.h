#ifndef OSH_TRANSPORT_ION_STEP_H
#define OSH_TRANSPORT_ION_STEP_H

/*
 * Internal header: declares the per-slot ion step function used by the ion
 * transport loop.  Not part of the public API.
 */

#include <stddef.h> /* size_t */

#include "common/osh_step_segment.h"
#include "openshieldhit/geometry.h"
#include "openshieldhit/status.h"

struct osh_particle_pool;
struct osh_gemca_runtime;
struct osh_transport_context;
struct osh_transport_profile;
struct osh_material_runtime;
struct osh_scoring_runtime;
struct osh_score_target;
struct osh_rng;

/**
 * @brief Advance one ion pool slot by a single CSDA transport step.
 *
 * @details
 * Reads the particle state from @p pool[slot] and executes one transport step:
 * CSDA energy loss (using rho × ds for the current medium), optional
 * Highland/Molière MCS via random hinge, optional Bohr straggling, and scoring.
 * Writes the updated state back to the pool.
 *
 * The step_segments array currently contains exactly one current-medium segment:
 * analytic zones use the zone boundary distance and CT voxel zones use the
 * current voxel-exit distance.  Physics limits (CSDA, theta) may shorten the
 * step before that boundary is reached.
 *
 * The random-hinge treatment follows Fippel & Soukup (Med Phys 2004).
 *
 * @par References
 * Fippel M, Soukup M. A Monte Carlo dose calculation algorithm for proton
 * therapy. Med Phys. 2004;31(8):2263-2273. doi:10.1118/1.1769631.
 *
 * On physics termination (energy cutoff, geometry exit, blackhole)
 * pool->e[slot] is zeroed and OSH_OK is returned; the slot is collected by
 * the next osh_particle_pool_compact() call.  Only allocation or I/O errors
 * from scoring propagate as non-OK status.
 *
 * @p prof is the caller's (per-worker) profile for the nuclear-event counters;
 * pass NULL to disable, or the worker's own profile so concurrent workers never
 * race on shared counters.  It is read-only diagnostics: it never touches the
 * RNG or physics state.
 *
 * @p target is the caller-owned deposit target (accumulator set + traversal
 * scratch) that scoring writes into; a NULL @p target (or NULL fields within it)
 * falls back to @p score_rt's shared master views, reproducing the single-worker
 * behaviour bit-for-bit.  A replica or parallel worker passes its own private set.
 */
enum osh_status osh_transport_ion_step(struct osh_particle_pool *pool,
                                       size_t slot,
                                       struct osh_zone_ref const *zone_ref,
                                       struct osh_step_segment const *step_segments,
                                       size_t n_step_segments,
                                       struct osh_gemca_runtime const *geom_rt,
                                       struct osh_transport_context *transport_ctx,
                                       struct osh_transport_profile *prof,
                                       struct osh_material_runtime const *material_rt,
                                       struct osh_scoring_runtime *score_rt,
                                       struct osh_score_target const *target,
                                       struct osh_rng *rng);

#endif /* OSH_TRANSPORT_ION_STEP_H */
