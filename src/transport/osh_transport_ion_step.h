#ifndef OSH_TRANSPORT_ION_STEP_H
#define OSH_TRANSPORT_ION_STEP_H

/*
 * Internal header: declares the per-slot ion step function used by the ion
 * transport loop.  Not part of the public API.
 */

#include <stddef.h> /* size_t */

#include "common/osh_step_segment.h"
#include "openshieldhit/status.h"

struct osh_particle_pool;
struct osh_gemca_runtime;
struct osh_transport_context;
struct osh_material_runtime;
struct osh_scoring_runtime;
struct osh_rng;

/**
 * @brief Advance one ion pool slot by a single CSDA transport step.
 *
 * @details
 * Reads the particle state from @p pool[slot] and executes one transport step:
 * CSDA energy loss (using the per-segment areal density Σ(rho_i × ds_i)),
 * optional Highland/Molière MCS via random hinge, optional Bohr straggling,
 * and per-segment scoring.  Writes the updated state back to the pool.
 *
 * The step_segments array describes the constant-density pieces the step may
 * traverse.  For analytic (non-voxelised) zones the caller supplies exactly one
 * synthetic segment {ds = boundary_ds, rho = zone_rho}.  For CT voxel zones the
 * caller supplies the Jacobs traversal output from dist_voxel_body_rt(), which
 * may contain up to OSH_STEP_SEGMENTS_MAX segments.  Physics limits (CSDA,
 * theta) may shorten the step so that fewer than n_step_segments are actually
 * traversed; segments beyond the actual step length are ignored.
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
 */
enum osh_status osh_transport_ion_step(struct osh_particle_pool *pool,
                                       size_t slot,
                                       size_t zone_idx,
                                       struct osh_step_segment const *step_segments,
                                       size_t n_step_segments,
                                       struct osh_gemca_runtime const *geom_rt,
                                       struct osh_transport_context *transport_ctx,
                                       struct osh_material_runtime const *material_rt,
                                       struct osh_scoring_runtime *score_rt,
                                       struct osh_rng *rng);

#endif /* OSH_TRANSPORT_ION_STEP_H */
