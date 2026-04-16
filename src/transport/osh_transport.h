#ifndef OSH_TRANSPORT_H
#define OSH_TRANSPORT_H

#include "common/osh_coord.h"
#include "common/osh_ray.h"
#include "common/osh_rc.h"
#include "common/osh_step.h"
#include "gemca/runtime/osh_gemca_runtime.h"

#ifdef __cplusplus
extern "C" {
#endif

struct beam_workspace;
struct osh_material_runtime;
struct osh_scoring_runtime;

/*
 * struct step and struct position are defined in common/osh_step.h so that
 * domain modules (scoring, gemca) can receive them without depending on
 * transport/.  This header re-exports them via the include above.
 */

/**
 * @brief Run the minimal straight-line CSDA transport loop.
 *
 * @details
 * This is the first end-to-end transport slice for OpenShieldHIT:
 *   - primaries are sampled from beam/
 *   - the top-level driver currently seeds the ion family only
 *   - zones and boundary distances come from gemca/runtime/
 *   - energy loss is computed from material CSDA range tables
 *   - scoring is applied step-by-step through scoring/runtime
 *
 * Internally, transport is being prepared for a scheduler-owned outer loop
 * with one queue or pool per particle family (ions, neutrons, photons,
 * electrons, ...).  The current implementation still executes only the ion
 * family, but the dispatch seam now lives in transport/ rather than inside
 * the ion kernel.
 *
 * Physics included:
 *   - CSDA energy loss (residual-range tables)
 *   - multiple Coulomb scattering (Highland/Molière, random-hinge method)
 *   - Gaussian energy straggling (Bohr variance)
 *
 * Not yet implemented:
 *   - nuclear interactions
 *   - secondaries
 *
 * DELTAE is treated as a maximum fractional energy-loss step criterion in
 * material. Boundary-limited steps are truncated and the exit energy is
 * recovered from the residual CSDA range.
 *
 * The caller is responsible for calling osh_gemca_runtime_setup() before this
 * function and osh_gemca_runtime_free() after it returns.
 */
enum osh_status osh_transport_run_minimal(struct beam_workspace const *beam,
                                          struct gemca_runtime const *geom_rt,
                                          struct osh_material_runtime const *tables,
                                          struct osh_scoring_runtime *scoring);

#ifdef __cplusplus
}
#endif

#endif /* OSH_TRANSPORT_H */
