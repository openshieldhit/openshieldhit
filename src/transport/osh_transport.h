#ifndef OSH_TRANSPORT_H
#define OSH_TRANSPORT_H

#include <stddef.h> /* size_t */

#include "common/osh_coord.h"
#include "common/osh_ray.h"
#include "common/osh_rc.h"
#include "common/osh_step.h"
#include "gemca/runtime/osh_gemca_runtime.h"

#ifdef __cplusplus
extern "C" {
#endif

struct osh_beam_runtime;
struct osh_material_runtime;
struct osh_scoring_runtime;

/*
 * struct step and struct position are defined in common/osh_step.h so that
 * domain modules (scoring, gemca) can receive them without depending on
 * transport/.  This header re-exports them via the include above.
 */

enum osh_transport_mcs_mode {
    OSH_TRANSPORT_MCS_OFF = 0,
    OSH_TRANSPORT_MCS_GAUSSIAN = 1,
    OSH_TRANSPORT_MCS_MOLIERE = 2
};

enum osh_transport_straggling_mode {
    OSH_TRANSPORT_STRAGGLING_OFF = 0,
    OSH_TRANSPORT_STRAGGLING_GAUSSIAN = 1,
    OSH_TRANSPORT_STRAGGLING_VAVILOV = 2
};

/**
 * @brief Immutable transport control parameters extracted from the beam configuration.
 *
 * @details
 * This struct carries the scalar knobs that govern how ions are transported.
 * It is populated by the caller before the transport call, so that the
 * transport layer never needs to include beam headers or touch the cold input
 * workspace directly.
 *
 * Transport model selections are stored as transport-owned enums rather than
 * beam-specific values, keeping the transport/runtime boundary explicit while
 * preserving enough information to reject unsupported models.
 *
 * Units follow the same conventions as the parsed input configuration:
 *   - energies in MeV or MeV/nucleon as noted
 *   - deltae and demin are per-step fractions / absolute floors
 */
struct osh_transport_params {
    size_t nstat;         /**< Total number of primary histories to transport. */
    float deltae;         /**< Max fractional energy loss per CSDA substep [dimensionless]. */
    float demin;          /**< Min energy loss per material substep [MeV/nucleon]. */
    float tcut;           /**< Lower ion energy cutoff [MeV/nucleon]. */
    int rndseed;          /**< Base RNG seed. */
    int rndoffset;        /**< RNG seed offset (added on top of rndseed). */
    char mcs_mode;        /**< enum osh_transport_mcs_mode value. */
    char straggling_mode; /**< enum osh_transport_straggling_mode value. */
};

/**
 * @brief Per-run transport context: immutable knobs plus mutable run state.
 *
 * @details
 * The params member is configured once before transport starts.  The mutable
 * fields hold warning-once bits and other run-lifetime state that must not be
 * stored in function-local statics if transport is later threaded.
 */
struct osh_transport_context {
    struct osh_transport_params params;
    char warned_boundary_demin_override;
};

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
enum osh_status osh_transport_run_minimal(struct osh_transport_context *transport_ctx,
                                          struct osh_beam_runtime *beam_rt,
                                          struct osh_gemca_runtime const *geom_rt,
                                          struct osh_material_runtime const *tables,
                                          struct osh_scoring_runtime *scoring);

#ifdef __cplusplus
}
#endif

#endif /* OSH_TRANSPORT_H */
