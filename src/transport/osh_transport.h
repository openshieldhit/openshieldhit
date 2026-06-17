#ifndef OSH_TRANSPORT_H
#define OSH_TRANSPORT_H

#include <stddef.h> /* size_t */

#include "common/osh_coord.h"
#include "common/osh_ray.h"
#include "common/osh_step.h"
#include "gemca/runtime/osh_gemca_runtime.h"
#include "openshieldhit/diag.h"
#include "openshieldhit/status.h"

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
    size_t nstat;           /**< Total number of primary histories to transport. */
    size_t pool_capacity;   /**< Live-history pool size; 0 selects the compiled default.
                                 Tunes the cache/parallelism trade-off only: per-history
                                 RNG streams make the physics each history sees identical
                                 across capacities (scored output matches up to
                                 floating-point reduction order in scoring). */
    float deltae;           /**< Max fractional energy loss per CSDA substep [dimensionless]. */
    float demin;            /**< Min energy loss per material substep [MeV/nucleon]. */
    float tcut;             /**< Lower ion energy cutoff [MeV/nucleon]. */
    int rndseed;            /**< Base RNG seed (RNDSEED); fixes the whole run's streams. */
    int rndoffset;          /**< Global history-index base (RNDOFFSET): added to every
                                 history index before seeding, so disjoint ranges (e.g.
                                 one per process/MPI rank) give disjoint, non-overlapping
                                 streams.  Not a value added to rndseed. */
    char mcs_mode;          /**< enum osh_transport_mcs_mode value. */
    char straggling_mode;   /**< enum osh_transport_straggling_mode value. */
    char nuclear_inelastic; /**< Non-zero to enable inelastic nuclear reactions (Tripathi). */
    char nuclear_elastic;   /**< Non-zero to enable pp elastic scattering. */
};

/**
 * @brief Wall-clock phase timers and event counters for one transport run.
 *
 * @details
 * Filled by the wavefront loop when profiling is enabled (the context's
 * profile pointer is non-NULL).  The five phase timers decompose the
 * transport wall time:
 *
 *   fill_s      pool refill from the beam source
 *   zone_ref_s  batched zone-ref lookup (zone + current HU/material)
 *   distance_s  batched current-medium boundary distance query
 *   step_s      per-particle physics step loop (osh_transport_ion_step)
 *   compact_s   dead-slot compaction and progress bookkeeping
 *
 * Profiling reads only the monotonic clock and pre-existing counters; it
 * never touches the RNG streams or any physics state, so instrumented runs
 * are bit-identical to un-instrumented ones.
 */
struct osh_transport_profile {
    double fill_s;                     /**< Time spent refilling the pool [s]. */
    double zone_ref_s;                 /**< Time in the zone-ref batch query [s]. */
    double distance_s;                 /**< Time in the boundary-distance batch query [s]. */
    double step_s;                     /**< Time in the per-particle step loop [s]. */
    double compact_s;                  /**< Time in pool compaction [s]. */
    double total_s;                    /**< Total transport wall time [s]. */
    unsigned long long steps;          /**< Total transport steps taken. */
    unsigned long long iterations;     /**< Wavefront loop iterations. */
    unsigned long long nuclear_events; /**< Nuclear interactions sampled (any kind). */
    unsigned long long secondaries;    /**< Secondaries produced by nuclear events. */
};

/**
 * @brief Per-run transport context: immutable knobs plus mutable run state.
 *
 * @details
 * The params member is configured once before transport starts.  The mutable
 * fields hold warning-once bits and other run-lifetime state that must not be
 * stored in function-local statics if transport is later threaded.
 */
struct osh_nuclear_handler;
struct osh_fragment_pool;
struct osh_neutron_pool;

struct osh_transport_context {
    struct osh_transport_params params;
    struct osh_diag_sink const *diag;                  /**< Borrowed; NULL silences transport diagnostics. */
    struct osh_nuclear_handler const *nuclear_handler; /**< Borrowed; NULL disables handler. */
    struct osh_fragment_pool *fragment_pool;           /**< Borrowed; residual fragments for future breakup. */
    struct osh_neutron_pool *neutron_pool;             /**< Borrowed; neutrons routed here instead of CSDA pool. */
    struct osh_transport_profile *profile;             /**< Borrowed; NULL disables phase timers/counters. */
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
 *   - nuclear interactions (pp elastic, abrasion + Fermi break-up secondaries)
 *
 * Not yet implemented:
 *   - de-excitation of heavy nuclear residues, A > 16 (evaporation/SMM)
 *
 * DELTAE is treated as a maximum fractional energy-loss step criterion in
 * material. Boundary-limited steps are truncated and the exit energy is
 * recovered from the residual CSDA range.
 *
 * The caller is responsible for calling osh_gemca_compile() before this
 * function and osh_gemca_runtime_free() after it returns.
 */
enum osh_status osh_transport_run_minimal(struct osh_transport_context *transport_ctx,
                                          struct osh_beam_runtime *beam_rt,
                                          struct osh_gemca_runtime const *geom_rt,
                                          struct osh_material_runtime const *material_rt,
                                          struct osh_scoring_runtime *score_rt);

#ifdef __cplusplus
}
#endif

#endif /* OSH_TRANSPORT_H */
