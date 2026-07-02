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
    OSH_TRANSPORT_MCS_GAUSSIAN = 1, /**< Highland Gaussian core (no tail). */
    OSH_TRANSPORT_MCS_MOLIERE = 2,  /**< Full Bethe-Molière (core + Rutherford tail). */
    OSH_TRANSPORT_MCS_WENTZEL = 3   /**< Reserved: Geant4 Wentzel-VI/Urban screened model (not implemented). */
};

enum osh_transport_straggling_mode {
    OSH_TRANSPORT_STRAGGLING_OFF = 0,
    OSH_TRANSPORT_STRAGGLING_GAUSSIAN = 1,
    OSH_TRANSPORT_STRAGGLING_VAVILOV = 2,
    OSH_TRANSPORT_STRAGGLING_URBAN = 3 /* Reserved (G4 Urbán); not implemented. */
};

/* Overridable from the build line (e.g. -DOSH_TRANSPORT_POOL_CAPACITY=256)
 * so the benchmark harness can sweep capacities without editing sources. */
#ifndef OSH_TRANSPORT_POOL_CAPACITY
#define OSH_TRANSPORT_POOL_CAPACITY 4096u
#endif

/** Default lower neutron transport cutoff [MeV] used when NEUTRLCUT <= 0. */
#define OSH_TRANSPORT_NEUTRON_CUTOFF_DEFAULT_MEV 1.0e-3f

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
    size_t nstat;                      /**< Total number of primary histories to transport. */
    size_t pool_capacity;              /**< Live-history pool size; 0 selects the compiled default.
                                            Tunes the cache/parallelism trade-off only: per-history
                                            RNG streams make the physics each history sees identical
                                            across capacities (scored output matches up to
                                            floating-point reduction order in scoring). */
    float deltae;                      /**< Max fractional energy loss per CSDA substep [dimensionless]. */
    float demin;                       /**< Min energy loss per material substep [MeV/nucleon]. */
    float tcut;                        /**< Lower ion energy cutoff [MeV/nucleon]. */
    float ncut;                        /**< Lower neutron energy cutoff [MeV]; <=0 uses transport default. */
    int rndseed;                       /**< Base RNG seed (RNDSEED); fixes the whole run's streams. */
    int rndoffset;                     /**< Global history-index base (RNDOFFSET): added to every
                                            history index before seeding, so disjoint ranges (e.g.
                                            one per process/MPI rank) give disjoint, non-overlapping
                                            streams.  Not a value added to rndseed. */
    char mcs_mode;                     /**< enum osh_transport_mcs_mode value. */
    char straggling_mode;              /**< enum osh_transport_straggling_mode value. */
    char nuclear_inelastic;            /**< Non-zero to enable inelastic nuclear reactions (Tripathi). */
    char nuclear_elastic;              /**< Non-zero to enable pp elastic scattering. */
    char nuclear_neutron_ion_feedback; /**< Non-zero to reschedule the ion family after the
                                            neutron pass, enabling a coupled ion ↔ neutron loop.
                                            0 (default): one ion pass + one neutron pass (SH12A
                                            mode); scheduler exits after the neutron pass.
                                            1: scheduler marks the ion family as having work
                                            again once the ion pool holds entries from (n,x)
                                            reactions, repeating until both pools empty.
                                            Not yet active: osh_transport_neutron.c does not
                                            push to the ion pool yet, so even with flag=1 the
                                            ion family is never rescheduled. */
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
 *
 * A profile is per-worker mutable state: each transport worker accumulates into
 * its own instance on the hot path (carried on @ref osh_worker_context), never
 * into a shared one, so concurrent workers do not race on these counters.  The
 * driver folds the per-worker profiles into the run's master profile after the
 * workers finish, with osh_transport_profile_merge().  In the single-worker
 * build the lone worker writes straight into the master, so there is nothing to
 * merge and the values are identical to the un-parallelised path.
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
 * @brief Fold one worker's transport profile into a destination profile.
 *
 * @details
 * The reduce step that combines per-worker profiles after a parallel run,
 * mirroring osh_scoring_accumulator_merge() for scoring tallies.  Counters and
 * per-phase timers are summed: the phase @c *_s fields therefore report
 * aggregate work-in-phase across workers (their sum can exceed the elapsed wall
 * time, which is the point — it exposes parallel efficiency).
 *
 * @c total_s is the exception.  Workers run concurrently, so elapsed time is the
 * span of the longest worker, not the sum of their spans; it is combined by
 * **maximum**, not addition.  A parallel driver that measures a single
 * wall-clock span around the whole parallel region may overwrite @c total_s with
 * that span afterwards; until then the max of the per-worker spans is a faithful
 * lower-bound proxy.  Merging one worker into a zeroed destination reproduces
 * that worker's values exactly (sum with zero; max with zero), so the
 * single-worker path is bit-identical whether or not the merge is used.
 *
 * @param[in,out] dst  Destination profile (accumulates @p src).  No-op if NULL.
 * @param[in]     src  Source profile to fold in.  No-op if NULL.
 */
void osh_transport_profile_merge(struct osh_transport_profile *dst, struct osh_transport_profile const *src);

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
struct osh_particle_pool;
struct osh_run_control;

struct osh_transport_context {
    struct osh_transport_params params;
    struct osh_diag_sink const *diag;                  /**< Borrowed; NULL silences transport diagnostics. */
    struct osh_nuclear_handler const *nuclear_handler; /**< Borrowed; NULL disables handler. */
    struct osh_fragment_pool *fragment_pool;           /**< Borrowed; residual fragments for future breakup. */
    struct osh_neutron_pool *neutron_pool;             /**< Borrowed; neutrons routed here instead of CSDA pool. */
    struct osh_particle_pool *ion_pool;                /**< Borrowed; pre-allocated ion wavefront pool. */
    struct osh_zone_ref *zone_refs;                    /**< Transport scratch: zone/material per slot. */
    double *dist_batch;                                /**< Transport scratch: boundary distance per slot. */
    size_t scratch_capacity;                           /**< Number of entries in zone_refs and dist_batch. */
    struct osh_transport_profile *profile;             /**< Borrowed master/destination profile; NULL disables
                                                            phase timers/counters.  Written by a worker (the lone
                                                            serial worker points its own profile here) or by the
                                                            driver merging per-worker profiles — not on the hot
                                                            path through this pointer. */
    struct osh_run_control const *run_control;         /**< Borrowed; clean-stop / wall-budget policy (issue #192).
                                                            NULL = run to completion, no early stop.  Consulted only
                                                            by the primary (ion) loop to halt new-primary injection;
                                                            the family scheduler still drains every banked secondary
                                                            so a partial result stays family-exact (issue #195). */
    size_t completed_primaries;                        /**< Out: primaries fully transported by the last run.  Equals
                                                            params.nstat unless a clean stop drained the run early; the
                                                            driver copies it into completed_nstat so output normalises
                                                            by the true count. */
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
 * Internally, transport uses a scheduler-owned outer loop with one pool per
 * particle family (ions, neutrons, photons, electrons, ...).  The dispatch
 * seam lives in transport/ rather than inside each family kernel, so that
 * per-family kernels are pure functions of their pool and can be offloaded or
 * parallelised independently in future refactors.
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
