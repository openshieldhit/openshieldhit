#ifndef OPENSHIELDHIT_BEAM_H
#define OPENSHIELDHIT_BEAM_H

#include <stddef.h>
#include <stdint.h>

#include "openshieldhit/beam_defs.h"
#include "openshieldhit/diag.h"
#include "openshieldhit/status.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @file beam.h
 * @brief Public cold-data beam model for `osh_core`.
 *
 * @details
 * This header defines the beam input structs that library users and frontend
 * adapters fill before calling @ref osh_beam_workspace_prepare().
 *
 * Design rule:
 * - Types here are the public cold model.
 * - Derived transport/setup state lives in core-private prepared/runtime
 *   structs and must not be written by callers.
 * - Frontends may populate these structs from files, JSON, WASM bindings, or
 *   any other input representation, but `osh_core` owns the prepare step.
 */

struct ripple_filter;
struct parlev;
struct osh_beam_prepared;

/**
 * @brief Public species identity used by the cold beam model.
 *
 * @details
 * This is intentionally lighter than `struct particle`: it carries only the
 * user-visible identity needed to describe the primary species. The full
 * particle constants are resolved later by core prepare/runtime code.
 */
struct osh_beam_species {
    int pdg;        /**< PDG Monte Carlo particle code; takes precedence over z/a when non-zero. */
    unsigned int z; /**< Atomic number (proton number). */
    unsigned int a; /**< Mass number (nucleon count). */
};

/**
 * @brief One cold beam spot description.
 *
 * @details
 * All values here are user/input-facing beam parameters. Prepared quantities
 * such as affine transforms or cumulative SOBP weights are stored outside the
 * public struct during @ref osh_beam_workspace_prepare().
 *
 * Coordinate convention for p[]:
 *   p[0], p[1] — physical beam-start position (x, y) in the beam-local PZALIGN
 *                frame [cm], measured from the beam axis at the beam-entrance
 *                plane (BEAMPOS z). These are NOT isocenter coordinates.
 *   p[2]       — beam-start z in PZALIGN [cm]; negative means upstream of the
 *                isocenter (isocenter is at z = 0 in PZALIGN).
 *
 * When loading from a USECBEAM spotlist file the app layer converts from the
 * file's isocenter convention to these physical beam-start coordinates before
 * populating this struct (see osh_app_osh.c).  The cold struct itself is
 * coordinate-agnostic: it always holds physical beam-start positions so that
 * the beam model and any other consumer can use it without knowing about
 * upstream file-format conventions.
 */
struct osh_beam_spot {
    double p[3];             /**< Physical beam-start position [cm]; see coordinate convention above. */
    double size[2];          /**< Lateral 1σ half-widths [cm] for x (0) and y (1).  Sign encodes
                                  profile: positive = Gaussian σ, negative = uniform half-width,
                                  zero = pencil beam.  Matches the BEAMSIGMA sign convention. */
    double div[2];           /**< Angular divergence [mrad] for x (0) and y (1); zero when unused. */
    double cor[2];           /**< Emittance correlation between position and divergence for x (0) and
                                  y (1); dimensionless, range [−1, 1]; zero when unused. */
    double t0;               /**< Nominal kinetic energy [MeV]; or [MeV/nucleon] when t0_per_nucleon
                                  is set.  Zero when p0 is used instead. */
    double tsigma;           /**< Energy spread 1σ [MeV or MeV/nucleon]; zero for monoenergetic beam. */
    double p0;               /**< Nominal momentum [MeV/c]; used when t0 is not set. */
    double psigma;           /**< Momentum spread 1σ [MeV/c]; zero when unused. */
    double wt;               /**< Relative spot weight (arbitrary units); normalised across all spots
                                  during prepare. */
    unsigned int spot_id;    /**< Spot identifier from the input file; informational only. */
    unsigned int layer_id;   /**< Energy layer identifier from the input file; informational only. */
    char shape;              /**< Lateral profile shape: 'G' = Gaussian, 'U' = uniform. */
    char tsigma_type;        /**< Origin of energy spread: 0 = from TMAX0/spot file,
                                  1 = from EXTSPEC table. */
    char t0_per_nucleon;     /**< Non-zero when t0 and tsigma are in MeV/nucleon (ions);
                                  converted to total MeV during prepare. */
    char tsigma_per_nucleon; /**< Non-zero when tsigma is per-nucleon; converted during prepare. */
};

/**
 * @brief Beam parameters shared by all spots.
 *
 * sad[]:       Source-to-Axis Distance [cm] for x (index 0) and y (index 1).
 *              Always a positive distance from the virtual upstream point source
 *              to the isocenter, regardless of beam direction.  Zero when SAD is
 *              not active.
 * use_sad:     1 if fan-out correction is active (finite SAD), 0 for parallel beam.
 * sad_was_set: 1 if BEAMSAD was explicitly present in the input file (even if
 *              set to INF).  Used to suppress the "USECBEAM without BEAMSAD"
 *              warning when the user intentionally chose parallel delivery.
 */
struct osh_beam_shared {
    double sad[2];    /**< Source-to-Axis Distance [cm] for x (0) and y (1); always positive.
                           Zero when SAD is not active.  See sad_was_set for intent. */
    double focus;     /**< Focus position along the beam axis [cm]; used with BEAMDIV.
                           Zero when angular divergence is not focused. */
    double theta;     /**< Beam direction polar angle [rad]; 0 = along +z. */
    double phi;       /**< Beam direction azimuthal angle [rad]. */
    char use_div;     /**< Non-zero when angular divergence (div[] in spots) is active. */
    char use_sad;     /**< 1 = finite SAD, fan-out correction applied; 0 = parallel beam. */
    char sad_was_set; /**< 1 if BEAMSAD was explicitly set in the input (even if INF).
                           Suppresses the "USECBEAM without BEAMSAD" warning for intentional
                           parallel delivery. */
};

/**
 * @brief Phase-space source descriptor.
 *
 * @details
 * For large external phase-space sources such as MCPL, the core runtime will
 * stream data in chunks later. The cold model therefore stores only the
 * source description and lightweight caches, not the full particle content.
 */
struct osh_beam_phsp {
    double *p[3];                            /**< Position arrays [cm]: p[0]=x, p[1]=y, p[2]=z; each length len. */
    double *d[3];                            /**< Direction cosine arrays: d[0]=u, d[1]=v, d[2]=w; each length len. */
    double *e;                               /**< Kinetic energy array [MeV]; length len. */
    double *wt;                              /**< Statistical weight array; length len. */
    int32_t *pdg;                            /**< PDG particle code array; length len. */
    size_t len;                              /**< Number of particles in this phase-space snapshot. */
    char *fname;                             /**< Source filename for streaming/reload; owned by this struct. */
    int32_t _cached_pdg;                     /**< Internal cache: PDG of the dominant species. */
    struct osh_beam_species _cached_species; /**< Internal cache: resolved species for _cached_pdg. */
    char _cached_species_valid;              /**< Non-zero when the cached species is up to date. */
};

/**
 * @brief Public cold beam workspace.
 *
 * @details
 * This is the input object populated by applications or frontend adapters.
 * The `prepared` pointer is reserved for core-internal prepared state and
 * must be treated as opaque by callers.
 */
struct osh_beam_workspace {
    struct osh_beam_phsp *phsp;         /**< Phase-space source; mutually exclusive with spots.
                                             NULL when not used. */
    struct osh_beam_spot *spots;        /**< Owned array of nspots cold spot descriptions. */
    struct osh_beam_species primary;    /**< Primary particle species shared by all spots. */
    struct osh_beam_shared shared;      /**< Beam parameters common to all spots (SAD, direction). */
    struct osh_beam_prepared *prepared; /**< Internal prepared state; owned by core, opaque to callers. */
    struct ripple_filter *rifi;         /**< Ripple filter descriptor; NULL when not used
                                             (feature not yet implemented). */
    struct parlev *parlev;              /**< Parallel lever optics descriptor; NULL when not used
                                             (feature not yet implemented). */
    size_t nspots;                      /**< Number of spots in spots[]; must be >= 1 when spots is set. */
    char has_primary;                   /**< Non-zero if primary species was explicitly set by the caller. */
    size_t nstat;                       /**< Total primary histories to simulate (NSTAT). */
    size_t nsave;                       /**< Save interval in histories; 0 = write only at end (NSTAT step). */
    double wall_budget_s;               /**< Wall-time budget [s] from the MAXTIME card; 0 = unlimited.
                                             A CLI --max-time overrides this value. */
    double dump_every_s;                /**< Periodic partial-result dump time cadence [s] from the DUMPEVERY
                                             card; 0 = off.  A CLI --dump-every overrides this value.  The
                                             count-cadence equivalent is the NSTAT save step (@ref nsave). */
    int rndseed;                        /**< RNG seed; same seed gives bit-for-bit reproducible results (RNDSEED). */
    int rndoffset;                      /**< RNG stream offset for producing independent parallel runs
                                             from the same seed. */
    float tcut;                         /**< Primary kinetic energy cutoff [MeV/nucleon]; primaries
                                             below this are killed during transport, and it doubles as
                                             the lower bound of the initial-sampling truncation window
                                             when tcut_upper is set (TCUT0). */
    float tcut_upper;                   /**< Upper bound [MeV/nucleon] of the initial-sampling truncation
                                             window (TCUT0 upper argument); 0 = unset, meaning the primary
                                             energy is drawn from the plain (untruncated) Gaussian N(t0,
                                             tsigma^2).  Unlike tcut, this does not affect in-flight
                                             transport — it only bounds the energy primaries are born with. */
    float pcut;                         /**< Secondary proton kinetic energy cutoff [MeV]. */
    float ncut;                         /**< Neutron energy cutoff [MeV]; <=0 uses the transport
                                             default (NEUTRLCUT). */
    float deltae;                       /**< Maximum relative energy loss per step, e.g. 0.005 = 0.5%
                                             (DELTAE). */
    float demin;                        /**< Minimum kinetic energy step size [MeV/nucleon]; prevents
                                             extremely small steps near the Bragg peak (DEMIN). */
    char straggl;                       /**< Energy straggling model: 0=off, 1=Gaussian, 2=Vavilov (STRAGG). */
    char scatter;                       /**< Multiple Coulomb scattering model: 0=off, 1=Gaussian/Rossi-Greisen,
                                             2=Molière (MSCAT). */
    char nuclear_inelastic;             /**< Non-zero to enable inelastic nuclear reactions (NUCRE 1 or 3). */
    char nuclear_elastic;               /**< Non-zero to enable pp nuclear elastic scattering (NUCRE 1 or 2). */
    char emtrans;                       /**< Non-zero to enable electromagnetic transport corrections. */
    char apcorr;                        /**< Non-zero to apply AP correction factor (APCORR). */
    char beam_mode;                     /**< Internal beam delivery mode flag. */
    char makeln;                        /**< Non-zero to write primary phase-space to file (MAKELN). */
    char neutrfast;                     /**< Non-zero for fast (analog) neutron transport. */
};

/**
 * @brief Allocate a beam workspace with core defaults.
 */
enum osh_status osh_beam_workspace_create(struct osh_beam_workspace **wb_out);

/**
 * @brief Build or rebuild internal prepared state from a cold beam workspace.
 *
 * @details
 * `osh_beam_workspace_create()` initializes `prepared` to NULL. Calling this
 * function resolves derived quantities and allocates the internal prepared
 * state. Repeated calls are allowed after cold-data edits; the previous
 * prepared state is released and rebuilt in place.
 *
 * @param[in,out] wb    Beam workspace to validate and finalize.
 * @param[in]     diag  Borrowed diagnostics sink for setup messages, or NULL.
 */
enum osh_status osh_beam_workspace_prepare(struct osh_beam_workspace *wb, struct osh_diag_sink const *diag);

/**
 * @brief Replace the owned beam spot array with caller-provided cold spot data.
 *
 * @details
 * The input array is deep-copied into the workspace. Existing workspace spots
 * are released and replaced on success. This function performs only lightweight
 * structural checks; heavier physics-derived validation remains part of
 * @ref osh_beam_workspace_prepare().
 *
 * Callers typically use this after parsing a beam source description from any
 * frontend format. Both single-spot beams and external spot lists should enter
 * the workspace through this function.
 *
 * @param[in,out] wb      Beam workspace to update.
 * @param[in]     spots   Spot array to copy, length @p nspots.
 * @param[in]     nspots  Number of spots; must be >= 1.
 *
 * @returns OSH_OK on success, or an error code on invalid input or allocation
 *          failure.
 */
enum osh_status osh_beam_spots_set(struct osh_beam_workspace *wb, struct osh_beam_spot const *spots, size_t nspots);

/**
 * @brief Free a beam workspace allocated by @ref osh_beam_workspace_create.
 */
enum osh_status osh_beam_workspace_free(struct osh_beam_workspace *wb);
enum osh_status osh_beam_phsp_free(struct osh_beam_phsp *phsp);

/**
 * @brief Print a concise beam workspace summary through a diagnostics sink.
 *
 * @param[in] wb    Workspace to print.
 * @param[in] diag  Borrowed diagnostics sink for summary output, or NULL.
 */
void osh_beam_print(struct osh_beam_workspace const *wb, struct osh_diag_sink const *diag);

/**
 * @brief Print one beam spot summary through a diagnostics sink.
 *
 * @param[in] spot  Spot to print.
 * @param[in] diag  Borrowed diagnostics sink for summary output, or NULL.
 */
void osh_beam_print_spot(struct osh_beam_spot const *spot, struct osh_diag_sink const *diag);

#ifdef __cplusplus
}
#endif

#endif /* OPENSHIELDHIT_BEAM_H */
