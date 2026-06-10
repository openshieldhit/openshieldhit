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
    int pdg;
    unsigned int z;
    unsigned int a;
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
    double p[3];
    double size[2];
    double div[2];
    double cor[2];
    double t0;
    double tsigma;
    double p0;
    double psigma;
    double wt;
    unsigned int spot_id;
    unsigned int layer_id;
    char shape;
    char tsigma_type;
    char t0_per_nucleon;
    char tsigma_per_nucleon;
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
    double sad[2];
    double focus;
    double theta;
    double phi;
    char use_div;
    char use_sad;
    char sad_was_set;
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
    double *p[3];
    double *d[3];
    double *e;
    double *wt;
    int32_t *pdg;
    size_t len;
    char *fname;
    int32_t _cached_pdg;
    struct osh_beam_species _cached_species;
    char _cached_species_valid;
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
    struct osh_beam_phsp *phsp;
    struct osh_beam_spot *spots;
    struct osh_beam_species primary;
    struct osh_beam_shared shared;
    struct osh_beam_prepared *prepared; /* internal prepared state; owned by core */
    struct ripple_filter *rifi;
    struct parlev *parlev;
    size_t nspots;
    char has_primary;
    size_t nstat;
    size_t nsave;
    int rndseed;
    int rndoffset;
    float tcut;
    float pcut;
    float ncut;
    float deltae;
    float demin;
    char straggl;
    char scatter;
    char nuclear_inelastic; /**< Non-zero to enable inelastic nuclear reactions. */
    char nuclear_elastic;   /**< Non-zero to enable pp elastic scattering. */
    char emtrans;
    char apcorr;
    char beam_mode;
    char makeln;
    char neutrfast;
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
