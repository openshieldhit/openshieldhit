#ifndef OPENSHIELDHIT_MATERIAL_H
#define OPENSHIELDHIT_MATERIAL_H

#include <stddef.h>

#include "openshieldhit/diag.h"
#include "openshieldhit/status.h"
#include "openshieldhit/voxel.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Physical state of a material.
 */
enum osh_material_state { OSH_MATERIAL_STATE_UNSET = 0, OSH_MATERIAL_STATE_CONDENSED = 1, OSH_MATERIAL_STATE_GAS = 2 };

enum { OSH_MATERIAL_INDEX_BLACKHOLE = 0, OSH_MATERIAL_INDEX_VACUUM = 1, OSH_MATERIAL_INDEX_FIRST_USER = 2 };

#define OSH_MATERIAL_NAME_BLACKHOLE "blackhole"
#define OSH_MATERIAL_NAME_VACUUM "vacuum"

struct osh_material_element {
    double atom_count;             /* Relative atom count, < 0 until derived for mass-fraction input. */
    double mass_fraction;          /* Mass fraction, < 0 until derived for atom-count input. */
    double mean_excitation_energy; /* Element-level mean excitation override [eV], < 0 if unset. */
    size_t lineno;                 /* Input line where this element was defined. */
    unsigned int z;                /* Atomic number Z parsed from NUCLID/ELEMENT cards. */
    unsigned int a;                /* Mass number A; 0 means natural element, >0 means explicit isotope. */
};

/**
 * @brief Tabulated stopping-power override stored inside one material.
 *
 * @details
 * The arrays referenced by this struct are owned by the containing
 * @ref osh_material and, transitively, by the owning
 * @ref osh_material_workspace. They remain valid only while that workspace
 * remains alive and unchanged by the library. Callers may read these arrays
 * but must not free, reallocate, or retain them beyond the owning workspace
 * lifetime.
 */
struct osh_material_dedx_override {
    double *energy_mev_per_u;   /* [npoints], strictly increasing kinetic energy per nucleon [MeV/u]. */
    double *dedx_mev_cm2_per_g; /* [npoints], mass stopping power [MeV cm^2/g]. */
    unsigned int projectile_z;  /* Projectile atomic number Z. */
    size_t npoints;
};

/**
 * @brief Material definition stored inside a material workspace.
 *
 * @details
 * All pointer members of this struct are owned by the containing
 * @ref osh_material_workspace. Their contents remain valid only while that
 * workspace remains alive and unchanged by the library. API consumers may
 * read these fields but must not free or reallocate them individually.
 */
struct osh_material {
    struct osh_material_element *elements;
    struct osh_material_dedx_override *dedx_overrides;
    char *name;

    double rho; /* Density [g/cm^3], < 0 if unset. */
    /* Material-level mean excitation energy [eV], < 0 if unset; fallback for uncovered dE/dx. */
    double mean_excitation_energy;
    float rgba[4]; /* Rendering color [0,1]: red, green, blue, alpha. */

    size_t nelements;
    size_t ndedx_overrides;
    size_t lineno; /* Input line where this material was defined. */
    size_t index;  /* Dense internal index: 0 blackhole, 1 vacuum, >=2 user-defined. */

    int icru_id; /* ICRU id, 0 if unset. */
    int state;   /* enum osh_material_state value. */
};

/**
 * @brief Heap-allocated container for all cold material definitions.
 *
 * @details
 * This workspace owns the storage referenced by its pointer members and by the
 * nested pointer members of each contained @ref osh_material. Callers must
 * release the workspace as a whole via @ref osh_material_workspace_free()
 * rather than freeing individual fields such as @ref materials, @ref wdir,
 * @ref fname, or nested arrays/strings.
 */
struct osh_material_workspace {
    struct osh_material *materials;
    char *wdir;
    char *fname;
    size_t nmaterials;
    int hu_table_type; /**< OSH_HU_TABLE_* selector (from openshieldhit/voxel.h); 0 (NONE) for non-CT runs. */
};

/**
 * @brief Allocate an empty material workspace and initialize reserved entries.
 *
 * @details
 * Initializes the built-in reserved materials:
 * - index 0: blackhole
 * - index 1: vacuum
 *
 * Application-layer parsers may then append user materials before calling
 * @ref osh_material_workspace_prepare().
 *
 * @param[out] wm_out  Receives the allocated workspace on success.
 *
 * @returns OSH_OK on success, or an OSH_E* code on failure.
 */
enum osh_status osh_material_workspace_create(struct osh_material_workspace **wm_out);

/**
 * @brief Validate and complete a parsed material workspace in place.
 *
 * @details
 * Runs the setup-time material assembly needed before transport-table
 * generation:
 * - validates raw parsed cards
 * - expands ICRU convenience definitions
 * - derives complementary composition fields
 * - resolves and validates mean excitation energies
 *
 * This function does not parse files and does not build transport tables.
 *
 * @param[in,out] wm    Workspace to validate and complete.
 * @param[in]     diag  Borrowed diagnostics sink for setup messages, or NULL.
 *
 * @returns OSH_OK on success, or an OSH_E* code on failure.
 */
enum osh_status osh_material_workspace_prepare(struct osh_material_workspace *wm, struct osh_diag_sink const *diag);

/**
 * @brief Set or replace one material-owned dE/dx override curve.
 *
 * @details
 * The curve is keyed by projectile atomic number `Z` and is deep-copied into
 * the material. Repeated calls with the same projectile replace the existing
 * curve. The runtime material compile step uses these overrides where present
 * and falls back to Bethe elsewhere.
 *
 * @param[in,out] mat                  Material to update.
 * @param[in]     z_projectile         Projectile atomic number `Z` (> 0).
 * @param[in]     energy_mev_per_u     Energy grid [MeV/u], length `n_points`.
 * @param[in]     dedx_mev_cm2_per_g   Mass stopping powers [MeV cm^2/g], length `n_points`.
 * @param[in]     n_points             Number of samples, must be >= 2.
 *
 * @returns OSH_OK on success, or an error code on invalid input or allocation
 *          failure.
 */
enum osh_status osh_material_dedx_set(struct osh_material *mat,
                                      unsigned int z_projectile,
                                      double const *energy_mev_per_u,
                                      double const *dedx_mev_cm2_per_g,
                                      size_t n_points);

/**
 * @brief Remove all material-owned dE/dx overrides from one material.
 *
 * @param[in,out] mat  Material to clear. Safe to call with NULL.
 */
void osh_material_dedx_clear(struct osh_material *mat);

/**
 * @brief Release a material workspace and all owned resources.
 *
 * @param[in] wm  Workspace to release. Safe to call with NULL.
 *
 * @returns OSH_OK.
 */
enum osh_status osh_material_workspace_free(struct osh_material_workspace *wm);

/**
 * @brief Look up a material by internal dense index.
 *
 * @param[in] wm     Material workspace to search.
 * @param[in] index  Internal dense material index.
 *
 * @returns Pointer to the matching material, or NULL if no match exists.
 */
struct osh_material const *osh_material_by_index(struct osh_material_workspace const *wm, size_t index);

/**
 * @brief Look up a material by user-facing MATERIAL name.
 *
 * @param[in] wm    Material workspace to search.
 * @param[in] name  Material name from mat.dat.
 *
 * @returns Pointer to the matching material, or NULL if no match exists.
 */
struct osh_material const *osh_material_by_name(struct osh_material_workspace const *wm, char const *name);

/**
 * @brief Print a concise material workspace summary through a diagnostics sink.
 *
 * @param[in] wm    Workspace to print.
 * @param[in] diag  Borrowed diagnostics sink for summary output, or NULL.
 */
void osh_material_print(struct osh_material_workspace const *wm, struct osh_diag_sink const *diag);

#ifdef __cplusplus
}
#endif

#endif /* OPENSHIELDHIT_MATERIAL_H */
