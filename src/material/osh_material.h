#ifndef OSH_MATERIAL_H
#define OSH_MATERIAL_H

#include <stddef.h>

#include "common/osh_logger.h"
#include "common/osh_rc.h"
#include "material/osh_material_icru.h" /* enum osh_material_state */

#ifdef __cplusplus
extern "C" {
#endif

enum { OSH_MATERIAL_INDEX_BLACKHOLE = 0, OSH_MATERIAL_INDEX_VACUUM = 1, OSH_MATERIAL_INDEX_FIRST_USER = 2 };

#define OSH_MATERIAL_NAME_BLACKHOLE "blackhole"
#define OSH_MATERIAL_NAME_VACUUM "vacuum"

struct material_element {
    double atom_count;             /* Relative atom count, < 0 until derived for mass-fraction input. */
    double mass_fraction;          /* Mass fraction, < 0 until derived for atom-count input. */
    double mean_excitation_energy; /* Element-level mean excitation override [eV], < 0 if unset. */
    size_t lineno;                 /* Input line where this element was defined. */
    unsigned int z;                /* Atomic number Z parsed from NUCLID/ELEMENT cards. */
    unsigned int a;                /* Mass number A; 0 means natural element, >0 means explicit isotope. */
};

struct material {
    struct material_element *elements;
    char *name;
    char *dedx_table_path; /* External dE/dx table path, owned; NULL if unset. */

    double rho; /* Density [g/cm^3], < 0 if unset. */
    /* Material-level mean excitation energy [eV], < 0 if unset; fallback for uncovered dE/dx. */
    double mean_excitation_energy;
    float rgba[4]; /* Rendering color [0,1]: red, green, blue, alpha. */

    size_t nelements;
    size_t lineno; /* Input line where this material was defined. */
    size_t index;  /* Dense internal index: 0 blackhole, 1 vacuum, >=2 user-defined. */

    int icru_id; /* ICRU id, 0 if unset. */
    int state;   /* enum osh_material_state value. */
};

struct material_workspace {
    struct material *materials;
    char *wdir;
    char *fname;
    size_t nmaterials;
};

/**
 * @brief Allocate, parse, and validate a material workspace from a mat.dat path.
 *
 * @details
 * This loader parses raw material definitions and performs the setup-time
 * assembly needed before transport-table generation: ICRU convenience
 * definitions are expanded, scalar user overrides are preserved, and explicit
 * element compositions are completed with mass fractions or relative atom
 * counts. It does not calculate stopping-power tables or derive optical depths.
 * Relative paths used by future material cards are resolved against the parsed
 * file's directory, which is retained in material_workspace::wdir.
 *
 * User materials are defined either by explicit element composition cards or by
 * an ICRU material reference. ICRU/default tables fill unset scalar properties,
 * while explicit user overrides such as RHO, STATE, MATERIALI/MIVALUE/MIAV, and
 * LOADDEDX remain authoritative.
 *
 * Mean excitation energy (MEE) precedence is:
 *  - single-element material: material-level and element-level MEE are
 *    equivalent; a known material MEE is copied to the sole element
 *  - compound with explicit element MEE values only: element MEE values are
 *    authoritative, and material MEE is derived by Bragg additivity
 *  - compound with only material-level MEE: material MEE is authoritative;
 *    element MEE defaults are filled from ICRU data and uniformly rescaled so
 *    their Bragg average matches the known material MEE
 *  - compound with both material-level MEE and explicit element-level MEE:
 *    invalid input and rejected during raw validation
 *
 * @param[in]  path    Path to the material input file.
 * @param[in]  lg      Logger for diagnostics; NULL uses the global default logger.
 * @param[out] wm_out  Receives the allocated workspace on success.
 *
 * @returns OSH_OK on success, or an OSH_E* code on failure.
 */
enum osh_status
osh_material_setup_from_path(char const *path, struct osh_logger *lg, struct material_workspace **wm_out);

/**
 * @brief Release a material workspace and all owned resources.
 *
 * @param[in] wm  Workspace to release. Safe to call with NULL.
 *
 * @returns OSH_OK.
 */
enum osh_status osh_material_workspace_free(struct material_workspace *wm);

/**
 * @brief Look up a material by internal dense index.
 *
 * @param[in] wm     Material workspace to search.
 * @param[in] index  Internal dense material index.
 *
 * @returns Pointer to the matching material, or NULL if no match exists.
 */
struct material const *osh_material_by_index(struct material_workspace const *wm, size_t index);

/**
 * @brief Look up a material by user-facing MATERIAL name.
 *
 * @param[in] wm    Material workspace to search.
 * @param[in] name  Material name from mat.dat.
 *
 * @returns Pointer to the matching material, or NULL if no match exists.
 */
struct material const *osh_material_by_name(struct material_workspace const *wm, char const *name);

/**
 * @brief Print a concise material workspace summary through the logger.
 *
 * @param[in] wm  Workspace to print.
 */
void osh_material_print(struct material_workspace const *wm);

#ifdef __cplusplus
}
#endif

#endif /* OSH_MATERIAL_H */
