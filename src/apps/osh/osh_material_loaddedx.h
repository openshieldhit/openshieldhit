#ifndef OSH_APP_OSH_MATERIAL_LOADDEDX_H
#define OSH_APP_OSH_MATERIAL_LOADDEDX_H

#include <stddef.h>

#include "openshieldhit/diag.h"
#include "openshieldhit/status.h"

/**
 * @brief Legacy SHIELD-HIT style LOADDEDX table loaded from text.
 *
 * @details
 * This is an app-layer import helper for the `LOADDEDX` card parser. The file
 * format stores one shared energy-per-nucleon grid and one stopping-power
 * column per projectile:
 *
 *   energy  col_1 ... col_N
 *
 * Columns are interpreted as a contiguous ion sequence starting at `Z = 1`.
 *
 * This helper is intentionally private to the `src/apps/osh/` parser layer.
 * Core material setup consumes only material-owned in-memory dE/dx overrides.
 *
 * The loaded table is stored in projectile-major order:
 *   `[projectile_idx][energy_idx]` with energy as the innermost index.
 */
struct osh_material_loaddedx_table {
    double *energy_grid;        /* [nenergy], kinetic energy per nucleon [MeV/u]. */
    float *mass_stopping_power; /* [nprojectiles][nenergy], [MeV cm^2/g]. */
    unsigned int *projectile_z; /* [nprojectiles], contiguous Z = 1..nprojectiles. */
    size_t nprojectiles;
    size_t nenergy;
};

/**
 * @brief Load a LOADDEDX text table from @p path.
 *
 * @details
 * The file contains rows of whitespace-separated numbers. Each row holds one
 * kinetic-energy-per-nucleon value [MeV/u] followed by one stopping-power
 * column per projectile [MeV cm^2/g]. At least 18 projectile columns are
 * required and the energy grid must be strictly increasing.
 *
 * On success the table owns its allocated arrays; release with
 * @ref osh_material_loaddedx_table_free.
 *
 * @param[in]  path   Path to the LOADDEDX text file.
 * @param[out] table  Receives the loaded table. Must point to a
 *                    zero-initialised struct (e.g. after @c memset).
 *
 * @returns OSH_OK on success, OSH_EIO on file errors, OSH_EPARSE on format
 *          errors, OSH_ENOMEM on allocation failure.
 */
enum osh_status osh_material_loaddedx_table_load(char const *path,
                                                 struct osh_diag_sink const *diag,
                                                 struct osh_material_loaddedx_table *table);

/**
 * @brief Release all arrays owned by @p table and zero-reset the struct.
 *
 * @param[in,out] table  Table to free. Safe to call with NULL.
 */
void osh_material_loaddedx_table_free(struct osh_material_loaddedx_table *table);

/**
 * @brief Read the projectile atomic number Z at @p projectile_idx.
 *
 * @param[in]  table          Loaded LOADDEDX table.
 * @param[in]  projectile_idx Zero-based projectile index (< table->nprojectiles).
 * @param[out] z_out          Receives the atomic number Z (>= 1). May be NULL.
 *
 * @returns OSH_OK on success, OSH_EINVAL if @p table is NULL or
 *          @p projectile_idx is out of range.
 */
enum osh_status osh_material_loaddedx_projectile_z(struct osh_material_loaddedx_table const *table,
                                                   size_t projectile_idx,
                                                   unsigned int *z_out);

/**
 * @brief Compute the flat array offset for a (projectile, energy) pair.
 *
 * @param[in] table          Loaded LOADDEDX table.
 * @param[in] projectile_idx Zero-based projectile index.
 * @param[in] energy_idx     Zero-based energy index.
 *
 * @returns Flat index into @c table->mass_stopping_power.
 */
static inline size_t
osh_material_loaddedx_index(struct osh_material_loaddedx_table const *table, size_t projectile_idx, size_t energy_idx) {
    return projectile_idx * table->nenergy + energy_idx;
}

/**
 * @brief Return the mass stopping power [MeV cm^2/g] for a (projectile, energy) pair.
 *
 * @param[in] table          Loaded LOADDEDX table.
 * @param[in] projectile_idx Zero-based projectile index.
 * @param[in] energy_idx     Zero-based energy index.
 *
 * @returns Mass stopping power in MeV cm^2/g.
 */
static inline float osh_material_loaddedx_mass_stopping_power(struct osh_material_loaddedx_table const *table,
                                                              size_t projectile_idx,
                                                              size_t energy_idx) {
    return table->mass_stopping_power[osh_material_loaddedx_index(table, projectile_idx, energy_idx)];
}

#endif /* OSH_APP_OSH_MATERIAL_LOADDEDX_H */
