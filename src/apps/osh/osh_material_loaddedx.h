#ifndef OSH_APP_OSH_MATERIAL_LOADDEDX_H
#define OSH_APP_OSH_MATERIAL_LOADDEDX_H

#include <stddef.h>

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
 */
struct osh_material_loaddedx_table {
    double *energy_grid;        /* [nenergy], kinetic energy per nucleon [MeV/u]. */
    float *mass_stopping_power; /* [nprojectiles][nenergy], [MeV cm^2/g]. */
    unsigned int *projectile_z; /* [nprojectiles], contiguous Z = 1..nprojectiles. */
    size_t nprojectiles;
    size_t nenergy;
};

enum osh_status osh_material_loaddedx_table_load(char const *path, struct osh_material_loaddedx_table *table);
void osh_material_loaddedx_table_free(struct osh_material_loaddedx_table *table);
enum osh_status osh_material_loaddedx_projectile_z(struct osh_material_loaddedx_table const *table,
                                                   size_t projectile_idx,
                                                   unsigned int *z_out);

static inline size_t
osh_material_loaddedx_index(struct osh_material_loaddedx_table const *table, size_t projectile_idx, size_t energy_idx) {
    return projectile_idx * table->nenergy + energy_idx;
}

static inline float osh_material_loaddedx_mass_stopping_power(struct osh_material_loaddedx_table const *table,
                                                              size_t projectile_idx,
                                                              size_t energy_idx) {
    return table->mass_stopping_power[osh_material_loaddedx_index(table, projectile_idx, energy_idx)];
}

#endif /* OSH_APP_OSH_MATERIAL_LOADDEDX_H */
