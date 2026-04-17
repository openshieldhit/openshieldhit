#ifndef OSH_MATERIAL_LOADDEDX_H
#define OSH_MATERIAL_LOADDEDX_H

#include <stddef.h>

#include "openshieldhit/status.h"

/**
 * @brief External legacy dE/dx table loaded from a LOADDEDX text file.
 *
 * @details
 * The legacy/libdEdx text files used by SHIELD-HIT style inputs store one
 * shared energy-per-nucleon grid and one stopping-power column for each ion.
 * Comments may vary between files, but the numeric payload follows a fixed
 * pattern:
 *
 *   energy  col_1 ... col_N
 *
 * where energy is in [MeV/nucleon] and stopping powers are in [MeV cm^2/g].
 * Columns are interpreted as a continuous ion sequence starting at Z = 1. The
 * representative isotope A for each column is taken from the default isotope
 * list in particle/.
 *
 * This is source material data attached to a material definition, not the final
 * transport runtime table. Runtime builders may resample it onto a separate
 * canonical energy grid later.
 *
 * The loaded table is stored in projectile-major order:
 *
 *   [projectile_idx][energy_idx]
 *
 * with energy as the innermost index.
 */
struct osh_material_loaddedx_table {
    double *energy_grid;        /* [nenergy], kinetic energy per nucleon [MeV/nucleon]. */
    float *mass_stopping_power; /* [nprojectiles][nenergy], [MeV cm^2/g]. */
    unsigned int *projectile_z; /* [nprojectiles], contiguous Z = 1..nprojectiles. */
    unsigned int *projectile_a; /* [nprojectiles], default isotope mass numbers. */
    size_t nprojectiles;
    size_t nenergy;
};

enum osh_status osh_material_loaddedx_table_load(char const *path, struct osh_material_loaddedx_table *table);
void osh_material_loaddedx_table_free(struct osh_material_loaddedx_table *table);
enum osh_status osh_material_loaddedx_projectile_za(struct osh_material_loaddedx_table const *table,
                                                    size_t projectile_idx,
                                                    unsigned int *z_out,
                                                    unsigned int *a_out);

static inline size_t
osh_material_loaddedx_index(struct osh_material_loaddedx_table const *table, size_t projectile_idx, size_t energy_idx) {
    return projectile_idx * table->nenergy + energy_idx;
}

static inline float osh_material_loaddedx_mass_stopping_power(struct osh_material_loaddedx_table const *table,
                                                              size_t projectile_idx,
                                                              size_t energy_idx) {
    return table->mass_stopping_power[osh_material_loaddedx_index(table, projectile_idx, energy_idx)];
}

#endif /* OSH_MATERIAL_LOADDEDX_H */
