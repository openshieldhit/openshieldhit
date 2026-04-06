#ifndef OSH_MATERIAL_RUNTIME_H
#define OSH_MATERIAL_RUNTIME_H

#include <stddef.h>

/**
 * @brief Initial runtime layout sketch for material-dependent transport tables.
 *
 * @details
 * This is intentionally a data-shape placeholder, not yet a constructed API.
 * The parsed/assembled material workspace is user-facing setup data; transport
 * should eventually use dense cache-friendly arrays indexed by:
 *
 *   material_idx, projectile_idx, energy_idx
 *
 * Tables are mass-normalized where possible so voxel/CT geometries can apply a
 * local density per step without rebuilding material tables for every density.
 */
struct osh_material_runtime_tables {
    double *energy_grid; /* Kinetic-energy grid [MeV or MeV/n], length nenergy. */

    float *default_rho;            /* [nmaterials], reference density [g/cm^3]. */
    float *mean_excitation_energy; /* [nmaterials], material-level MEE [eV]. */
    float *mass_stopping_power;    /* [nmaterials][nprojectiles][nenergy], MeV cm^2/g. */
    float *range_mass_thickness;   /* [nmaterials][nprojectiles][nenergy], g/cm^2. */
    float *interaction_coeff_mass; /* [nmaterials][nprojectiles][nenergy], cm^2/g. */

    size_t nmaterials;
    size_t nprojectiles;
    size_t nenergy;
};

#endif /* OSH_MATERIAL_RUNTIME_H */
