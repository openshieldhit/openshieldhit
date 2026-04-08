#ifndef OSH_MATERIAL_RUNTIME_H
#define OSH_MATERIAL_RUNTIME_H

#include <stddef.h>

/**
 * @brief Runtime layout sketch for material-dependent atomic transport tables.
 *
 * @details
 * This is intentionally still a data-shape placeholder, not yet a constructed
 * API. The parsed/assembled material workspace is setup data; transport should
 * eventually use dense row-major arrays indexed by:
 *
 *   material_idx, projectile_idx, energy_idx
 *
 * The storage order is therefore:
 *
 *   [material][projectile][energy]
 *
 * with energy as the innermost index. Interpolation in stopping-power and
 * range tables is then a local contiguous read even when particles in a batch
 * drift apart in energy after the first transport step.
 *
 * Tables are mass-normalized where possible so voxel/CT geometries can apply a
 * local density per step without rebuilding material tables for every density.
 *
 * Large runtime tables are expected to use float storage for cache and GPU
 * friendliness. Table construction, Bethe evaluation, and range integration
 * may still use double internally before results are stored. The shared energy
 * grid is kept as double for now because interpolation coordinates span the
 * full transport energy range.
 *
 * This struct should own mutable storage while tables are being built. If we
 * want a read-only transport view later, that should be a separate const view
 * type rather than making the owning pointers const here.
 */
struct osh_material_runtime_tables {
    double *energy_grid;           /* [nenergy], kinetic energy per nucleon [MeV/nucleon]. */
    float *default_rho;            /* [nmaterials], reference density [g/cm^3]. */
    float *mean_excitation_energy; /* [nmaterials], material-level MEE [eV]. */
    float *mass_stopping_power;    /* [nmaterials][nprojectiles][nenergy], [MeV cm^2/g]. */
    float *range_mass_thickness;   /* [nmaterials][nprojectiles][nenergy], [g/cm^2]. */

    size_t nmaterials;
    size_t nprojectiles;
    size_t nenergy;
};

#endif /* OSH_MATERIAL_RUNTIME_H */
