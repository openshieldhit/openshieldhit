#ifndef OSH_MATERIAL_RUNTIME_H
#define OSH_MATERIAL_RUNTIME_H

#include <assert.h>
#include <math.h>
#include <stddef.h>

#include "common/osh_hd.h"
#include "openshieldhit/geometry.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Material-owned runtime layout for material-dependent atomic tables.
 *
 * @details
 * This struct is the compiled CPU/GPU-facing boundary produced by the prepare
 * layer and consumed by the hot transport kernel. The parsed material workspace
 * remains setup data; transport uses only dense row-major arrays indexed by:
 *
 *   material_idx, projectile_idx, energy_idx
 *
 * Storage order is:
 *
 *   [material][projectile][energy]
 *
 * with energy innermost so interpolation reads stay contiguous even after a
 * particle batch drifts apart in energy. Large tables use float storage for
 * cache and GPU friendliness, while build-time evaluation and integration may
 * still use double internally before results are stored.
 *
 * Tables are mass-normalized so voxel/CT geometries can apply a local density
 * per step without rebuilding material tables for every density.
 */
struct osh_material_runtime {
    double log_emin;             /* log(energy_grid[0]) for O(1) index. */
    double inv_dlog;             /* 1 / (log(emax/emin) / (nenergy-1)). */
    double emin;                 /* Minimum energy [MeV/nucleon]. */
    double emax;                 /* Maximum energy [MeV/nucleon]. */
    float *mass_stopping_power;  /* [nmaterials][nprojectiles][nenergy] [MeV cm^2/g]. */
    float *range_csda;           /* [nmaterials][nprojectiles][nenergy] [g/cm^2]. */
    unsigned int *projectile_z;  /* [nprojectiles]. */
    unsigned int *projectile_a;  /* [nprojectiles]. */
    double *projectile_mass_mev; /* [nprojectiles] nuclear rest mass [MeV/c^2]. */

    /*
     * Per-material scalars used by multiple-scattering and energy-straggling
     * models.  These are derived from the element composition at prepare time
     * and stored here so the hot transport kernel does not need to touch the
     * material_workspace (which is setup-only data).
     *
     * rho         — material density [g/cm^3].
     *               Used to convert geometric step lengths to areal density in
     *               the hot transport kernel. Zero for vacuum and blackhole.
     * z_mean      — effective atomic number (mass-fraction-weighted sum of Z_i).
     *               Used by the Hubert effective-charge formula in osh_physics_bethe_z_eff().
     * z_over_a    — effective Z/A [mol/g] (sum_i w_i * Z_i / A_i).
     *               Used by the Bohr straggling variance in osh_physics_strag_sigma().
     * rad_length  — radiation length X₀ [g/cm²] (PDG approximate formula, mixture rule).
     *               Used by the Highland MCS formula in osh_physics_highland_theta0().
     *               Zero for vacuum and blackhole (skips scattering in transport).
     */
    float *rho;        /* [nmaterials] density [g/cm^3]. */
    float *z_mean;     /* [nmaterials] effective atomic number. */
    float *z_over_a;   /* [nmaterials] Z/A [mol/g]. */
    float *rad_length; /* [nmaterials] radiation length X0 [g/cm^2]. */
    /* Bethe-Molière (mode 2) per-medium constants; see compute_moliere_constants(). */
    float *moliere_chic2;    /* [nmaterials] χ_c² coefficient 0.157·Σ w_i Z_i(Z_i+1)/A_i [MeV² cm²/g];
                                χ_c² = moliere_chic2 · z_eff² · d / (pβ)². Zero for vacuum/blackhole. */
    float *moliere_screen_z; /* [nmaterials] effective screening Z for the Molière screening angle χ_a. */
    float *hu_rho_lut;       /* [2601] HU→density [g/cm³] indexed by hu+1000; NULL for non-CT runs. */

    size_t nmaterials;
    size_t nprojectiles;
    size_t nenergy;
};

/**
 * @brief Default number of energy grid points for the runtime table.
 *
 * 500 log-uniform points over [0.025, 1000] MeV/nucleon give ~1.4% relative
 * spacing, well within the ~4 significant digits of the source data.
 */
enum { OSH_MATERIAL_RUNTIME_NENERGY = 500 };

/** Default lower energy bound [MeV/nucleon], matching LOADDEDX file range. */
#define OSH_MATERIAL_RUNTIME_EMIN 0.025

/** Default upper energy bound [MeV/nucleon], matching LOADDEDX file range. */
#define OSH_MATERIAL_RUNTIME_EMAX 1000.0

/**
 * @brief Look up mass stopping power at an arbitrary energy.
 *
 * @details
 * O(1) index computation via log(E), then linear interpolation between
 * adjacent grid points. Energy is clamped to [emin, emax] rather than
 * extrapolating.
 *
 * @param[in] tables       Runtime tables built by the prepare layer.
 * @param[in] mat_idx      Dense material index (from zone->material_idx).
 * @param[in] proj_idx     Projectile column index (0-based, Z = proj_idx+1).
 * @param[in] e_per_nuc    Kinetic energy per nucleon [MeV/nucleon].
 *
 * @returns Mass stopping power [MeV cm^2/g].
 */
OSH_HD static inline double osh_material_runtime_sp_lookup(struct osh_material_runtime const *tables,
                                                    size_t mat_idx,
                                                    size_t proj_idx,
                                                    double e_per_nuc) {
    double x;
    size_t idx;
    double frac;
    size_t base;

    if (e_per_nuc <= tables->emin) {
        e_per_nuc = tables->emin;
    } else if (e_per_nuc >= tables->emax) {
        e_per_nuc = tables->emax;
    }

    x = (log(e_per_nuc) - tables->log_emin) * tables->inv_dlog;
    idx = (size_t) x;
    if (idx >= tables->nenergy - 1u) {
        idx = tables->nenergy - 2u;
        x = (double) idx + 1.0;
    }
    frac = x - (double) idx;
    base = (mat_idx * tables->nprojectiles + proj_idx) * tables->nenergy;
    return (double) tables->mass_stopping_power[base + idx] * (1.0 - frac)
           + (double) tables->mass_stopping_power[base + idx + 1u] * frac;
}

/**
 * @brief Look up CSDA range at an arbitrary energy.
 *
 * @details
 * Uses the same log-grid indexing and linear interpolation strategy as
 * osh_material_runtime_sp_lookup(). Energy is clamped to [emin, emax].
 *
 * @param[in] tables       Runtime tables built by the prepare layer.
 * @param[in] mat_idx      Dense material index.
 * @param[in] proj_idx     Projectile column index.
 * @param[in] e_per_nuc    Kinetic energy per nucleon [MeV/nucleon].
 *
 * @returns CSDA range [g/cm^2].
 *
 * @note
 * Transport currently calls this in a scalar fashion.  A natural follow-up is
 * a batch/SIMD sibling for lanes sharing the same material/projectile column,
 * so the ion stepper can amortize the paired forward and inverse range lookups.
 */
OSH_HD static inline double osh_material_runtime_range_lookup(struct osh_material_runtime const *tables,
                                                        size_t mat_idx,
                                                        size_t proj_idx,
                                                        double e_per_nuc) {
    double x;
    size_t idx;
    double frac;
    size_t base;

    if (e_per_nuc <= tables->emin) {
        e_per_nuc = tables->emin;
    } else if (e_per_nuc >= tables->emax) {
        e_per_nuc = tables->emax;
    }

    x = (log(e_per_nuc) - tables->log_emin) * tables->inv_dlog;
    idx = (size_t) x;
    if (idx >= tables->nenergy - 1u) {
        idx = tables->nenergy - 2u;
        x = (double) idx + 1.0;
    }
    frac = x - (double) idx;
    base = (mat_idx * tables->nprojectiles + proj_idx) * tables->nenergy;
    return (double) tables->range_csda[base + idx] * (1.0 - frac) + (double) tables->range_csda[base + idx + 1u] * frac;
}

/**
 * @brief Return the density [g/cm³] for the zone described by @p zr.
 *
 * @details
 * For analytic zones (@p zr->has_hu == 0) returns the ASSIGNMAT density
 * tables->rho[zr->material_idx].  For voxel zones (@p zr->has_hu != 0)
 * returns the current voxel density tables->hu_rho_lut[zr->hu + 1000].
 *
 * @param[in] tables  Material runtime tables.
 * @param[in] zr      Zone reference filled by osh_gemca_runtime_get_zone_ref_batch().
 *
 * @returns Density [g/cm³].
 */
OSH_HD static inline double osh_material_runtime_get_rho(struct osh_material_runtime const *tables,
                                                   struct osh_zone_ref const *zr) {
    if (zr->has_hu) {
        int hu;
        assert(tables->hu_rho_lut != NULL);
        if (!tables->hu_rho_lut) {
            return 0.0;
        }
        hu = zr->hu < -1000 ? -1000 : (zr->hu > 1600 ? 1600 : (int) zr->hu);
        return (double) tables->hu_rho_lut[hu + 1000];
    }
    return (double) tables->rho[zr->material_idx];
}

#ifdef __cplusplus
}
#endif

#endif /* OSH_MATERIAL_RUNTIME_H */
