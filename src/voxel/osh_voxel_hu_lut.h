#ifndef _OSH_VOXEL_HU_LUT_H
#define _OSH_VOXEL_HU_LUT_H

#include <stdint.h>

#include "openshieldhit/material.h"
#include "openshieldhit/status.h"
#include "openshieldhit/voxel.h"

/*
 * Internal API for HU-to-material and HU-to-density calibration.
 *
 * Each supported CT calibration scheme provides three operations:
 *  1. Register the scheme's tissue bins as material workspace entries.
 *  2. Build the HU→bin-index LUT used by the geometry runtime (GEMCA).
 *  3. Build the HU→density LUT used by the material runtime (transport).
 *
 * Callers (geometry compile, material compile, app parser) select the
 * appropriate pair of build functions based on the hu_table_type field
 * stored on the cold workspaces.
 */

/* ---- Schneider 2000 (24 bins) -------------------------------------------- */

/**
 * @brief Register Schneider 2000 tissue bins as material workspace entries.
 *
 * @details
 * Appends 24 material definitions (named "schneider_00" … "schneider_23")
 * to @p wm.  Called by the application parser when the HUTABLE card selects
 * the Schneider scheme.  The resulting materials are then compiled into
 * transport tables by the material runtime.
 *
 * @param[in,out] wm  Material workspace to extend.
 *
 * @returns OSH_OK on success, OSH_ENOMEM on allocation failure.
 */
enum osh_status osh_voxel_register_schneider_materials(struct osh_material_workspace *wm);

/**
 * @brief Map a raw HU value to a Schneider 2000 density [g/cm³].
 *
 * @details
 * Implements the piecewise-linear Schneider 2000 equations (Eqs. 20–24).
 * This function is called once per LUT entry by
 * @ref osh_voxel_build_hu_rho_lut_schneider2000 and is also available for
 * direct point queries.
 *
 * @param[in] hu  Hounsfield unit in [-1000, 1600]; clamped if out of range.
 *
 * @returns Density [g/cm³].
 */
float osh_voxel_hu2rho_schneider2000(int16_t hu);

/**
 * @brief Map a raw HU value to its Schneider 2000 bin index.
 *
 * @param[in] hu  Hounsfield unit; clamped if out of range.
 *
 * @returns Bin index in [0, 23].
 */
int osh_voxel_hu2idx_schneider2000(int16_t hu);

/**
 * @brief Fill a [2601] LUT with Schneider 2000 bin indices.
 *
 * @details
 * Fills @p lut so that @c lut[hu+1000] is the material bin index for
 * Hounsfield value @c hu.  Used by the geometry runtime compile step.
 *
 * @param[out] lut  Caller-allocated array of OSH_VOXEL_HU_LUT_SIZE entries.
 */
void osh_voxel_build_hu_bin_lut_schneider2000(uint8_t lut[OSH_VOXEL_HU_LUT_SIZE]);

/**
 * @brief Fill a [2601] LUT with Schneider 2000 densities [g/cm³].
 *
 * @details
 * Fills @p lut so that @c lut[hu+1000] is the density for Hounsfield value
 * @c hu.  Used by the material runtime compile step.
 *
 * @param[out] lut  Caller-allocated array of OSH_VOXEL_HU_LUT_SIZE entries.
 */
void osh_voxel_build_hu_rho_lut_schneider2000(float lut[OSH_VOXEL_HU_LUT_SIZE]);

/* ---- Permatassari 2020 (40 bins) ----------------------------------------- */

/**
 * @brief Register Permatassari 2020 tissue bins as material workspace entries.
 *
 * @details
 * Appends 40 material definitions (named "permatassari2020_00" … "_39") to
 * @p wm.  Called by the application parser when the HUTABLE card selects the
 * Permatassari scheme.
 *
 * @param[in,out] wm  Material workspace to extend.
 *
 * @returns OSH_OK on success, OSH_ENOMEM on allocation failure.
 */
enum osh_status osh_voxel_register_permatassari_materials(struct osh_material_workspace *wm);

/**
 * @brief Map a raw HU value to a Permatassari 2020 density [g/cm³].
 *
 * @details
 * Uses a linear density model: @c rho = density_factor[bin] * (1000 + hu).
 * The @p bin is typically supplied from the HU→bin LUT for efficiency.
 *
 * @param[in] hu   Hounsfield unit; clamped if out of range.
 * @param[in] bin  Calibration bin index; clamped if out of range.
 *
 * @returns Density [g/cm³], minimum 0.0001.
 */
float osh_voxel_hu2rho_permatassari2020(int16_t hu, int bin);

/**
 * @brief Fill a [2601] LUT with Permatassari 2020 bin indices.
 *
 * @param[out] lut  Caller-allocated array of OSH_VOXEL_HU_LUT_SIZE entries.
 */
void osh_voxel_build_hu_bin_lut_permatassari2020(uint8_t lut[OSH_VOXEL_HU_LUT_SIZE]);

/**
 * @brief Fill a [2601] LUT with Permatassari 2020 densities [g/cm³].
 *
 * @param[out] lut  Caller-allocated array of OSH_VOXEL_HU_LUT_SIZE entries.
 */
void osh_voxel_build_hu_rho_lut_permatassari2020(float lut[OSH_VOXEL_HU_LUT_SIZE]);

/* ---- WEPL conversions (not used in transport path) ----------------------- */

/**
 * @brief Convert a Hounsfield unit to a water-equivalent path-length ratio.
 *
 * @details
 * Provided as an optional utility for external callers (e.g. range probing,
 * treatment planning cross-checks).  Not used in the main Monte Carlo
 * transport path; stopping-power scaling via the density LUT is used there
 * instead.
 *
 * @p alg selects the conversion algorithm:
 * - 1: Minohara 1993
 * - 2: Jacob 1996
 * - 3: Geiss 1999 (stub, returns 0)
 *
 * @param[in] hu   Hounsfield unit in [-1000, 4000]; returns 0 outside range.
 * @param[in] alg  Algorithm selector (1–3).
 *
 * @returns WEPL ratio (dimensionless, in units of 1/1000 — multiply by 1000
 *          to get percent), or 0 for unknown @p alg.
 */
float osh_voxel_hu2wepl(int16_t hu, char alg);

#endif /* _OSH_VOXEL_HU_LUT_H */
