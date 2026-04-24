#ifndef OPENSHIELDHIT_VOXEL_H
#define OPENSHIELDHIT_VOXEL_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @file voxel.h
 * @brief Public constants shared by the voxel geometry and material subsystems.
 *
 * @details
 * Voxel (CT) support spans two modules: the geometry engine (GEMCA) and the
 * material runtime.  Both need the HU calibration table selector and the LUT
 * size without depending on each other, so these constants are placed in a
 * shared public header.
 *
 * The HU value range covered by the LUT is [-1000, 1600], giving 2601 entries
 * when indexed by (HU + 1000).
 */

/**
 * @brief Number of entries in every HU lookup table.
 *
 * @details
 * Covers the range [-1000, 1600] inclusive.  Indexed as @c lut[hu + 1000].
 */
#define OSH_VOXEL_HU_LUT_SIZE 2601

/**
 * @name HU calibration table type selectors
 *
 * @brief Identify which CT-to-material calibration scheme is active.
 *
 * @details
 * Stored on both @ref osh_geometry_workspace and @ref osh_material_workspace
 * so each compile step can build its own LUT independently.
 *
 * @p OSH_HU_TABLE_NONE is the default (0) so that zero-initialised workspaces
 * correctly indicate a non-CT run.
 *
 * @{
 */

/** No HU table registered; geometry contains no voxel bodies. */
#define OSH_HU_TABLE_NONE 0

/**
 * Schneider W et al., Phys. Med. Biol. 45:459–478 (2000).
 * 24 tissue bins.
 */
#define OSH_HU_TABLE_SCHNEIDER 1

/**
 * Permatassari et al., Phys. Med. Biol. 65:ab9702 (2020), Method C.
 * 40 tissue bins.
 */
#define OSH_HU_TABLE_PERMATASSARI 2

/** @} */

#ifdef __cplusplus
}
#endif

#endif /* OPENSHIELDHIT_VOXEL_H */
