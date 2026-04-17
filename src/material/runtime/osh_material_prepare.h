#ifndef OSH_MATERIAL_PREPARE_H
#define OSH_MATERIAL_PREPARE_H

#include "material/osh_material.h"
#include "material/runtime/osh_material_runtime.h"
#include "openshieldhit/status.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Build the runtime stopping-power and range tables.
 *
 * @details
 * For each user material (indices >= OSH_MATERIAL_INDEX_FIRST_USER):
 *  - If the material has a LOADDEDX path: the source table is loaded, and
 *    stopping powers are resampled onto the log-uniform runtime grid using
 *    log-log interpolation.
 *  - Otherwise: the analytic Bethe model is evaluated at each grid point.
 *
 * The projectile list is the union of:
 *  - the caller-requested Z = 1..z_max range
 *  - all contiguous projectile columns present in any LOADDEDX source table
 *
 * Representative projectile mass numbers come from the particle default
 * isotope map so LOADDEDX import and Bethe fallback use the same projectile
 * convention. Materials whose external table covers fewer projectile columns
 * than the global runtime set are completed with Bethe fallback for the
 * remaining higher-Z projectiles.
 *
 * All materials use the same projectile list and energy grid, so material
 * indices remain directly usable as array offsets.
 *
 * @param[in]  wm       Completed material workspace.
 * @param[in]  z_max    Maximum projectile atomic number to include.
 * @param[out] tables   Receives the allocated runtime tables on success.
 *
 * @returns OSH_OK on success, or an error code on failure.
 */
enum osh_status
osh_material_prepare(struct osh_material_workspace const *wm, unsigned int z_max, struct osh_material_runtime *tables);

/**
 * @brief Release a runtime tables struct and all owned memory.
 *
 * @param[in] tables  Tables to release.  Safe to call on a zero-initialised struct.
 */
void osh_material_runtime_free(struct osh_material_runtime *tables);

#ifdef __cplusplus
}
#endif

#endif /* OSH_MATERIAL_PREPARE_H */
