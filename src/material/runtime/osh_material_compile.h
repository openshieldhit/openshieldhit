#ifndef OSH_MATERIAL_COMPILE_H
#define OSH_MATERIAL_COMPILE_H

#include "material/osh_material.h"
#include "material/runtime/osh_material_runtime.h"
#include "openshieldhit/logger.h"
#include "openshieldhit/status.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Build the runtime stopping-power and range tables.
 *
 * @details
 * For each user material (indices >= OSH_MATERIAL_INDEX_FIRST_USER):
 *  - If the material has dE/dx overrides: those source curves are resampled
 *    onto the log-uniform runtime grid using log-log interpolation.
 *  - Otherwise: the analytic Bethe model is evaluated at each grid point.
 *
 * The projectile list is the union of:
 *  - the caller-requested Z = 1..z_max range
 *  - all projectile Z values present in any material dE/dx override
 *
 * Representative projectile mass numbers come from the particle default
 * isotope map for all projectiles. Materials without an override for one
 * runtime projectile column are completed with Bethe fallback for that
 * projectile.
 *
 * All materials use the same projectile list and energy grid, so material
 * indices remain directly usable as array offsets.
 *
 * @param[in]  wm       Completed material workspace.
 * @param[in]  z_max    Maximum projectile atomic number to include.
 * @param[in]  diag     Borrowed diagnostics sink for compile messages, or NULL.
 * @param[out] tables   Receives the allocated runtime tables on success.
 *
 * @returns OSH_OK on success, or an error code on failure.
 */
enum osh_status osh_material_compile(struct osh_material_workspace const *wm,
                                     unsigned int z_max,
                                     struct osh_diag_sink const *diag,
                                     struct osh_material_runtime *tables);

/**
 * @brief Release a runtime tables struct and all owned memory.
 *
 * @param[in] tables  Tables to release.  Safe to call on a zero-initialised struct.
 */
void osh_material_runtime_free(struct osh_material_runtime *tables);

#ifdef __cplusplus
}
#endif

#endif /* OSH_MATERIAL_COMPILE_H */
