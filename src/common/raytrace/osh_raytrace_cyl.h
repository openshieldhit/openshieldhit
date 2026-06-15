#ifndef OSH_RAYTRACE_CYL_H
#define OSH_RAYTRACE_CYL_H

#include <stddef.h>

#include "common/raytrace/osh_raytrace.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Compute all voxel crossings for a ray segment through a cylindrical (R,Z) grid.
 *
 * @details
 * Implements the Jacobs (1998) alpha-parameterisation adapted for a cylindrical grid.
 * R-shell boundaries are found via quadratic circle intersection; Z-plane boundaries
 * use the same linear DDA as the Cartesian implementation.
 *
 * The grid descriptor reuses @ref osh_raytrace_grid with the following field convention:
 *
 *   origin[0], spacing[0], n[0]  →  r_min, dr, nr
 *   origin[2], spacing[2], n[2]  →  z_min, dz, nz
 *   origin[1], spacing[1], n[1]  →  unused (set 0, 0, 1 by the caller)
 *
 * Flat index: @c ir + @c n[0]*iz  (R varies fastest within each Z slice).
 *
 * @note @c vol_inv in the returned crossings is NOT filled by this function.
 *       The caller (scoring layer) fills it after traversal.
 *
 * @param[in]  grid      Grid descriptor following the CYL field convention above.
 * @param[in]  p         Ray start point [cm], in world (X,Y,Z) coordinates.
 * @param[in]  v         Unit direction vector (|v| == 1).
 * @param[in]  ds        Total step length [cm].
 * @param[out] crossings Caller-allocated array; capacity >= 2*n[0] + n[2].
 * @param[out] n_out     Number of crossings written.
 *
 * @returns 1 if one or more crossings were found, 0 if the ray misses the grid.
 *
 * @see Jacobs F et al., J Comput Inf Technol. 1998;6(1):89-94.
 *      https://hrcak.srce.hr/150245
 */
int osh_raytrace_cyl_traverse(struct osh_raytrace_grid const *grid,
                              double const p[3],
                              double const v[3],
                              double ds,
                              struct osh_voxel_crossing *crossings,
                              size_t *n_out);

#ifdef __cplusplus
}
#endif

#endif /* OSH_RAYTRACE_CYL_H */
