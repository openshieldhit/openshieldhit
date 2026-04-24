#ifndef OSH_RAYTRACE_H
#define OSH_RAYTRACE_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Default flat-index layout for grids that do not explicitly request tiling. */
#define OSH_RAYTRACE_GRID_TILE_ORDER_DEFAULT ((uint8_t) 0u)

/**
 * @brief Uniform voxel grid descriptor for raytrace traversal.
 *
 * @details
 * Describes a regular Cartesian grid by its lower corner, voxel spacing, and
 * voxel count along each axis.  The grid spans:
 *   [origin[i],  origin[i] + n[i]*spacing[i]]  for i = 0,1,2.
 *
 * Voxel index convention is controlled by tile_order (OSH_VOXEL_ORDER_*):
 *   0 (ROW_MAJOR): idx = ix + n[0]*(iy + n[1]*iz)
 *   8 (MORTON8):   8×8×8 Morton-tiled — see osh_voxel_order.h
 * The DCM/VOX loader's compiled default is OSH_VOXEL_LAYOUT_DEFAULT; callers
 * may still set grid->tile_order explicitly per grid if needed. Other grid
 * producers should set tile_order deterministically; use
 * OSH_RAYTRACE_GRID_TILE_ORDER_DEFAULT for row-major grids.
 *
 * Both the CT transport grid and scoring mesh grids use this descriptor.
 * Callers construct it from their own geometry representation and pass it
 * to osh_raytrace_traverse(); the algorithm sees no domain types.
 */
struct osh_raytrace_grid {
    double origin[3];   /* position of the corner of voxel [0,0,0] [cm] */
    double spacing[3];  /* voxel size along each axis [cm] */
    size_t n[3];        /* number of voxels along each axis */
    uint8_t tile_order; /* OSH_VOXEL_ORDER_*; default row-major via OSH_RAYTRACE_GRID_TILE_ORDER_DEFAULT */
};

/**
 * @brief One voxel entered by a ray, with the path length spent inside it.
 */
struct osh_voxel_crossing {
    size_t idx;      /* flat voxel index in the layout given by grid->tile_order */
    double path_len; /* track length through this voxel [cm] */
};

/**
 * @brief Locate the voxel containing point @p p.
 *
 * @details
 * Uses the same flat-index convention as @ref osh_raytrace_traverse.  The grid
 * lower boundary is inclusive and the upper boundary is exclusive:
 * [origin[i], origin[i] + n[i]*spacing[i]).
 *
 * @param[in]  grid     Grid descriptor.
 * @param[in]  p        Point [cm].
 * @param[out] idx_out  Flat voxel index in @p grid storage order.
 *
 * @returns 1 if @p p is inside the grid, 0 otherwise.
 */
int osh_raytrace_locate(struct osh_raytrace_grid const *grid, double const p[3], size_t *idx_out);

/**
 * @brief Compute only the first voxel crossing for a ray segment through a grid.
 *
 * @details
 * Traverses the ray starting at @p p, travelling along unit direction @p v, for
 * at most @p ds.  The returned crossing is identical to the first entry that
 * @ref osh_raytrace_traverse would emit, but this helper avoids walking the
 * rest of the grid and does not require a caller-owned crossing array.
 *
 * @param[in]  grid  Grid descriptor.
 * @param[in]  p     Ray start point [cm].
 * @param[in]  v     Unit direction vector (|v| == 1).
 * @param[in]  ds    Total step length [cm].
 * @param[out] out   First voxel crossing.
 *
 * @returns 1 if a crossing was found, 0 if the ray segment misses the grid.
 */
int osh_raytrace_first_crossing(struct osh_raytrace_grid const *grid,
                                double const p[3],
                                double const v[3],
                                double ds,
                                struct osh_voxel_crossing *out);

/**
 * @brief Compute all voxel crossings for a ray segment through a grid.
 *
 * @details
 * Traverses the ray starting at p, travelling along unit direction v, for a
 * total distance ds.  For each voxel entered, one osh_voxel_crossing entry is
 * written recording the voxel index and the path length inside it.  The sum
 * of all path_len values equals the clipped length of the ray inside the grid
 * (which may be less than ds if the ray enters or exits the grid mid-step).
 *
 * The caller must pre-allocate crossings[] with capacity at least
 * n[0] + n[1] + n[2], which is a safe upper bound on the number of voxels
 * a straight ray can cross in one pass through the grid.
 *
 * Which algorithm is compiled in is selected at build time via the
 * OSH_RAYTRACE_ALGORITHM CMake cache variable (default: SIDDON).
 *
 * @param[in]  grid      Grid descriptor.
 * @param[in]  p         Ray start point [cm].
 * @param[in]  v         Unit direction vector (|v| == 1).
 * @param[in]  ds        Total step length [cm].
 * @param[out] crossings Caller-allocated array; capacity >= n[0]+n[1]+n[2].
 * @param[out] n_out     Number of crossings written.
 *
 * @returns 1 if one or more crossings were found, 0 if the ray misses the
 *          grid entirely.
 */
int osh_raytrace_traverse(struct osh_raytrace_grid const *grid,
                          double const p[3],
                          double const v[3],
                          double ds,
                          struct osh_voxel_crossing *crossings,
                          size_t *n_out);

#ifdef __cplusplus
}
#endif

#endif /* OSH_RAYTRACE_H */
