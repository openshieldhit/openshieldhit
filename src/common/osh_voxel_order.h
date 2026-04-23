#ifndef OSH_VOXEL_ORDER_H
#define OSH_VOXEL_ORDER_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/** tile_order value: standard row-major layout, idx = ix + nx*(iy + ny*iz). */
#define OSH_VOXEL_ORDER_ROW_MAJOR 0u

/**
 * tile_order value: 8×8×8 Morton-tiled layout.
 *
 * The array is divided into Tx×Ty×Tz tiles of 512 voxels each, where
 * Tx = ceil(nx/8) etc.  Tiles are stored in row-major order; voxels within
 * each tile follow a 3-way Morton (Z-order) curve so that spatially nearby
 * voxels are close in memory regardless of which axis a ray travels along.
 *
 * The allocated array has Tx*Ty*Tz*512 elements; slots for padding voxels
 * (beyond the actual nx/ny/nz boundary) are zero-initialised and never
 * accessed by a traversal that respects the grid bounds.
 */
#define OSH_VOXEL_ORDER_MORTON8 8u

/*
 * Build-time default voxel layout used by DCM/VOX setup when the caller does
 * not override tile_order explicitly. Configure via CMake cache variable:
 *   -DOSH_VOXEL_LAYOUT=ROW_MAJOR   (baseline, default)
 *   -DOSH_VOXEL_LAYOUT=MORTON8
 *
 * TODO: add axis-permuted row-major layouts (tile_order 1..7) once their
 * index contract and selection heuristic are ready.
 */
#ifndef OSH_VOXEL_LAYOUT_DEFAULT
#define OSH_VOXEL_LAYOUT_DEFAULT OSH_VOXEL_ORDER_ROW_MAJOR
#endif

#if (OSH_VOXEL_LAYOUT_DEFAULT != OSH_VOXEL_ORDER_ROW_MAJOR) && (OSH_VOXEL_LAYOUT_DEFAULT != OSH_VOXEL_ORDER_MORTON8)
#error "OSH_VOXEL_LAYOUT_DEFAULT must be OSH_VOXEL_ORDER_ROW_MAJOR or OSH_VOXEL_ORDER_MORTON8"
#endif

/*
 * 3-way bit-interleave LUTs for 3-bit intra-tile coordinates (0..7).
 *   x bits → positions 0,3,6 of the 9-bit Morton index
 *   y bits → positions 1,4,7
 *   z bits → positions 2,5,8
 *
 * Invariant: osh_voxel_m3x[i] | osh_voxel_m3y[j] | osh_voxel_m3z[k] for all i,j,k in
 * [0,7] yields all 512 distinct values 0..511 exactly once.
 */
static uint16_t const osh_voxel_m3x[8] = {0, 1, 8, 9, 64, 65, 72, 73};
static uint16_t const osh_voxel_m3y[8] = {0, 2, 16, 18, 128, 130, 144, 146};
static uint16_t const osh_voxel_m3z[8] = {0, 4, 32, 36, 256, 260, 288, 292};

/**
 * @brief Compute the Morton-tiled flat index for voxel (ix, iy, iz).
 *
 * @param ix,iy,iz  Voxel coordinates (must be within the grid bounds).
 * @param Tx,Ty     Tile counts along x and y: ceil(nx/8), ceil(ny/8).
 *                  Pre-compute these once outside the traversal loop.
 */
static inline size_t osh_voxel_tile_idx(size_t ix, size_t iy, size_t iz, size_t Tx, size_t Ty) {
    size_t tile_id = (ix >> 3u) + Tx * ((iy >> 3u) + Ty * (iz >> 3u));
    size_t intra = (size_t) osh_voxel_m3x[ix & 7u] | (size_t) osh_voxel_m3y[iy & 7u] | (size_t) osh_voxel_m3z[iz & 7u];
    return tile_id * 512u + intra;
}

/**
 * @brief Reorder a voxel array from row-major to the requested layout.
 *
 * Allocates a new buffer containing all voxels of an nx×ny×nz volume in the
 * layout selected by tile_order.  The returned element count (*out_len) may
 * exceed nx*ny*nz for tiled layouts due to boundary padding; padding slots
 * are zero-initialised.
 *
 * Returns NULL on allocation failure or unknown tile_order.  The caller owns
 * the returned buffer and must free() it.
 *
 * @param src           Source array in row-major order (ix fastest).
 * @param nx,ny,nz      Volume dimensions.
 * @param element_size  Size of each element in bytes.
 * @param tile_order    OSH_VOXEL_ORDER_* constant.
 * @param out_len       Receives the number of elements in the returned buffer.
 */
void *osh_voxel_reorder(
    void const *src, size_t nx, size_t ny, size_t nz, size_t element_size, uint8_t tile_order, size_t *out_len);

#ifdef __cplusplus
}
#endif

#endif /* OSH_VOXEL_ORDER_H */
