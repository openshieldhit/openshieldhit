#include "common/osh_voxel_order.h"

#include <stdlib.h>
#include <string.h>

void *osh_voxel_reorder(
    void const *src, size_t nx, size_t ny, size_t nz, size_t element_size, uint8_t tile_order, size_t *out_len) {
    size_t total;
    void *dst;
    size_t Tx, Ty, Tz;
    size_t ix, iy, iz;

    if (!src || !out_len || element_size == 0u) {
        return NULL;
    }

    if (tile_order == OSH_VOXEL_ORDER_ROW_MAJOR) {
        total = nx * ny * nz;
        dst = malloc(total * element_size);
        if (!dst) {
            return NULL;
        }
        memcpy(dst, src, total * element_size);
        *out_len = total;
        return dst;
    }

    if (tile_order != OSH_VOXEL_ORDER_MORTON8) {
        return NULL;
    }

    Tx = (nx + 7u) >> 3u;
    Ty = (ny + 7u) >> 3u;
    Tz = (nz + 7u) >> 3u;
    total = Tx * Ty * Tz * 512u;

    /* calloc so padding slots beyond the volume boundary are zero-initialised. */
    dst = calloc(total, element_size);
    if (!dst) {
        return NULL;
    }

    for (iz = 0u; iz < nz; iz++) {
        for (iy = 0u; iy < ny; iy++) {
            for (ix = 0u; ix < nx; ix++) {
                size_t si = ix + nx * (iy + ny * iz);
                size_t di = osh_voxel_tile_idx(ix, iy, iz, Tx, Ty);
                memcpy((char *) dst + di * element_size, (char const *) src + si * element_size, element_size);
            }
        }
    }

    *out_len = total;
    return dst;
}
