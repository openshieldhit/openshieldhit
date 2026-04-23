#include "common/osh_voxel_order.h"

#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int mul_size_overflow(size_t a, size_t b, size_t *out) {
    if (!out) {
        return 1;
    }
    if (a != 0u && b > SIZE_MAX / a) {
        return 1;
    }
    *out = a * b;
    return 0;
}

void *osh_voxel_reorder(
    void const *src, size_t nx, size_t ny, size_t nz, size_t element_size, uint8_t tile_order, size_t *out_len) {
    size_t total;
    size_t nbytes;
    size_t tmp;
    void *dst;
    size_t Tx, Ty, Tz;
    size_t ix, iy, iz;

    if (!out_len) {
        return NULL;
    }
    *out_len = 0u;
    if (!src || element_size == 0u) {
        return NULL;
    }

    if (tile_order == OSH_VOXEL_ORDER_ROW_MAJOR) {
        if (mul_size_overflow(nx, ny, &tmp) || mul_size_overflow(tmp, nz, &total)
            || mul_size_overflow(total, element_size, &nbytes)) {
            return NULL;
        }
        dst = malloc(nbytes);
        if (!dst) {
            return NULL;
        }
        memcpy(dst, src, nbytes);
        *out_len = total;
        return dst;
    }

    if (tile_order != OSH_VOXEL_ORDER_MORTON8) {
        return NULL;
    }

    if (nx > SIZE_MAX - 7u || ny > SIZE_MAX - 7u || nz > SIZE_MAX - 7u) {
        return NULL;
    }
    Tx = (nx + 7u) >> 3u;
    Ty = (ny + 7u) >> 3u;
    Tz = (nz + 7u) >> 3u;
    if (mul_size_overflow(Tx, Ty, &tmp) || mul_size_overflow(tmp, Tz, &tmp) || mul_size_overflow(tmp, 512u, &total)
        || mul_size_overflow(total, element_size, &nbytes)) {
        return NULL;
    }

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
