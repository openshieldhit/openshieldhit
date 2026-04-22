#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_voxel_order.h"
#include "test_assert.h"

static void test_morton_lut_completeness(void);
static void test_reorder_roundtrip(void);
static void test_reorder_nonpower8_boundary(void);
static void test_row_major_baseline(void);

int main(void) {
    test_morton_lut_completeness();
    test_reorder_roundtrip();
    test_reorder_nonpower8_boundary();
    test_row_major_baseline();
    return 0;
}

static void test_morton_lut_completeness(void) {
    /* Verify all 512 intra-tile indices are produced exactly once. */
    uint8_t seen[512];
    size_t i, j, k, idx;

    memset(seen, 0, sizeof(seen));
    for (i = 0u; i < 8u; i++) {
        for (j = 0u; j < 8u; j++) {
            for (k = 0u; k < 8u; k++) {
                idx = (size_t) _osh_m3x[i] | (size_t) _osh_m3y[j] | (size_t) _osh_m3z[k];
                ASSERT_TRUE(idx < 512u);
                ASSERT_TRUE(seen[idx] == 0u);
                seen[idx] = 1u;
            }
        }
    }
    for (i = 0u; i < 512u; i++) {
        ASSERT_TRUE(seen[i] == 1u);
    }
}

static void test_reorder_roundtrip(void) {
    /* 16×16×16 volume: reorder to Morton-8, verify every voxel. */
    size_t const nx = 16u, ny = 16u, nz = 16u;
    size_t const n = nx * ny * nz;
    int16_t *src;
    int16_t *dst;
    size_t out_len;
    size_t Tx, Ty;
    size_t ix, iy, iz;
    size_t i;

    src = (int16_t *) malloc(n * sizeof(int16_t));
    ASSERT_TRUE(src != NULL);
    for (i = 0u; i < n; i++) {
        src[i] = (int16_t) (i % 2000u - 1000);
    }

    dst = (int16_t *) osh_voxel_reorder(src, nx, ny, nz, sizeof(int16_t), OSH_VOXEL_ORDER_MORTON8, &out_len);
    ASSERT_TRUE(dst != NULL);

    Tx = (nx + 7u) >> 3u; /* 2 */
    Ty = (ny + 7u) >> 3u; /* 2 */
    ASSERT_TRUE(out_len == Tx * ((ny + 7u) >> 3u) * ((nz + 7u) >> 3u) * 512u);

    for (iz = 0u; iz < nz; iz++) {
        for (iy = 0u; iy < ny; iy++) {
            for (ix = 0u; ix < nx; ix++) {
                size_t si = ix + nx * (iy + ny * iz);
                size_t di = osh_voxel_tile_idx(ix, iy, iz, Tx, Ty);
                ASSERT_TRUE(dst[di] == src[si]);
            }
        }
    }

    free(src);
    free(dst);
}

static void test_reorder_nonpower8_boundary(void) {
    /* Non-multiples of 8: 10×11×9. Checks boundary tiles are handled correctly. */
    size_t const nx = 10u, ny = 11u, nz = 9u;
    size_t const n = nx * ny * nz;
    int16_t *src;
    int16_t *dst;
    size_t out_len;
    size_t Tx, Ty;
    size_t ix, iy, iz;
    size_t i;

    src = (int16_t *) malloc(n * sizeof(int16_t));
    ASSERT_TRUE(src != NULL);
    for (i = 0u; i < n; i++) {
        src[i] = (int16_t) i;
    }

    dst = (int16_t *) osh_voxel_reorder(src, nx, ny, nz, sizeof(int16_t), OSH_VOXEL_ORDER_MORTON8, &out_len);
    ASSERT_TRUE(dst != NULL);

    Tx = (nx + 7u) >> 3u; /* 2 */
    Ty = (ny + 7u) >> 3u; /* 2 */
    ASSERT_TRUE(out_len == Tx * ((ny + 7u) >> 3u) * ((nz + 7u) >> 3u) * 512u);

    for (iz = 0u; iz < nz; iz++) {
        for (iy = 0u; iy < ny; iy++) {
            for (ix = 0u; ix < nx; ix++) {
                size_t si = ix + nx * (iy + ny * iz);
                size_t di = osh_voxel_tile_idx(ix, iy, iz, Tx, Ty);
                ASSERT_TRUE(dst[di] == src[si]);
            }
        }
    }

    free(src);
    free(dst);
}

static void test_row_major_baseline(void) {
    /* tile_order=ROW_MAJOR: output is a plain copy, same element order. */
    size_t const nx = 4u, ny = 3u, nz = 2u;
    size_t const n = nx * ny * nz;
    int16_t *src;
    int16_t *dst;
    size_t out_len;
    size_t i;

    src = (int16_t *) malloc(n * sizeof(int16_t));
    ASSERT_TRUE(src != NULL);
    for (i = 0u; i < n; i++) {
        src[i] = (int16_t) i;
    }

    dst = (int16_t *) osh_voxel_reorder(src, nx, ny, nz, sizeof(int16_t), OSH_VOXEL_ORDER_ROW_MAJOR, &out_len);
    ASSERT_TRUE(dst != NULL);
    ASSERT_TRUE(out_len == n);
    for (i = 0u; i < n; i++) {
        ASSERT_TRUE(dst[i] == src[i]);
    }

    free(src);
    free(dst);
}
