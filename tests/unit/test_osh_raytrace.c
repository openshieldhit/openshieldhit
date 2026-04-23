#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "common/osh_voxel_order.h"
#include "common/raytrace/osh_raytrace.h"

#define ASSERT_TRUE(cond)                                                                                              \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            fprintf(stderr, "ASSERT FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__);                                 \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

static void test_diagonal_edge_and_corner_hits_do_not_emit_zero_length(void) {
    struct osh_raytrace_grid grid = {
        .origin = {0.0, 0.0, 0.0},
        .spacing = {1.0, 1.0, 1.0},
        .n = {2u, 2u, 2u},
    };
    double inv_sqrt3 = 1.0 / sqrt(3.0);
    double p[3] = {0.25, 0.25, 0.25};
    double v[3] = {inv_sqrt3, inv_sqrt3, inv_sqrt3};
    struct osh_voxel_crossing crossings[6];
    size_t n = 0u;
    int rc;
    size_t i;
    double sum = 0.0;

    rc = osh_raytrace_traverse(&grid, p, v, 4.0, crossings, &n);
    ASSERT_TRUE(rc == 1);
    ASSERT_TRUE(n == 2u);
    ASSERT_TRUE(crossings[0].idx == 0u);
    ASSERT_TRUE(crossings[1].idx == 7u);
    for (i = 0; i < n; ++i) {
        ASSERT_TRUE(crossings[i].path_len > 0.0);
        sum += crossings[i].path_len;
    }
    ASSERT_TRUE(fabs(sum - (1.75 * sqrt(3.0))) < 1.0e-10);
}

static void test_morton_tile_order_emits_morton_indices(void) {
    struct osh_raytrace_grid grid = {
        .origin = {0.0, 0.0, 0.0},
        .spacing = {1.0, 1.0, 1.0},
        .n = {8u, 8u, 8u},
        .tile_order = OSH_VOXEL_ORDER_MORTON8,
    };
    double p[3] = {1.25, 2.25, 3.25};
    double v[3] = {1.0, 0.0, 0.0};
    struct osh_voxel_crossing crossings[24];
    size_t n = 0u;
    int rc;

    rc = osh_raytrace_traverse(&grid, p, v, 3.0, crossings, &n);
    ASSERT_TRUE(rc == 1);
    ASSERT_TRUE(n == 4u);
    ASSERT_TRUE(crossings[0].idx == osh_voxel_tile_idx(1u, 2u, 3u, 1u, 1u));
    ASSERT_TRUE(crossings[1].idx == osh_voxel_tile_idx(2u, 2u, 3u, 1u, 1u));
    ASSERT_TRUE(crossings[2].idx == osh_voxel_tile_idx(3u, 2u, 3u, 1u, 1u));
    ASSERT_TRUE(crossings[3].idx == osh_voxel_tile_idx(4u, 2u, 3u, 1u, 1u));
}

int main(void) {
    test_diagonal_edge_and_corner_hits_do_not_emit_zero_length();
    test_morton_tile_order_emits_morton_indices();
    return 0;
}
