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

static void test_locate_uses_grid_storage_order(void) {
    struct osh_raytrace_grid row_grid = {
        .origin = {0.0, 0.0, 0.0},
        .spacing = {1.0, 1.0, 1.0},
        .n = {3u, 3u, 3u},
    };
    struct osh_raytrace_grid morton_grid = {
        .origin = {0.0, 0.0, 0.0},
        .spacing = {1.0, 1.0, 1.0},
        .n = {8u, 8u, 8u},
        .tile_order = OSH_VOXEL_ORDER_MORTON8,
    };
    double p[3] = {1.5, 1.5, 1.5};
    double outside[3] = {3.0, 1.5, 1.5};
    size_t idx = 0u;

    ASSERT_TRUE(osh_raytrace_locate(&row_grid, p, &idx) == 1);
    ASSERT_TRUE(idx == 13u);
    ASSERT_TRUE(osh_raytrace_locate(&row_grid, outside, &idx) == 0);

    ASSERT_TRUE(osh_raytrace_locate(&morton_grid, p, &idx) == 1);
    ASSERT_TRUE(idx == osh_voxel_tile_idx(1u, 1u, 1u, 1u, 1u));
}

static void test_first_crossing_matches_traverse_first_entry(void) {
    struct osh_raytrace_grid grid = {
        .origin = {0.0, 0.0, 0.0},
        .spacing = {1.0, 1.0, 1.0},
        .n = {3u, 3u, 3u},
    };
    double p[3] = {0.25, 1.5, 1.5};
    double v[3] = {1.0, 0.0, 0.0};
    struct osh_voxel_crossing first;
    struct osh_voxel_crossing crossings[9];
    size_t n = 0u;
    int rc_first;
    int rc_all;

    rc_first = osh_raytrace_first_crossing(&grid, p, v, 4.0, &first);
    rc_all = osh_raytrace_traverse(&grid, p, v, 4.0, crossings, &n);

    ASSERT_TRUE(rc_first == 1);
    ASSERT_TRUE(rc_all == 1);
    ASSERT_TRUE(n > 0u);
    ASSERT_TRUE(first.idx == crossings[0].idx);
    ASSERT_TRUE(fabs(first.path_len - crossings[0].path_len) < 1.0e-12);
}

static void test_first_crossing_clips_outside_start(void) {
    struct osh_raytrace_grid grid = {
        .origin = {0.0, 0.0, 0.0},
        .spacing = {1.0, 1.0, 1.0},
        .n = {3u, 3u, 3u},
    };
    double p[3] = {-0.5, 1.5, 1.5};
    double v[3] = {1.0, 0.0, 0.0};
    struct osh_voxel_crossing first;

    ASSERT_TRUE(osh_raytrace_first_crossing(&grid, p, v, 4.0, &first) == 1);
    ASSERT_TRUE(first.idx == 12u);
    ASSERT_TRUE(fabs(first.path_len - 1.0) < 1.0e-12);
}

static void test_first_crossing_boundary_with_negative_direction(void) {
    struct osh_raytrace_grid grid = {
        .origin = {0.0, 0.0, 0.0},
        .spacing = {1.0, 1.0, 1.0},
        .n = {3u, 3u, 3u},
    };
    double p[3] = {1.0, 1.5, 1.5};
    double v[3] = {-1.0, 0.0, 0.0};
    struct osh_voxel_crossing first;

    ASSERT_TRUE(osh_raytrace_first_crossing(&grid, p, v, 2.0, &first) == 1);
    ASSERT_TRUE(first.idx == 12u);
    ASSERT_TRUE(fabs(first.path_len - 1.0) < 1.0e-12);
}

int main(void) {
    test_diagonal_edge_and_corner_hits_do_not_emit_zero_length();
    test_morton_tile_order_emits_morton_indices();
    test_locate_uses_grid_storage_order();
    test_first_crossing_matches_traverse_first_entry();
    test_first_crossing_clips_outside_start();
    test_first_crossing_boundary_with_negative_direction();
    return 0;
}
