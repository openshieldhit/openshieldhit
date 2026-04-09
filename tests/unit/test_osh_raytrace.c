#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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

int main(void) {
    test_diagonal_edge_and_corner_hits_do_not_emit_zero_length();
    return 0;
}
