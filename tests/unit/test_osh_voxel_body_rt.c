#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_coord.h"
#include "common/osh_ray.h"
#include "gemca/osh_gemca2.h"
#include "gemca/osh_gemca2_defines.h"
#include "gemca/runtime/osh_gemca_runtime.h"
#include "gemca/runtime/osh_gemca_runtime_voxel.h"
#include "gemca/voxel/osh_gemca2_voxel_hu.h"
#include "test_assert.h"

/*
 * Test fixture: a 3×3×3 voxel grid, 1 cm voxels, origin at (0,0,0).
 * Row-major flat index: idx = ix + 3*(iy + 3*iz).
 *
 * A ray travelling along +X through the center row (iy=1, iz=1) crosses
 * voxel indices 12, 13, 14 in order.
 */

#define GRID_N 3
#define GRID_N_VOX ((size_t) GRID_N * (size_t) GRID_N * (size_t) GRID_N)

static uint8_t s_bin_lut[2601];
static float s_rho_lut[2601];

static void build_runtime(struct osh_gemca_runtime *rt, struct gemca_rt_body *body, int16_t *hu, size_t n_vox) {
    static double const identity[16] = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
    size_t i;

    memset(rt, 0, sizeof(*rt));
    memset(body, 0, sizeof(*body));

    memcpy(body->t, identity, sizeof(identity));
    body->ct_grid.origin[0] = 0.0;
    body->ct_grid.origin[1] = 0.0;
    body->ct_grid.origin[2] = 0.0;
    body->ct_grid.spacing[0] = 1.0;
    body->ct_grid.spacing[1] = 1.0;
    body->ct_grid.spacing[2] = 1.0;
    body->ct_grid.n[0] = GRID_N;
    body->ct_grid.n[1] = GRID_N;
    body->ct_grid.n[2] = GRID_N;
    body->ct_grid.tile_order = 0; /* row-major */
    body->type = OSH_GEMCA_BODY_VOX;

    for (i = 0; i < n_vox; i++) {
        hu[i] = 0; /* default: water-ish, bin 5 */
    }
    body->hu = hu;

    osh_gemca_voxel_build_hu_lut(s_bin_lut);
    osh_gemca_voxel_build_rho_lut_schneider2000(s_rho_lut);
    rt->hu_bin_lut = s_bin_lut;
    rt->hu_rho_lut = s_rho_lut;
    rt->bodies = body;
    rt->nbodies = 1;
}

/* Ray along +X through y=1.5, z=1.5 (center of iy=1, iz=1 row). */
static struct ray make_x_ray(void) {
    struct ray r;
    r.p[0] = -0.5; /* start just outside left face */
    r.p[1] = 1.5;
    r.p[2] = 1.5;
    r.cp[0] = 1.0;
    r.cp[1] = 0.0;
    r.cp[2] = 0.0;
    r.system = OSH_COORD_UNIVERSE;
    return r;
}

/* ---- Tests ---------------------------------------------------------------- */

static void test_ray_misses_returns_infinity(void) {
    struct osh_gemca_runtime rt;
    struct gemca_rt_body body;
    int16_t hu[GRID_N * GRID_N * GRID_N];
    struct gemca_rt_voxel_segment segs[16];
    struct ray r;
    size_t n_segs;
    int bin;
    double ds;

    build_runtime(&rt, &body, hu, GRID_N_VOX);

    /* Ray passes entirely beside the grid (y outside [0,3]). */
    r = make_x_ray();
    r.p[1] = 10.0; /* far outside */

    ds = dist_voxel_body_rt(&rt, 0, &r, segs, 16, &n_segs, &bin);
    if (isinf(OSH_GEMCA_INFINITY)) {
        ASSERT_TRUE(isinf(ds));
        ASSERT_TRUE(ds > 0.0);
    } else {
        ASSERT_TRUE(ds == OSH_GEMCA_INFINITY);
    }
    ASSERT_TRUE(n_segs == 0);
    ASSERT_TRUE(bin == -1);
}

static void test_uniform_grid_all_same_bin(void) {
    struct osh_gemca_runtime rt;
    struct gemca_rt_body body;
    int16_t hu[GRID_N * GRID_N * GRID_N];
    struct gemca_rt_voxel_segment segs[16];
    struct ray r;
    size_t n_segs;
    size_t i;
    int bin;
    double ds;
    double expected_rho;

    build_runtime(&rt, &body, hu, GRID_N_VOX);
    /* All HU=0 → bin 5, rho ≈ 1.018 g/cm³. */

    r = make_x_ray();
    ds = dist_voxel_body_rt(&rt, 0, &r, segs, 16, &n_segs, &bin);

    ASSERT_TRUE(n_segs == 3);
    ASSERT_TRUE(bin == 5);

    /* Total distance = 3 × 1.0 cm. */
    ASSERT_TRUE(fabs(ds - 3.0) < 1e-9);

    expected_rho = (double) s_rho_lut[0 + 1000];
    for (i = 0; i < n_segs; i++) {
        ASSERT_TRUE(fabs(segs[i].ds - 1.0) < 1e-9);
        ASSERT_TRUE(fabs(segs[i].rho - expected_rho) < 1e-6);
    }
}

static void test_bin_change_stops_traversal(void) {
    struct osh_gemca_runtime rt;
    struct gemca_rt_body body;
    int16_t hu[GRID_N * GRID_N * GRID_N];
    struct gemca_rt_voxel_segment segs[16];
    struct ray r;
    size_t n_segs;
    int bin;
    double ds;

    build_runtime(&rt, &body, hu, GRID_N_VOX);

    /* Make the second X-voxel (ix=1, iy=1, iz=1, flat idx=13) a different bin.
     * HU=1600 → bin 23 (dense bone). */
    hu[13] = 1600;

    r = make_x_ray();
    ds = dist_voxel_body_rt(&rt, 0, &r, segs, 16, &n_segs, &bin);

    /* Only crosses first voxel (ix=0); stops at bin change. */
    ASSERT_TRUE(n_segs == 1);
    ASSERT_TRUE(bin == 5);
    ASSERT_TRUE(fabs(ds - 1.0) < 1e-9);
}

static void test_segs_cap_limits_output(void) {
    struct osh_gemca_runtime rt;
    struct gemca_rt_body body;
    int16_t hu[GRID_N * GRID_N * GRID_N];
    struct gemca_rt_voxel_segment segs[1]; /* cap = 1 */
    struct ray r;
    size_t n_segs;
    int bin;
    double ds;

    build_runtime(&rt, &body, hu, GRID_N_VOX);
    /* All same bin → would produce 3 segments, but cap = 1. */

    r = make_x_ray();
    ds = dist_voxel_body_rt(&rt, 0, &r, segs, 1, &n_segs, &bin);

    ASSERT_TRUE(n_segs == 1);
    ASSERT_TRUE(bin == 5);
    ASSERT_TRUE(fabs(ds - 1.0) < 1e-9); /* only the first voxel's distance */
    ASSERT_TRUE(fabs(segs[0].ds - 1.0) < 1e-9);
}

static void test_null_segs_distance_only(void) {
    struct osh_gemca_runtime rt;
    struct gemca_rt_body body;
    int16_t hu[GRID_N * GRID_N * GRID_N];
    struct ray r;
    size_t n_segs;
    int bin;
    double ds;

    build_runtime(&rt, &body, hu, GRID_N_VOX);

    r = make_x_ray();
    ds = dist_voxel_body_rt(&rt, 0, &r, NULL, 0, &n_segs, &bin);

    ASSERT_TRUE(bin == 5);
    ASSERT_TRUE(fabs(ds - 3.0) < 1e-9);
}

static void test_oversized_grid_returns_infinity(void) {
    struct osh_gemca_runtime rt;
    struct gemca_rt_body body;
    int16_t hu[GRID_N * GRID_N * GRID_N];
    struct gemca_rt_voxel_segment segs[16];
    struct ray r;
    size_t n_segs;
    int bin;
    double ds;

    build_runtime(&rt, &body, hu, GRID_N_VOX);
    body.ct_grid.n[0] = 1024u;
    body.ct_grid.n[1] = 1024u;
    body.ct_grid.n[2] = 1024u; /* 3072 crossings required > _CROSSINGS_CAP */

    r = make_x_ray();
    ds = dist_voxel_body_rt(&rt, 0, &r, segs, 16, &n_segs, &bin);

    if (isinf(OSH_GEMCA_INFINITY)) {
        ASSERT_TRUE(isinf(ds));
        ASSERT_TRUE(ds > 0.0);
    } else {
        ASSERT_TRUE(ds == OSH_GEMCA_INFINITY);
    }
    ASSERT_TRUE(n_segs == 0);
    ASSERT_TRUE(bin == -1);
}

static void test_rho_matches_lut(void) {
    struct osh_gemca_runtime rt;
    struct gemca_rt_body body;
    int16_t hu[GRID_N * GRID_N * GRID_N];
    struct gemca_rt_voxel_segment segs[16];
    struct ray r;
    size_t n_segs;
    size_t i;
    int bin;

    build_runtime(&rt, &body, hu, GRID_N_VOX);

    /* Use a non-zero HU so the density is clearly non-trivial. */
    for (i = 0; i < GRID_N_VOX; i++) {
        hu[i] = 200;
    }

    r = make_x_ray();
    dist_voxel_body_rt(&rt, 0, &r, segs, 16, &n_segs, &bin);

    ASSERT_TRUE(n_segs == 3);
    for (i = 0; i < n_segs; i++) {
        double expected = (double) s_rho_lut[200 + 1000];
        ASSERT_TRUE(fabs(segs[i].rho - expected) < 1e-6);
    }
}

int main(void) {
    test_ray_misses_returns_infinity();
    test_uniform_grid_all_same_bin();
    test_bin_change_stops_traversal();
    test_segs_cap_limits_output();
    test_null_segs_distance_only();
    test_oversized_grid_returns_infinity();
    test_rho_matches_lut();
    return 0;
}
