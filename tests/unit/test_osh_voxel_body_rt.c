#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_coord.h"
#include "common/osh_ray.h"
#include "common/osh_step_segment.h"
#include "common/osh_voxel_order.h"
#include "gemca/osh_gemca2.h"
#include "gemca/osh_gemca2_defines.h"
#include "gemca/runtime/osh_gemca_runtime.h"
#include "gemca/runtime/osh_gemca_runtime_voxel.h"
#include "test_assert.h"
#include "voxel/osh_voxel_hu_lut.h"

/*
 * Test fixture: a 3×3×3 voxel grid, 1 cm voxels, origin at (0,0,0).
 * Row-major flat index: idx = ix + 3*(iy + 3*iz).
 *
 * A ray travelling along +X through the center row (iy=1, iz=1) crosses
 * voxel indices 12, 13, 14 in order.
 */

#define GRID_N 3
#define GRID_N_VOX ((size_t) GRID_N * (size_t) GRID_N * (size_t) GRID_N)
#define MORTON8_TILE_VOX 512u

static uint8_t s_bin_lut[2601];

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

    osh_voxel_build_hu_bin_lut_schneider2000(s_bin_lut);
    rt->hu_bin_lut = s_bin_lut;
    rt->bodies = body;
    rt->nbodies = 1;
}

static void set_plane(struct gemca_rt_surface *surface, int type, double p0, double p1) {
    memset(surface, 0, sizeof(*surface));
    surface->type = type;
    surface->p[0] = p0;
    surface->p[1] = p1;
}

static void attach_box_voxel_zone(struct osh_gemca_runtime *rt,
                                  struct gemca_rt_body *body,
                                  struct gemca_rt_surface surfaces[6],
                                  struct gemca_rt_zone *zone,
                                  struct gemca_rt_insn *insn) {
    set_plane(&surfaces[0], OSH_GEMCA_SURF_PLANEX, -1.0, 0.0);
    set_plane(&surfaces[1], OSH_GEMCA_SURF_PLANEX, 1.0, -3.0);
    set_plane(&surfaces[2], OSH_GEMCA_SURF_PLANEY, -1.0, 0.0);
    set_plane(&surfaces[3], OSH_GEMCA_SURF_PLANEY, 1.0, -3.0);
    set_plane(&surfaces[4], OSH_GEMCA_SURF_PLANEZ, -1.0, 0.0);
    set_plane(&surfaces[5], OSH_GEMCA_SURF_PLANEZ, 1.0, -3.0);

    body->surf_begin = 0u;
    body->nsurfs = 6;
    body->coord = OSH_COORD_UNIVERSE;

    memset(insn, 0, sizeof(*insn));
    insn->op = GEMCA_RT_PUSH_VOXEL_BODY;
    insn->operand = 0;

    memset(zone, 0, sizeof(*zone));
    zone->insns = insn;
    zone->ninsns = 1;
    zone->material_idx = 99u;
    zone->voxel_body_idx = 0;

    rt->surfaces = surfaces;
    rt->nsurfaces = 6u;
    rt->zones = zone;
    rt->nzones = 1u;
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
    struct osh_step_segment step_segments[16];
    struct ray r;
    size_t n_step_segments;
    int bin;
    double ds;

    build_runtime(&rt, &body, hu, GRID_N_VOX);

    /* Ray passes entirely beside the grid (y outside [0,3]). */
    r = make_x_ray();
    r.p[1] = 10.0; /* far outside */

    ds = dist_voxel_body_rt(&rt, 0, &r, step_segments, 16, &n_step_segments, &bin);
    if (isinf(OSH_GEMCA_INFINITY)) {
        ASSERT_TRUE(isinf(ds));
        ASSERT_TRUE(ds > 0.0);
    } else {
        ASSERT_TRUE(ds == OSH_GEMCA_INFINITY);
    }
    ASSERT_TRUE(n_step_segments == 0);
    ASSERT_TRUE(bin == -1);
}

static void test_uniform_grid_all_same_bin(void) {
    struct osh_gemca_runtime rt;
    struct gemca_rt_body body;
    int16_t hu[GRID_N * GRID_N * GRID_N];
    struct osh_step_segment step_segments[16];
    struct ray r;
    size_t n_step_segments;
    int bin;
    double ds;

    build_runtime(&rt, &body, hu, GRID_N_VOX);
    /* All HU=0 → bin 5. */

    r = make_x_ray();
    ds = dist_voxel_body_rt(&rt, 0, &r, step_segments, 16, &n_step_segments, &bin);

    ASSERT_TRUE(n_step_segments == 1);
    ASSERT_TRUE(bin == 5);

    ASSERT_TRUE(fabs(ds - 1.0) < 1e-9);
    ASSERT_TRUE(fabs(step_segments[0].ds - 1.0) < 1e-9);
    ASSERT_TRUE(step_segments[0].rho == 0.0); /* rho no longer filled here */
}

static void test_bin_change_stops_traversal(void) {
    struct osh_gemca_runtime rt;
    struct gemca_rt_body body;
    int16_t hu[GRID_N * GRID_N * GRID_N];
    struct osh_step_segment step_segments[16];
    struct ray r;
    size_t n_step_segments;
    int bin;
    double ds;

    build_runtime(&rt, &body, hu, GRID_N_VOX);

    /* Make the second X-voxel (ix=1, iy=1, iz=1, flat idx=13) a different bin.
     * HU=1600 → bin 23 (dense bone). */
    hu[13] = 1600;

    r = make_x_ray();
    ds = dist_voxel_body_rt(&rt, 0, &r, step_segments, 16, &n_step_segments, &bin);

    /* Only crosses first voxel (ix=0); stops at the current voxel boundary. */
    ASSERT_TRUE(n_step_segments == 1);
    ASSERT_TRUE(bin == 5);
    ASSERT_TRUE(fabs(ds - 1.0) < 1e-9);
}

static void test_step_segments_cap_limits_output(void) {
    struct osh_gemca_runtime rt;
    struct gemca_rt_body body;
    int16_t hu[GRID_N * GRID_N * GRID_N];
    struct osh_step_segment step_segments[1]; /* cap = 1 */
    struct ray r;
    size_t n_step_segments;
    int bin;
    double ds;

    build_runtime(&rt, &body, hu, GRID_N_VOX);

    r = make_x_ray();
    ds = dist_voxel_body_rt(&rt, 0, &r, step_segments, 1, &n_step_segments, &bin);

    ASSERT_TRUE(n_step_segments == 1);
    ASSERT_TRUE(bin == 5);
    ASSERT_TRUE(fabs(ds - 1.0) < 1e-9); /* only the first voxel's distance */
    ASSERT_TRUE(fabs(step_segments[0].ds - 1.0) < 1e-9);
}

static void test_null_segs_distance_only(void) {
    struct osh_gemca_runtime rt;
    struct gemca_rt_body body;
    int16_t hu[GRID_N * GRID_N * GRID_N];
    struct ray r;
    size_t n_step_segments;
    int bin;
    double ds;

    build_runtime(&rt, &body, hu, GRID_N_VOX);

    r = make_x_ray();
    ds = dist_voxel_body_rt(&rt, 0, &r, NULL, 0, &n_step_segments, &bin);

    ASSERT_TRUE(bin == 5);
    ASSERT_TRUE(n_step_segments == 1);
    ASSERT_TRUE(fabs(ds - 1.0) < 1e-9);
}

static void test_oversized_grid_returns_infinity(void) {
    struct osh_gemca_runtime rt;
    struct gemca_rt_body body;
    int16_t hu[GRID_N * GRID_N * GRID_N];
    struct osh_step_segment step_segments[16];
    struct ray r;
    size_t n_step_segments;
    int bin;
    double ds;

    build_runtime(&rt, &body, hu, GRID_N_VOX);
    body.ct_grid.n[0] = 1024u;
    body.ct_grid.n[1] = 1024u;
    body.ct_grid.n[2] = 1024u; /* 3072 crossings required > _CROSSINGS_CAP */

    r = make_x_ray();
    ds = dist_voxel_body_rt(&rt, 0, &r, step_segments, 16, &n_step_segments, &bin);

    if (isinf(OSH_GEMCA_INFINITY)) {
        ASSERT_TRUE(isinf(ds));
        ASSERT_TRUE(ds > 0.0);
    } else {
        ASSERT_TRUE(ds == OSH_GEMCA_INFINITY);
    }
    ASSERT_TRUE(n_step_segments == 0);
    ASSERT_TRUE(bin == -1);
}

static void test_bin_assigned_for_nonzero_hu(void) {
    struct osh_gemca_runtime rt;
    struct gemca_rt_body body;
    int16_t hu[GRID_N * GRID_N * GRID_N];
    struct osh_step_segment step_segments[16];
    struct ray r;
    size_t n_step_segments;
    size_t i;
    int bin;

    build_runtime(&rt, &body, hu, GRID_N_VOX);

    /* Use a non-zero HU to verify bin assignment. */
    for (i = 0; i < GRID_N_VOX; i++) {
        hu[i] = 200;
    }

    r = make_x_ray();
    dist_voxel_body_rt(&rt, 0, &r, step_segments, 16, &n_step_segments, &bin);

    ASSERT_TRUE(n_step_segments == 1);
    ASSERT_TRUE(bin == (int) s_bin_lut[200 + 1000]);
    ASSERT_TRUE(step_segments[0].rho == 0.0); /* density lives in material_rt, not here */
}

static void test_zone_ref_uses_morton_voxel_index(void) {
    struct osh_gemca_runtime rt;
    struct gemca_rt_body body;
    struct gemca_rt_surface surfaces[6];
    struct gemca_rt_zone zone;
    struct gemca_rt_insn insn;
    int16_t hu[MORTON8_TILE_VOX];
    double x[1] = {1.5};
    double y[1] = {1.5};
    double z[1] = {1.5};
    double ux[1] = {1.0};
    double uy[1] = {0.0};
    double uz[1] = {0.0};
    struct osh_zone_ref zone_ref;
    size_t morton_idx;
    size_t i;

    build_runtime(&rt, &body, hu, MORTON8_TILE_VOX);
    body.ct_grid.tile_order = OSH_VOXEL_ORDER_MORTON8;
    attach_box_voxel_zone(&rt, &body, surfaces, &zone, &insn);

    for (i = 0; i < MORTON8_TILE_VOX; i++) {
        hu[i] = 0;
    }
    morton_idx = osh_voxel_tile_idx(1u, 1u, 1u, 1u, 1u);
    hu[morton_idx] = 1600;

    osh_gemca_runtime_get_zone_ref_batch(&rt, x, y, z, ux, uy, uz, 1u, &zone_ref);

    ASSERT_TRUE(zone_ref.zone_idx == 0u);
    ASSERT_TRUE(zone_ref.has_hu == 1);
    ASSERT_TRUE(zone_ref.hu == 1600);
    ASSERT_TRUE(zone_ref.material_idx == (size_t) s_bin_lut[1600 + 1000]);
}

int main(void) {
    test_ray_misses_returns_infinity();
    test_uniform_grid_all_same_bin();
    test_bin_change_stops_traversal();
    test_step_segments_cap_limits_output();
    test_null_segs_distance_only();
    test_oversized_grid_returns_infinity();
    test_bin_assigned_for_nonzero_hu();
    test_zone_ref_uses_morton_voxel_index();
    return 0;
}
