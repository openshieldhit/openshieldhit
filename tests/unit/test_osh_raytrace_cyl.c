#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/raytrace/osh_raytrace_cyl.h"
#include "test_assert.h"

/* Grid with r_min=0, dr=1, nr=5, z_min=0, dz=1, nz=5. */
static struct osh_raytrace_grid make_grid(void) {
    struct osh_raytrace_grid g;
    memset(&g, 0, sizeof(g));
    g.origin[0] = 0.0;
    g.spacing[0] = 1.0;
    g.n[0] = 5u; /* R */
    g.origin[1] = 0.0;
    g.spacing[1] = 0.0;
    g.n[1] = 1u; /* unused */
    g.origin[2] = 0.0;
    g.spacing[2] = 1.0;
    g.n[2] = 5u; /* Z */
    return g;
}

/* CYL max crossings = 2*nr + nz. */
#define MAX_CROSS 15u

static void test_axial_step_only_z_crossings(void) {
    /* Step parallel to Z axis: no R crossings expected. */
    struct osh_raytrace_grid grid = make_grid();
    /* p at (1.5, 0, 0.5): r=1.5 → ir=1;  z=0.5 → iz=0. */
    double p[3] = {1.5, 0.0, 0.5};
    double v[3] = {0.0, 0.0, 1.0};
    double ds = 3.0; /* goes from z=0.5 to z=3.5; crosses z=1,2,3 */
    struct osh_voxel_crossing crossings[MAX_CROSS];
    size_t n = 0u;
    int rc;
    size_t j;
    double sum = 0.0;

    rc = osh_raytrace_cyl_traverse(&grid, p, v, ds, crossings, &n);
    ASSERT_TRUE(rc == 1);
    ASSERT_TRUE(n == 4u); /* iz=0(0.5cm), iz=1(1cm), iz=2(1cm), iz=3(0.5cm) */

    /* All crossings must be in ir=1. */
    for (j = 0; j < n; ++j) {
        size_t ir = crossings[j].idx % grid.n[0];
        ASSERT_TRUE(ir == 1u);
        ASSERT_TRUE(crossings[j].path_len > 0.0);
        sum += crossings[j].path_len;
    }
    /* Total scored path must equal ds. */
    ASSERT_TRUE(fabs(sum - ds) < 1.0e-10);

    /* Verify iz sequence: 0, 1, 2, 3. */
    ASSERT_TRUE(crossings[0].idx == 1u + 5u * 0u);
    ASSERT_TRUE(crossings[1].idx == 1u + 5u * 1u);
    ASSERT_TRUE(crossings[2].idx == 1u + 5u * 2u);
    ASSERT_TRUE(crossings[3].idx == 1u + 5u * 3u);
}

static void test_diagonal_step_path_len_sum(void) {
    /* Step moves both outward in R and along Z: verify sum(path_len)==clipped ds. */
    struct osh_raytrace_grid grid = make_grid();
    /* p at (0.5, 0, 0.5): r=0.5→ir=0, z=0.5→iz=0.
     * v moves radially outward + along Z: (1/sqrt2, 0, 1/sqrt2). */
    double inv = 1.0 / sqrt(2.0);
    double p[3] = {0.5, 0.0, 0.5};
    double v[3] = {inv, 0.0, inv};
    double ds = 6.0; /* long enough to cross multiple R shells and Z slabs */
    struct osh_voxel_crossing crossings[MAX_CROSS];
    size_t n = 0u;
    int rc;
    size_t j;
    double sum = 0.0;

    rc = osh_raytrace_cyl_traverse(&grid, p, v, ds, crossings, &n);
    ASSERT_TRUE(rc == 1);
    ASSERT_TRUE(n > 0u);

    for (j = 0; j < n; ++j) {
        ASSERT_TRUE(crossings[j].path_len > 0.0);
        ASSERT_TRUE(crossings[j].idx < grid.n[0] * grid.n[2]);
        sum += crossings[j].path_len;
    }
    /* Grid exits at r=5 or z=5; path is clipped. At ds=6, the Z exit happens
     * first (z=0.5+6/sqrt(2)≈4.74 < 5) and r=0.5+6/sqrt(2)≈4.74 < 5, so the
     * full step stays inside.  Sum should equal ds. */
    ASSERT_TRUE(fabs(sum - ds) < 1.0e-9);
}

static void test_miss_outside_r_max(void) {
    struct osh_raytrace_grid grid = make_grid();
    /* Particle far outside r_max=5, moving even further away. */
    double p[3] = {10.0, 0.0, 2.5};
    double v[3] = {1.0, 0.0, 0.0};
    struct osh_voxel_crossing crossings[MAX_CROSS];
    size_t n = 99u;
    int rc;

    rc = osh_raytrace_cyl_traverse(&grid, p, v, 1.0, crossings, &n);
    ASSERT_TRUE(rc == 0);
    ASSERT_TRUE(n == 0u);
}

static void test_miss_inside_hollow_hole(void) {
    /* Grid with non-zero r_min (hollow cylinder). */
    struct osh_raytrace_grid grid = make_grid();
    grid.origin[0] = 5.0; /* r_min=5 */

    /* p at r=1 (inside hole), step along Z — should miss. */
    double p[3] = {1.0, 0.0, 2.0};
    double v[3] = {0.0, 0.0, 1.0};
    struct osh_voxel_crossing crossings[MAX_CROSS];
    size_t n = 99u;
    int rc;

    rc = osh_raytrace_cyl_traverse(&grid, p, v, 2.0, crossings, &n);
    ASSERT_TRUE(rc == 0);
    ASSERT_TRUE(n == 0u);
}

static void test_miss_outside_z_range(void) {
    struct osh_raytrace_grid grid = make_grid();
    /* Both start and end are below z_min=0. */
    double p[3] = {1.5, 0.0, -3.0};
    double v[3] = {0.0, 0.0, -1.0};
    struct osh_voxel_crossing crossings[MAX_CROSS];
    size_t n = 99u;
    int rc;

    rc = osh_raytrace_cyl_traverse(&grid, p, v, 2.0, crossings, &n);
    ASSERT_TRUE(rc == 0);
    ASSERT_TRUE(n == 0u);
}

static void test_vol_inv_initialised_zero(void) {
    /* Traversal must NOT fill vol_inv — it must remain 0.0 (as written). */
    struct osh_raytrace_grid grid = make_grid();
    double p[3] = {1.5, 0.0, 0.5};
    double v[3] = {0.0, 0.0, 1.0};
    struct osh_voxel_crossing crossings[MAX_CROSS];
    size_t n = 0u;
    size_t j;

    memset(crossings, 0x55, sizeof(crossings)); /* poison with non-zero */
    osh_raytrace_cyl_traverse(&grid, p, v, 2.0, crossings, &n);
    ASSERT_TRUE(n > 0u);
    for (j = 0; j < n; ++j) {
        ASSERT_TRUE(crossings[j].vol_inv == 0.0);
    }
}

static void test_turning_point_through_axis(void) {
    /* Ray whose R decreases to a minimum then increases: the genuinely hard case.
     * Grid: r_min=0, dr=1, nr=5, z_min=0, dz=1, nz=5 (from make_grid).
     * p=(-3, 0.5, 2.5): r=sqrt(9.25)≈3.04 → ir=3, iz=2.
     * v=(1,0,0): sweeps x from -3 to +3 (ds=6), minimum r=0.5 at x=0.
     * The same R shells ir=3,2,1 are each crossed TWICE (in then out).
     * Expected ir sequence: 3,2,1,0,1,2,3 (all iz=2). sum(path_len)==ds. */
    struct osh_raytrace_grid grid = make_grid();
    struct osh_voxel_crossing crossings[MAX_CROSS];
    double p[3] = {-3.0, 0.5, 2.5};
    double v[3] = {1.0, 0.0, 0.0};
    double ds = 6.0;
    size_t n = 0u;
    int rc;
    size_t j;
    double sum = 0.0;
    /* Expected ir sequence at iz=2 (idx = ir + 5*2 = ir + 10) */
    size_t expected_idx[7] = {13u, 12u, 11u, 10u, 11u, 12u, 13u};

    rc = osh_raytrace_cyl_traverse(&grid, p, v, ds, crossings, &n);
    ASSERT_TRUE(rc == 1);
    ASSERT_TRUE(n == 7u);

    for (j = 0; j < n; ++j) {
        ASSERT_TRUE(crossings[j].path_len > 0.0);
        ASSERT_TRUE(crossings[j].idx == expected_idx[j]);
        sum += crossings[j].path_len;
    }
    ASSERT_TRUE(fabs(sum - ds) < 1.0e-9);
}

static void test_hollow_start_inside_hole(void) {
    /* Hollow cylinder (r_min=2): particle starts inside the hole and exits radially.
     * Grid: r=[2,5] (nr=3, dr=1), z=[0,5] (nz=5, dz=1).
     * p=(0.5,0,2.5): r=0.5 < r_min=2, inside hole, iz=2.
     * v=(1,0,0): moving outward.  ds=4.0 → x goes from 0.5 to 4.5.
     * Path inside hole [x=0.5 to x=2] must NOT be scored.
     * Valid crossings: ir=0 [x=2..3], ir=1 [x=3..4], partial ir=2 [x=4..4.5].
     * Sum of path_len must equal 2.5 (not 4.0). */
    struct osh_raytrace_grid grid = make_grid();
    struct osh_voxel_crossing crossings[MAX_CROSS];
    double p[3] = {0.5, 0.0, 2.5};
    double v[3] = {1.0, 0.0, 0.0};
    double ds = 4.0;
    size_t n = 0u;
    int rc;
    size_t j;
    double sum = 0.0;

    grid.origin[0] = 2.0; /* r_min=2 */
    rc = osh_raytrace_cyl_traverse(&grid, p, v, ds, crossings, &n);
    ASSERT_TRUE(rc == 1);
    ASSERT_TRUE(n > 0u);
    for (j = 0; j < n; ++j) {
        ASSERT_TRUE(crossings[j].path_len > 0.0);
        sum += crossings[j].path_len;
    }
    ASSERT_TRUE(fabs(sum - 2.5) < 1.0e-9);
}

static void test_hollow_through_hole(void) {
    /* Hollow cylinder (r_min=2): particle traverses the hole mid-step.
     * Grid: r=[2,5] (nr=3, dr=1), z=[0,5] (nz=5, dz=1).
     * p=(-4.5,0,2.5), v=(1,0,0), ds=9.0 → x from -4.5 to 4.5.
     * Crosses: r=4 at x=-4 (ir=2→1), r=3 at x=-3 (ir=1→0), r=2 at x=-2 (enters hole),
     *           r=2 at x=2 (exits hole, ir=0), r=3 at x=3 (ir=0→1), r=4 at x=4 (ir=1→2).
     * Scored lengths: ir=2: 0.5+0.5=1, ir=1: 1+1=2, ir=0: 1+1=2 → total 5.0 cm. */
    struct osh_raytrace_grid grid = make_grid();
    struct osh_voxel_crossing crossings[MAX_CROSS];
    double p[3] = {-4.5, 0.0, 2.5};
    double v[3] = {1.0, 0.0, 0.0};
    double ds = 9.0;
    size_t n = 0u;
    int rc;
    size_t j;
    double sum = 0.0;

    grid.origin[0] = 2.0; /* r_min=2 */
    rc = osh_raytrace_cyl_traverse(&grid, p, v, ds, crossings, &n);
    ASSERT_TRUE(rc == 1);
    ASSERT_TRUE(n > 0u);
    for (j = 0; j < n; ++j) {
        ASSERT_TRUE(crossings[j].path_len > 0.0);
        sum += crossings[j].path_len;
    }
    /* 4 cm of the step is inside the hollow hole and must not be scored */
    ASSERT_TRUE(fabs(sum - 5.0) < 1.0e-9);
}

int main(int argc, char **argv) {
    static struct {
        char const *name;
        void (*fn)(void);
    } cases[] = {
        {"test_axial_step_only_z_crossings", test_axial_step_only_z_crossings},
        {"test_diagonal_step_path_len_sum", test_diagonal_step_path_len_sum},
        {"test_miss_outside_r_max", test_miss_outside_r_max},
        {"test_miss_inside_hollow_hole", test_miss_inside_hollow_hole},
        {"test_miss_outside_z_range", test_miss_outside_z_range},
        {"test_vol_inv_initialised_zero", test_vol_inv_initialised_zero},
        {"test_turning_point_through_axis", test_turning_point_through_axis},
        {"test_hollow_start_inside_hole", test_hollow_start_inside_hole},
        {"test_hollow_through_hole", test_hollow_through_hole},
    };

    size_t ncases = sizeof(cases) / sizeof(cases[0]);
    size_t i;

    if (argc > 1) {
        /* Run a single named test (used by register_named_c_tests CTest adapter). */
        for (i = 0; i < ncases; ++i) {
            if (strcmp(argv[1], cases[i].name) == 0) {
                cases[i].fn();
                return 0;
            }
        }
        fprintf(stderr, "Unknown test: %s\n", argv[1]);
        return 1;
    }

    for (i = 0; i < ncases; ++i) {
        cases[i].fn();
    }
    return 0;
}
