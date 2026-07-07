/*
 * Unit tests for the scoring-geometry runtime helpers (issue #245 / #251).
 *
 * These helpers were extracted from osh_scoring_step.c into
 * osh_scoring_geometry_runtime.c; real transport only exercises their happy
 * path, so the validation branches are covered here directly over hand-built
 * geometries:
 *   test_mesh_to_grid_ok        — X/Y/Z resolved by label, spacing/origin/n filled
 *   test_mesh_to_grid_errors    — NULL, wrong kind, wrong naxes, missing axis,
 *                                 non-positive nbins, non-positive extent
 *   test_cyl_to_grid_ok         — R/Z resolved by label, [1] slot left unused
 *   test_cyl_to_grid_errors     — NULL, wrong kind, wrong naxes, missing axis,
 *                                 non-positive nbins, non-positive extent
 *   test_zone_bin_index         — hit / miss / NULL / negative-zone paths
 *
 * Pure arithmetic over hand-built runtimes, so identical on every OS.
 */

#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "common/raytrace/osh_raytrace.h"
#include "openshieldhit/status.h"
#include "scoring/runtime/osh_scoring_defs.h"
#include "scoring/runtime/osh_scoring_geometry_runtime.h"
#include "scoring/runtime/osh_scoring_geometry_runtime_internal.h"
#include "test_assert.h"

/* Set one axis descriptor.  Labels are single characters ("X".."R"), so store the
 * char directly rather than pull in a string API. */
static void set_axis(struct osh_scoring_axis_runtime *ax, char label, double lo, double hi, int nbins) {
    memset(ax, 0, sizeof(*ax));
    ax->lo = lo;
    ax->hi = hi;
    ax->nbins = nbins;
    ax->label[0] = label;
    ax->label[1] = '\0';
}

/* ---- Mesh: happy path ----------------------------------------------------- */

static void test_mesh_to_grid_ok(void) {
    struct osh_scoring_axis_runtime axes[3];
    struct osh_scoring_geometry_runtime geo;
    struct osh_raytrace_grid grid;

    /* Declaration order deliberately not X,Y,Z: the builder must resolve by label. */
    set_axis(&axes[0], 'Z', 0.0, 6.0, 3); /* dz = 2 */
    set_axis(&axes[1], 'X', 0.0, 2.0, 2); /* dx = 1 */
    set_axis(&axes[2], 'Y', 0.0, 4.0, 2); /* dy = 2 */
    memset(&geo, 0, sizeof(geo));
    geo.geo_kind = OSH_SCORING_GEO_MESH;
    geo.axes = axes;
    geo.naxes = 3u;

    memset(&grid, 0, sizeof(grid));
    ASSERT_TRUE(osh_scoring_geometry_mesh_to_grid(&geo, &grid) == OSH_OK);

    ASSERT_TRUE(grid.origin[0] == 0.0 && grid.origin[1] == 0.0 && grid.origin[2] == 0.0);
    ASSERT_TRUE(grid.spacing[0] == 1.0 && grid.spacing[1] == 2.0 && grid.spacing[2] == 2.0);
    ASSERT_TRUE(grid.n[0] == 2u && grid.n[1] == 2u && grid.n[2] == 3u);
    ASSERT_TRUE(grid.tile_order == OSH_RAYTRACE_GRID_TILE_ORDER_DEFAULT);
}

/* ---- Mesh: every validation branch --------------------------------------- */

static void test_mesh_to_grid_errors(void) {
    struct osh_scoring_axis_runtime axes[3];
    struct osh_scoring_geometry_runtime geo;
    struct osh_raytrace_grid grid;

    set_axis(&axes[0], 'X', 0.0, 2.0, 2);
    set_axis(&axes[1], 'Y', 0.0, 4.0, 2);
    set_axis(&axes[2], 'Z', 0.0, 6.0, 3);
    memset(&geo, 0, sizeof(geo));
    geo.geo_kind = OSH_SCORING_GEO_MESH;
    geo.axes = axes;
    geo.naxes = 3u;

    /* NULL arguments. */
    ASSERT_TRUE(osh_scoring_geometry_mesh_to_grid(NULL, &grid) == OSH_EINVAL);
    ASSERT_TRUE(osh_scoring_geometry_mesh_to_grid(&geo, NULL) == OSH_EINVAL);

    /* Wrong geometry kind is ENOTSUP, not EINVAL. */
    geo.geo_kind = OSH_SCORING_GEO_CYL;
    ASSERT_TRUE(osh_scoring_geometry_mesh_to_grid(&geo, &grid) == OSH_ENOTSUP);
    geo.geo_kind = OSH_SCORING_GEO_MESH;

    /* A mesh must have exactly three axes. */
    geo.naxes = 2u;
    ASSERT_TRUE(osh_scoring_geometry_mesh_to_grid(&geo, &grid) == OSH_EINVAL);
    geo.naxes = 3u;

    /* Missing a required label (Z relabelled) fails the label lookup. */
    axes[2].label[0] = 'Q';
    ASSERT_TRUE(osh_scoring_geometry_mesh_to_grid(&geo, &grid) == OSH_EINVAL);
    axes[2].label[0] = 'Z';

    /* Non-positive bin count. */
    axes[1].nbins = 0;
    ASSERT_TRUE(osh_scoring_geometry_mesh_to_grid(&geo, &grid) == OSH_EINVAL);
    axes[1].nbins = 2;

    /* Non-positive extent (lo == hi) makes the spacing zero. */
    axes[0].hi = axes[0].lo;
    ASSERT_TRUE(osh_scoring_geometry_mesh_to_grid(&geo, &grid) == OSH_EINVAL);
}

/* ---- Cyl: happy path ------------------------------------------------------ */

static void test_cyl_to_grid_ok(void) {
    struct osh_scoring_axis_runtime axes[2];
    struct osh_scoring_geometry_runtime geo;
    struct osh_raytrace_grid grid;

    /* Z declared before R: order must not matter. */
    set_axis(&axes[0], 'Z', 0.0, 10.0, 5); /* dz = 2 */
    set_axis(&axes[1], 'R', 0.0, 3.0, 3);  /* dr = 1 */
    memset(&geo, 0, sizeof(geo));
    geo.geo_kind = OSH_SCORING_GEO_CYL;
    geo.axes = axes;
    geo.naxes = 2u;

    memset(&grid, 0, sizeof(grid));
    ASSERT_TRUE(osh_scoring_geometry_cyl_to_grid(&geo, &grid) == OSH_OK);

    /* [0] = radial, [2] = axial, [1] unused (0/0/1). */
    ASSERT_TRUE(grid.origin[0] == 0.0 && grid.origin[1] == 0.0 && grid.origin[2] == 0.0);
    ASSERT_TRUE(grid.spacing[0] == 1.0 && grid.spacing[1] == 0.0 && grid.spacing[2] == 2.0);
    ASSERT_TRUE(grid.n[0] == 3u && grid.n[1] == 1u && grid.n[2] == 5u);
    ASSERT_TRUE(grid.tile_order == OSH_RAYTRACE_GRID_TILE_ORDER_DEFAULT);
}

/* ---- Cyl: every validation branch ---------------------------------------- */

static void test_cyl_to_grid_errors(void) {
    struct osh_scoring_axis_runtime axes[2];
    struct osh_scoring_geometry_runtime geo;
    struct osh_raytrace_grid grid;

    set_axis(&axes[0], 'R', 0.0, 3.0, 3);
    set_axis(&axes[1], 'Z', 0.0, 10.0, 5);
    memset(&geo, 0, sizeof(geo));
    geo.geo_kind = OSH_SCORING_GEO_CYL;
    geo.axes = axes;
    geo.naxes = 2u;

    ASSERT_TRUE(osh_scoring_geometry_cyl_to_grid(NULL, &grid) == OSH_EINVAL);
    ASSERT_TRUE(osh_scoring_geometry_cyl_to_grid(&geo, NULL) == OSH_EINVAL);

    geo.geo_kind = OSH_SCORING_GEO_MESH;
    ASSERT_TRUE(osh_scoring_geometry_cyl_to_grid(&geo, &grid) == OSH_ENOTSUP);
    geo.geo_kind = OSH_SCORING_GEO_CYL;

    geo.naxes = 3u;
    ASSERT_TRUE(osh_scoring_geometry_cyl_to_grid(&geo, &grid) == OSH_EINVAL);
    geo.naxes = 2u;

    axes[1].label[0] = 'Q'; /* Z relabelled -> label lookup fails */
    ASSERT_TRUE(osh_scoring_geometry_cyl_to_grid(&geo, &grid) == OSH_EINVAL);
    axes[1].label[0] = 'Z';

    axes[0].nbins = 0;
    ASSERT_TRUE(osh_scoring_geometry_cyl_to_grid(&geo, &grid) == OSH_EINVAL);
    axes[0].nbins = 3;

    axes[0].hi = axes[0].lo; /* dr = 0 */
    ASSERT_TRUE(osh_scoring_geometry_cyl_to_grid(&geo, &grid) == OSH_EINVAL);
}

/* ---- Zone bin index ------------------------------------------------------- */

static void test_zone_bin_index(void) {
    struct osh_scoring_geometry_runtime geo;
    size_t zones[3];
    size_t idx;

    zones[0] = 5u;
    zones[1] = 7u;
    zones[2] = 9u;
    memset(&geo, 0, sizeof(geo));
    geo.geo_kind = OSH_SCORING_GEO_ZONE;
    geo.zone_indices = zones;
    geo.nzone_indices = 3u;

    /* Hits map the transport zone id to its dense bin position. */
    idx = 99u;
    ASSERT_TRUE(osh_scoring_geometry_zone_bin_index(&geo, 5, &idx) == 1 && idx == 0u);
    idx = 99u;
    ASSERT_TRUE(osh_scoring_geometry_zone_bin_index(&geo, 7, &idx) == 1 && idx == 1u);
    idx = 99u;
    ASSERT_TRUE(osh_scoring_geometry_zone_bin_index(&geo, 9, &idx) == 1 && idx == 2u);

    /* A zone not in the selected list is a miss; idx is left untouched. */
    idx = 99u;
    ASSERT_TRUE(osh_scoring_geometry_zone_bin_index(&geo, 6, &idx) == 0 && idx == 99u);

    /* Guard inputs: NULL geo, NULL out, negative zone all return 0. */
    ASSERT_TRUE(osh_scoring_geometry_zone_bin_index(NULL, 5, &idx) == 0);
    ASSERT_TRUE(osh_scoring_geometry_zone_bin_index(&geo, 5, NULL) == 0);
    ASSERT_TRUE(osh_scoring_geometry_zone_bin_index(&geo, -1, &idx) == 0);
}

int main(void) {
    test_mesh_to_grid_ok();
    test_mesh_to_grid_errors();
    test_cyl_to_grid_ok();
    test_cyl_to_grid_errors();
    test_zone_bin_index();
    printf("All osh_scoring_geometry_runtime tests passed.\n");
    return 0;
}
