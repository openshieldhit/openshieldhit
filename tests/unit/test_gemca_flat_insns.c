/*
 * Unit test for the GPU-portable flat instruction store added to
 * struct osh_gemca_runtime by M1 (see src/gemca/runtime/README.md
 * "GPU migration path").
 *
 * For every fixture geometry we:
 *   - compile the cold workspace into a runtime,
 *   - assert rt->insn_begin[] is monotonically non-decreasing with the
 *     expected per-zone step sizes (insn_begin[j+1] - insn_begin[j]
 *     == zones[j].ninsns),
 *   - assert rt->ninsns_flat == insn_begin[nzones] == sum of all
 *     zones[j].ninsns,
 *   - assert the flat array's contents equal the per-zone arrays
 *     element-by-element for every zone.
 *
 * The CPU path keeps using zones[j].insns; this test guards the additive
 * GPU-portable copy against regressions in setup_zones().
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_app_osh.h"
#include "common/osh_coord.h"
#include "common/osh_ray.h"
#include "gemca/osh_gemca2.h"
#include "gemca/runtime/osh_gemca_runtime.h"
#include "openshieldhit/geometry.h"
#include "openshieldhit/voxel.h"

#define ASSERT_TRUE(cond)                                                                                              \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            fprintf(stderr, "ASSERT FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__);                                 \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

/*
 * A small table of fixture geometries scattered around the source tree.
 * Each entry is a path relative to OSH_PROJECT_SOURCE_DIR.  We compile each
 * one to a runtime and exercise the flat instruction-store invariants.
 */
static char const *FIXTURE_GEOMETRIES[] = {
    "examples/01_sdl_viewer/geo_RCC03.dat",
    "examples/01_sdl_viewer/geo.dat",
    "examples/02_bnct/geo_cell.dat",
    "tests/fixtures/test01/geo.dat",
};

static size_t const N_FIXTURES = sizeof(FIXTURE_GEOMETRIES) / sizeof(FIXTURE_GEOMETRIES[0]);

static void test_flat_insn_store_invariants(void) {
    size_t f;

    for (f = 0; f < N_FIXTURES; ++f) {
        struct osh_geometry_workspace *geom;
        struct osh_gemca_prepared *g;
        struct osh_gemca_runtime rt;
        char geo_path[512];
        size_t j;
        int expected_total;
        int prev_offset;

        geom = NULL;
        g = NULL;
        memset(&rt, 0, sizeof(rt));

        snprintf(geo_path, sizeof(geo_path), "%s/%s", OSH_PROJECT_SOURCE_DIR, FIXTURE_GEOMETRIES[f]);

        ASSERT_TRUE(osh_geometry_setup_from_path(geo_path, NULL, &geom) == OSH_OK);
        g = geom->prepared;
        ASSERT_TRUE(osh_gemca_compile(g, OSH_HU_TABLE_NONE, 0u, NULL, &rt) == OSH_OK);

        /* Sum each zone's ninsns; the flat buffer must hold exactly that many
         * instructions and insn_begin[nzones] must equal the same sum. */
        expected_total = 0;
        for (j = 0; j < rt.nzones; ++j) {
            expected_total += rt.zones[j].ninsns;
        }

        ASSERT_TRUE(rt.insns_flat != NULL);
        ASSERT_TRUE(rt.insn_begin != NULL);
        ASSERT_TRUE(rt.ninsns_flat == (size_t) expected_total);
        ASSERT_TRUE(rt.insn_begin[rt.nzones] == expected_total);

        /* insn_begin[] must be monotonically non-decreasing and each
         * per-zone slice must match the contents of zones[j].insns[]
         * element-by-element. */
        prev_offset = 0;
        for (j = 0; j < rt.nzones; ++j) {
            int const begin = rt.insn_begin[j];
            int const end = rt.insn_begin[j + 1u];
            int const len = end - begin;

            ASSERT_TRUE(begin == prev_offset);
            ASSERT_TRUE(len == rt.zones[j].ninsns);
            ASSERT_TRUE(len >= 0);

            if (len > 0) {
                ASSERT_TRUE(memcmp(&rt.insns_flat[begin],
                                   rt.zones[j].insns,
                                   (size_t) len * sizeof(struct gemca_rt_insn)) == 0);
            }

            prev_offset = end;
        }

        ASSERT_TRUE(prev_offset == expected_total);

        osh_gemca_runtime_free(&rt);
        osh_geometry_workspace_free(geom);
    }
}

int main(void) {
    test_flat_insn_store_invariants();
    return 0;
}
