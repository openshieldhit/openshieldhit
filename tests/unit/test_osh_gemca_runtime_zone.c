#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_geometry_parse.h"
#include "common/osh_coord.h"
#include "common/osh_ray.h"
#include "gemca/osh_gemca2.h"
#include "gemca/runtime/osh_gemca_runtime.h"
#include "openshieldhit/geometry.h"

#define ASSERT_TRUE(cond)                                                                                              \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            fprintf(stderr, "ASSERT FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__);                                 \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

#define N_ZONE_CASES 773u

static uint64_t test_rng_state = 0xa0761d6478bd642full;

static double next_unit(void) {
    test_rng_state = test_rng_state * 6364136223846793005ull + 1ull;
    return (double) (test_rng_state >> 11) * (1.0 / 9007199254740992.0);
}

static double next_signed(double scale) {
    return (2.0 * next_unit() - 1.0) * scale;
}

static void fill_zone_cases(double *x, double *y, double *z, double *ux, double *uy, double *uz) {
    size_t i;

    x[0] = 0.0;
    y[0] = 0.0;
    z[0] = 0.0;
    ux[0] = 1.0;
    uy[0] = 0.0;
    uz[0] = 0.0;
    x[1] = 1.0;
    y[1] = 0.0;
    z[1] = 0.0;
    ux[1] = -1.0;
    uy[1] = 0.0;
    uz[1] = 0.0;
    x[2] = 0.0;
    y[2] = 1.0;
    z[2] = 0.0;
    ux[2] = 0.0;
    uy[2] = -1.0;
    uz[2] = 0.0;
    x[3] = 0.0;
    y[3] = 0.0;
    z[3] = 1.0;
    ux[3] = 0.0;
    uy[3] = 0.0;
    uz[3] = -1.0;
    x[4] = 5.0;
    y[4] = 5.0;
    z[4] = 5.0;
    ux[4] = -1.0;
    uy[4] = -1.0;
    uz[4] = -1.0;

    for (i = 5u; i < N_ZONE_CASES; ++i) {
        x[i] = next_signed(20.0);
        y[i] = next_signed(20.0);
        z[i] = next_signed(20.0);
        ux[i] = next_signed(1.0);
        uy[i] = next_signed(1.0);
        uz[i] = next_signed(1.0);
    }
}

static void test_zone_batch_matches_scalar(void) {
    struct osh_geometry_workspace *geom;
    struct osh_gemca_prepared *g;
    struct osh_gemca_runtime rt;
    struct ray r;
    char geo_path[512];
    double x[N_ZONE_CASES];
    double y[N_ZONE_CASES];
    double z[N_ZONE_CASES];
    double ux[N_ZONE_CASES];
    double uy[N_ZONE_CASES];
    double uz[N_ZONE_CASES];
    size_t zones_batch[N_ZONE_CASES];
    size_t i;

    geom = NULL;
    g = NULL;
    memset(&rt, 0, sizeof(rt));

    snprintf(geo_path, sizeof(geo_path), "%s/examples/01_sdl_viewer/geo_RCC03.dat", OSH_PROJECT_SOURCE_DIR);

    ASSERT_TRUE(osh_geometry_parse_file(geo_path, &geom) == OSH_OK);
    ASSERT_TRUE(geom != NULL);
    ASSERT_TRUE(osh_geometry_workspace_prepare(geom) == 0);
    g = geom->prepared;
    ASSERT_TRUE(osh_gemca_runtime_setup(g, &rt) == OSH_OK);

    fill_zone_cases(x, y, z, ux, uy, uz);
    osh_gemca_runtime_get_zone_batch(&rt, x, y, z, ux, uy, uz, N_ZONE_CASES, zones_batch);

    r.system = OSH_COORD_UNIVERSE;
    for (i = 0; i < N_ZONE_CASES; ++i) {
        size_t zone_scalar;

        r.p[0] = x[i];
        r.p[1] = y[i];
        r.p[2] = z[i];
        r.cp[0] = ux[i];
        r.cp[1] = uy[i];
        r.cp[2] = uz[i];

        zone_scalar = osh_gemca_get_zone_index(g, &r);
        ASSERT_TRUE(zones_batch[i] == zone_scalar);
    }

    osh_gemca_runtime_free(&rt);
    osh_geometry_workspace_free(geom);
}

int main(void) {
    test_zone_batch_matches_scalar();
    return 0;
}
