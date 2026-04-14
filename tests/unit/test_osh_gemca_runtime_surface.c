#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_coord.h"
#include "common/osh_ray.h"
#include "gemca/osh_gemca2_calc_surface.h"
#include "gemca/osh_gemca2_defines.h"
#include "gemca/runtime/osh_gemca_runtime.h"

#define ASSERT_TRUE(cond)                                                                                              \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            fprintf(stderr, "ASSERT FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__);                                 \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

#define N_RANDOM_CASES 64u
#define N_EDGE_CASES 5u
#define N_CASES (N_EDGE_CASES + N_RANDOM_CASES)

static uint64_t test_rng_state = 0x9e3779b97f4a7c15ull;

static double next_unit(void) {
    test_rng_state = test_rng_state * 6364136223846793005ull + 1ull;
    return (double) (test_rng_state >> 11) * (1.0 / 9007199254740992.0);
}

static double next_signed(double scale) {
    return (2.0 * next_unit() - 1.0) * scale;
}

static void append_case(double *x,
                        double *y,
                        double *z,
                        double *ux,
                        double *uy,
                        double *uz,
                        size_t *idx,
                        double px,
                        double py,
                        double pz,
                        double vx,
                        double vy,
                        double vz) {
    x[*idx] = px;
    y[*idx] = py;
    z[*idx] = pz;
    ux[*idx] = vx;
    uy[*idx] = vy;
    uz[*idx] = vz;
    (*idx)++;
}

static void
fill_edge_cases(int type, double *x, double *y, double *z, double *ux, double *uy, double *uz, size_t *idx) {
    switch (type) {
    case OSH_GEMCA_SURF_SPHERE:
        append_case(x, y, z, ux, uy, uz, idx, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
        append_case(x, y, z, ux, uy, uz, idx, 3.0, 0.0, 0.0, -1.0, 0.0, 0.0);
        append_case(x, y, z, ux, uy, uz, idx, 2.0, 0.0, 0.0, -1.0, 0.0, 0.0);
        append_case(x, y, z, ux, uy, uz, idx, 2.0, 0.0, 0.0, 1.0, 0.0, 0.0);
        append_case(x, y, z, ux, uy, uz, idx, 2.0, 0.0, 0.0, 0.0, 1.0, 0.0);
        break;

    case OSH_GEMCA_SURF_ELLIPSOID:
        append_case(x, y, z, ux, uy, uz, idx, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
        append_case(x, y, z, ux, uy, uz, idx, 3.0, 0.0, 0.0, -1.0, 0.0, 0.0);
        append_case(x, y, z, ux, uy, uz, idx, 2.0, 0.0, 0.0, -1.0, 0.0, 0.0);
        append_case(x, y, z, ux, uy, uz, idx, 2.0, 0.0, 0.0, 1.0, 0.0, 0.0);
        append_case(x, y, z, ux, uy, uz, idx, 2.0, 0.0, 0.0, 0.0, 1.0, 0.0);
        break;

    case OSH_GEMCA_SURF_CYLZ:
        append_case(x, y, z, ux, uy, uz, idx, 0.0, 0.0, 5.0, 1.0, 0.0, 0.0);
        append_case(x, y, z, ux, uy, uz, idx, 3.0, 0.0, 0.0, -1.0, 0.0, 0.0);
        append_case(x, y, z, ux, uy, uz, idx, 2.0, 0.0, 7.0, -1.0, 0.0, 0.0);
        append_case(x, y, z, ux, uy, uz, idx, 2.0, 0.0, 7.0, 1.0, 0.0, 0.0);
        append_case(x, y, z, ux, uy, uz, idx, 2.0, 0.0, 7.0, 0.0, 0.0, 1.0);
        break;

    case OSH_GEMCA_SURF_ELLZ:
        append_case(x, y, z, ux, uy, uz, idx, 0.0, 0.0, 5.0, 1.0, 0.0, 0.0);
        append_case(x, y, z, ux, uy, uz, idx, 3.0, 0.0, 0.0, -1.0, 0.0, 0.0);
        append_case(x, y, z, ux, uy, uz, idx, 2.0, 0.0, 7.0, -1.0, 0.0, 0.0);
        append_case(x, y, z, ux, uy, uz, idx, 2.0, 0.0, 7.0, 1.0, 0.0, 0.0);
        append_case(x, y, z, ux, uy, uz, idx, 2.0, 0.0, 7.0, 0.0, 0.0, 1.0);
        break;

    case OSH_GEMCA_SURF_CONE:
        append_case(x, y, z, ux, uy, uz, idx, 0.0, 0.0, 4.0, 1.0, 0.0, 0.0);
        append_case(x, y, z, ux, uy, uz, idx, 3.0, 0.0, 4.0, -1.0, 0.0, 0.0);
        append_case(x, y, z, ux, uy, uz, idx, 2.0, 0.0, 4.0, -1.0, 0.0, 0.0);
        append_case(x, y, z, ux, uy, uz, idx, 2.0, 0.0, 4.0, 1.0, 0.0, 0.0);
        append_case(x, y, z, ux, uy, uz, idx, 2.0, 0.0, 4.0, 0.0, 1.0, 0.0);
        break;

    case OSH_GEMCA_SURF_PLANEX:
        append_case(x, y, z, ux, uy, uz, idx, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0);
        append_case(x, y, z, ux, uy, uz, idx, 3.0, 0.0, 0.0, -1.0, 0.0, 0.0);
        append_case(x, y, z, ux, uy, uz, idx, 2.0, 0.0, 0.0, -1.0, 0.0, 0.0);
        append_case(x, y, z, ux, uy, uz, idx, 2.0, 0.0, 0.0, 1.0, 0.0, 0.0);
        append_case(x, y, z, ux, uy, uz, idx, 2.0, 0.0, 0.0, 0.0, 1.0, 0.0);
        break;

    case OSH_GEMCA_SURF_PLANEY:
        append_case(x, y, z, ux, uy, uz, idx, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0);
        append_case(x, y, z, ux, uy, uz, idx, 0.0, 3.0, 0.0, 0.0, -1.0, 0.0);
        append_case(x, y, z, ux, uy, uz, idx, 0.0, 2.0, 0.0, 0.0, -1.0, 0.0);
        append_case(x, y, z, ux, uy, uz, idx, 0.0, 2.0, 0.0, 0.0, 1.0, 0.0);
        append_case(x, y, z, ux, uy, uz, idx, 0.0, 2.0, 0.0, 1.0, 0.0, 0.0);
        break;

    case OSH_GEMCA_SURF_PLANEZ:
        append_case(x, y, z, ux, uy, uz, idx, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0);
        append_case(x, y, z, ux, uy, uz, idx, 0.0, 0.0, 3.0, 0.0, 0.0, -1.0);
        append_case(x, y, z, ux, uy, uz, idx, 0.0, 0.0, 2.0, 0.0, 0.0, -1.0);
        append_case(x, y, z, ux, uy, uz, idx, 0.0, 0.0, 2.0, 0.0, 0.0, 1.0);
        append_case(x, y, z, ux, uy, uz, idx, 0.0, 0.0, 2.0, 1.0, 0.0, 0.0);
        break;

    case OSH_GEMCA_SURF_PLANE:
        append_case(x, y, z, ux, uy, uz, idx, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
        append_case(x, y, z, ux, uy, uz, idx, 4.0, 0.0, 0.0, -1.0, 0.0, 0.0);
        append_case(x, y, z, ux, uy, uz, idx, 3.0, 0.0, 0.0, -1.0, 0.0, 0.0);
        append_case(x, y, z, ux, uy, uz, idx, 3.0, 0.0, 0.0, 1.0, 0.0, 0.0);
        append_case(x, y, z, ux, uy, uz, idx, 3.0, 0.0, 0.0, 2.0, -1.0, 0.0);
        break;

    default:
        ASSERT_TRUE(0);
    }
}

static void fill_random_cases(double *x, double *y, double *z, double *ux, double *uy, double *uz, size_t *idx) {
    while (*idx < N_CASES) {
        append_case(x,
                    y,
                    z,
                    ux,
                    uy,
                    uz,
                    idx,
                    next_signed(5.0),
                    next_signed(5.0),
                    next_signed(5.0),
                    next_signed(1.0),
                    next_signed(1.0),
                    next_signed(1.0));
    }
}

static void compare_batch_with_scalar(int type, double const *params, int np) {
    struct gemca_rt_surface rt_sf;
    struct surface cold_sf;
    double cold_params[OSH_GEMCA_RT_SURF_NPAR];
    double x[N_CASES];
    double y[N_CASES];
    double z[N_CASES];
    double ux[N_CASES];
    double uy[N_CASES];
    double uz[N_CASES];
    int inside_batch[N_CASES];
    size_t i;
    size_t idx;

    memset(&rt_sf, 0, sizeof(rt_sf));
    memset(&cold_sf, 0, sizeof(cold_sf));
    memset(cold_params, 0, sizeof(cold_params));

    rt_sf.type = type;
    cold_sf.type = type;
    cold_sf.np = np;
    cold_sf.p = cold_params;

    for (i = 0; i < (size_t) np; ++i) {
        rt_sf.p[i] = params[i];
        cold_params[i] = params[i];
    }

    idx = 0u;
    fill_edge_cases(type, x, y, z, ux, uy, uz, &idx);
    ASSERT_TRUE(idx == N_EDGE_CASES);
    fill_random_cases(x, y, z, ux, uy, uz, &idx);
    ASSERT_TRUE(idx == N_CASES);

    osh_gemca_runtime_check_surface_batch(&rt_sf, x, y, z, ux, uy, uz, N_CASES, inside_batch);

    for (i = 0; i < N_CASES; ++i) {
        struct ray r;
        int inside_scalar;

        r.p[0] = x[i];
        r.p[1] = y[i];
        r.p[2] = z[i];
        r.cp[0] = ux[i];
        r.cp[1] = uy[i];
        r.cp[2] = uz[i];
        r.system = OSH_COORD_UNIVERSE;

        inside_scalar = osh_gemca2_check_surface_side(&cold_sf, &r);
        ASSERT_TRUE(inside_batch[i] == inside_scalar);
    }
}

static void test_surface_batch_matches_scalar(void) {
    double sphere[1] = {4.0};
    double ellipsoid[3] = {4.0, 9.0, 16.0};
    double cylz[1] = {4.0};
    double ellz[2] = {4.0, 9.0};
    double cone[2] = {0.0, 0.25};
    double planex[2] = {1.0, -2.0};
    double planey[2] = {1.0, -2.0};
    double planez[2] = {1.0, -2.0};
    double plane[4] = {1.0, 2.0, -1.0, -3.0};

    compare_batch_with_scalar(OSH_GEMCA_SURF_SPHERE, sphere, 1);
    compare_batch_with_scalar(OSH_GEMCA_SURF_ELLIPSOID, ellipsoid, 3);
    compare_batch_with_scalar(OSH_GEMCA_SURF_CYLZ, cylz, 1);
    compare_batch_with_scalar(OSH_GEMCA_SURF_ELLZ, ellz, 2);
    compare_batch_with_scalar(OSH_GEMCA_SURF_CONE, cone, 2);
    compare_batch_with_scalar(OSH_GEMCA_SURF_PLANEX, planex, 2);
    compare_batch_with_scalar(OSH_GEMCA_SURF_PLANEY, planey, 2);
    compare_batch_with_scalar(OSH_GEMCA_SURF_PLANEZ, planez, 2);
    compare_batch_with_scalar(OSH_GEMCA_SURF_PLANE, plane, 4);
}

int main(void) {
    test_surface_batch_matches_scalar();
    return 0;
}
