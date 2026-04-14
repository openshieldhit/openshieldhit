#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_coord.h"
#include "common/osh_ray.h"
#include "gemca/osh_gemca2.h"
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

#define N_RANDOM_CASES 96u
#define N_EDGE_CASES 6u
#define N_CASES (N_EDGE_CASES + N_RANDOM_CASES)

static uint64_t test_rng_state = 0xd1b54a32d192ed03ull;

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

static enum osh_status reference_transform_to_local(struct body const *b, struct ray const *r, struct ray *tr) {
    int i;
    int j;

    for (i = 0; i < 3; ++i) {
        tr->p[i] = r->p[i];
        tr->cp[i] = r->cp[i];
    }
    tr->system = (unsigned char) b->coord;

    switch (b->coord) {
    case OSH_COORD_UNIVERSE:
        break;

    case OSH_COORD_BCALIGN:
        for (i = 0; i < 3; ++i) {
            j = i * 4;
            tr->p[i] = r->p[i] + b->t[j + 3];
            tr->cp[i] = r->cp[i];
        }
        break;

    case OSH_COORD_BZALIGN:
        osh_ray_transform(r, tr, b->t);
        break;

    default:
        return OSH_ENOTSUP;
    }

    return OSH_OK;
}

static int reference_in_body(struct body const *b, struct ray const *r) {
    struct ray tr;
    int i;

    if (reference_transform_to_local(b, r, &tr) != OSH_OK) {
        return 0;
    }

    for (i = 0; i < b->nsurfs; ++i) {
        if (!osh_gemca2_check_surface_side(b->surfs[i], &tr)) {
            return 0;
        }
    }

    return 1;
}

static void init_cold_body_from_runtime(struct body *cold,
                                        struct gemca_rt_body const *rt_body,
                                        struct gemca_rt_surface const *rt_surfaces,
                                        struct surface *cold_surfaces,
                                        struct surface **cold_surface_ptrs,
                                        double params[][OSH_GEMCA_RT_SURF_NPAR]) {
    size_t i;
    size_t j;

    memset(cold, 0, sizeof(*cold));
    memcpy(cold->t, rt_body->t, sizeof(cold->t));
    cold->coord = rt_body->coord;
    cold->type = rt_body->type;
    cold->nsurfs = rt_body->nsurfs;
    cold->surfs = cold_surface_ptrs;

    for (i = 0; i < (size_t) rt_body->nsurfs; ++i) {
        struct gemca_rt_surface const *src;

        src = &rt_surfaces[rt_body->surf_begin + i];
        memset(&cold_surfaces[i], 0, sizeof(cold_surfaces[i]));
        cold_surfaces[i].type = src->type;
        cold_surfaces[i].np = OSH_GEMCA_RT_SURF_NPAR;
        cold_surfaces[i].p = params[i];
        for (j = 0; j < OSH_GEMCA_RT_SURF_NPAR; ++j) {
            params[i][j] = src->p[j];
        }
        cold_surface_ptrs[i] = &cold_surfaces[i];
    }
}

static void fill_sphere_cases(double *x,
                              double *y,
                              double *z,
                              double *ux,
                              double *uy,
                              double *uz) {
    size_t idx;

    idx = 0u;
    append_case(x, y, z, ux, uy, uz, &idx, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    append_case(x, y, z, ux, uy, uz, &idx, 3.0, 0.0, 0.0, -1.0, 0.0, 0.0);
    append_case(x, y, z, ux, uy, uz, &idx, 2.0, 0.0, 0.0, -1.0, 0.0, 0.0);
    append_case(x, y, z, ux, uy, uz, &idx, 2.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    append_case(x, y, z, ux, uy, uz, &idx, 0.0, 1.0, 1.0, 0.0, -1.0, 0.0);
    append_case(x, y, z, ux, uy, uz, &idx, -1.0, -1.0, -1.0, 0.2, 0.3, -0.4);

    while (idx < N_CASES) {
        append_case(x, y, z, ux, uy, uz, &idx,
                    next_signed(4.0), next_signed(4.0), next_signed(4.0),
                    next_signed(1.0), next_signed(1.0), next_signed(1.0));
    }
}

static void fill_bcalign_box_cases(double *x,
                                   double *y,
                                   double *z,
                                   double *ux,
                                   double *uy,
                                   double *uz) {
    size_t idx;
    double const tx = 1.5;
    double const ty = -0.5;
    double const tz = 2.0;

    idx = 0u;
    append_case(x, y, z, ux, uy, uz, &idx, -tx, -ty, -tz, 1.0, 0.0, 0.0);            /* local (0,0,0) */
    append_case(x, y, z, ux, uy, uz, &idx, 1.0 - tx, -ty, -tz, -1.0, 0.0, 0.0);      /* local on +X face */
    append_case(x, y, z, ux, uy, uz, &idx, -1.0 - tx, -ty, -tz, 1.0, 0.0, 0.0);      /* local on -X face */
    append_case(x, y, z, ux, uy, uz, &idx, 1.2 - tx, -ty, -tz, -1.0, 0.0, 0.0);      /* outside +X */
    append_case(x, y, z, ux, uy, uz, &idx, -tx, 1.1 - ty, -tz, 0.0, -1.0, 0.0);      /* outside +Y */
    append_case(x, y, z, ux, uy, uz, &idx, -tx, -ty, -1.0 - tz, 0.0, 0.0, 1.0);      /* local on -Z face */

    while (idx < N_CASES) {
        append_case(x, y, z, ux, uy, uz, &idx,
                    next_signed(4.0), next_signed(4.0), next_signed(4.0),
                    next_signed(1.0), next_signed(1.0), next_signed(1.0));
    }
}

static void fill_bzalign_capped_cyl_cases(double *x,
                                          double *y,
                                          double *z,
                                          double *ux,
                                          double *uy,
                                          double *uz) {
    size_t idx;

    idx = 0u;

    /* Inverse map for the chosen BZALIGN matrix:
     * local_x = world_y - 1
     * local_y = -world_x + 2
     * local_z = world_z - 0.5
     */
    append_case(x, y, z, ux, uy, uz, &idx, 2.0, 1.0, 1.5, 1.0, 0.0, 0.0);            /* local (0,0,1) */
    append_case(x, y, z, ux, uy, uz, &idx, 2.0, 2.0, 1.5, 0.0, -1.0, 0.0);           /* local (1,0,1) on cyl */
    append_case(x, y, z, ux, uy, uz, &idx, 2.0, 1.0, 0.5, 0.0, 0.0, 1.0);            /* local z=0 cap */
    append_case(x, y, z, ux, uy, uz, &idx, 2.0, 1.0, 2.5, 0.0, 0.0, -1.0);           /* local z=2 cap */
    append_case(x, y, z, ux, uy, uz, &idx, 2.0, 2.2, 1.5, 0.0, -1.0, 0.0);           /* outside cyl */
    append_case(x, y, z, ux, uy, uz, &idx, 2.0, 1.0, 2.8, 0.0, 0.0, -1.0);           /* above cap */

    while (idx < N_CASES) {
        append_case(x, y, z, ux, uy, uz, &idx,
                    next_signed(4.0), next_signed(4.0), next_signed(4.0),
                    next_signed(1.0), next_signed(1.0), next_signed(1.0));
    }
}

static void compare_body_batch_with_scalar(struct gemca_runtime const *rt,
                                           size_t body_idx,
                                           struct body const *cold_body,
                                           void (*fill_cases)(double *, double *, double *, double *, double *, double *)) {
    double x[N_CASES];
    double y[N_CASES];
    double z[N_CASES];
    double ux[N_CASES];
    double uy[N_CASES];
    double uz[N_CASES];
    int inside_batch[N_CASES];
    size_t i;

    fill_cases(x, y, z, ux, uy, uz);

    osh_gemca_runtime_check_body_batch(rt, body_idx, x, y, z, ux, uy, uz, N_CASES, inside_batch);

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

        inside_scalar = reference_in_body(cold_body, &r);
        ASSERT_TRUE(inside_batch[i] == inside_scalar);
    }
}

static void test_body_batch_matches_scalar(void) {
    struct gemca_runtime rt;
    struct gemca_rt_body rt_bodies[3];
    struct gemca_rt_surface rt_surfaces[10];
    struct body cold_bodies[3];
    struct surface cold_surfaces[3][6];
    struct surface *cold_surface_ptrs[3][6];
    double cold_params[3][6][OSH_GEMCA_RT_SURF_NPAR];

    memset(&rt, 0, sizeof(rt));
    memset(rt_bodies, 0, sizeof(rt_bodies));
    memset(rt_surfaces, 0, sizeof(rt_surfaces));
    memset(cold_bodies, 0, sizeof(cold_bodies));
    memset(cold_surfaces, 0, sizeof(cold_surfaces));
    memset(cold_surface_ptrs, 0, sizeof(cold_surface_ptrs));
    memset(cold_params, 0, sizeof(cold_params));

    rt.bodies = rt_bodies;
    rt.nbodies = 3u;
    rt.surfaces = rt_surfaces;
    rt.nsurfaces = 10u;

    /* Body 0: universe sphere, r^2 = 4 */
    rt_bodies[0].surf_begin = 0u;
    rt_bodies[0].nsurfs = 1;
    rt_bodies[0].coord = OSH_COORD_UNIVERSE;
    rt_surfaces[0].type = OSH_GEMCA_SURF_SPHERE;
    rt_surfaces[0].p[0] = 4.0;

    /* Body 1: BCALIGN axis-aligned box, local cube [-1,1]^3 */
    rt_bodies[1].surf_begin = 1u;
    rt_bodies[1].nsurfs = 6;
    rt_bodies[1].coord = OSH_COORD_BCALIGN;
    rt_bodies[1].t[3] = 1.5;
    rt_bodies[1].t[7] = -0.5;
    rt_bodies[1].t[11] = 2.0;

    rt_surfaces[1].type = OSH_GEMCA_SURF_PLANEX; rt_surfaces[1].p[0] = 1.0;  rt_surfaces[1].p[1] = -1.0;
    rt_surfaces[2].type = OSH_GEMCA_SURF_PLANEX; rt_surfaces[2].p[0] = -1.0; rt_surfaces[2].p[1] = -1.0;
    rt_surfaces[3].type = OSH_GEMCA_SURF_PLANEY; rt_surfaces[3].p[0] = 1.0;  rt_surfaces[3].p[1] = -1.0;
    rt_surfaces[4].type = OSH_GEMCA_SURF_PLANEY; rt_surfaces[4].p[0] = -1.0; rt_surfaces[4].p[1] = -1.0;
    rt_surfaces[5].type = OSH_GEMCA_SURF_PLANEZ; rt_surfaces[5].p[0] = 1.0;  rt_surfaces[5].p[1] = -1.0;
    rt_surfaces[6].type = OSH_GEMCA_SURF_PLANEZ; rt_surfaces[6].p[0] = -1.0; rt_surfaces[6].p[1] = -1.0;

    /* Body 2: BZALIGN capped cylinder */
    rt_bodies[2].surf_begin = 7u;
    rt_bodies[2].nsurfs = 3;
    rt_bodies[2].coord = OSH_COORD_BZALIGN;
    rt_bodies[2].t[0] = 0.0;  rt_bodies[2].t[1] = 1.0;  rt_bodies[2].t[2] = 0.0;  rt_bodies[2].t[3] = 1.0;
    rt_bodies[2].t[4] = -1.0; rt_bodies[2].t[5] = 0.0;  rt_bodies[2].t[6] = 0.0;  rt_bodies[2].t[7] = -2.0;
    rt_bodies[2].t[8] = 0.0;  rt_bodies[2].t[9] = 0.0;  rt_bodies[2].t[10] = 1.0; rt_bodies[2].t[11] = 0.5;
    rt_bodies[2].t[15] = 1.0;

    rt_surfaces[7].type = OSH_GEMCA_SURF_CYLZ;   rt_surfaces[7].p[0] = 1.0;
    rt_surfaces[8].type = OSH_GEMCA_SURF_PLANEZ; rt_surfaces[8].p[0] = 1.0;  rt_surfaces[8].p[1] = -2.0;
    rt_surfaces[9].type = OSH_GEMCA_SURF_PLANEZ; rt_surfaces[9].p[0] = -1.0; rt_surfaces[9].p[1] = 0.0;

    init_cold_body_from_runtime(&cold_bodies[0], &rt_bodies[0], rt_surfaces,
                                cold_surfaces[0], cold_surface_ptrs[0], cold_params[0]);
    init_cold_body_from_runtime(&cold_bodies[1], &rt_bodies[1], rt_surfaces,
                                cold_surfaces[1], cold_surface_ptrs[1], cold_params[1]);
    init_cold_body_from_runtime(&cold_bodies[2], &rt_bodies[2], rt_surfaces,
                                cold_surfaces[2], cold_surface_ptrs[2], cold_params[2]);

    compare_body_batch_with_scalar(&rt, 0u, &cold_bodies[0], fill_sphere_cases);
    compare_body_batch_with_scalar(&rt, 1u, &cold_bodies[1], fill_bcalign_box_cases);
    compare_body_batch_with_scalar(&rt, 2u, &cold_bodies[2], fill_bzalign_capped_cyl_cases);
}

int main(void) {
    test_body_batch_matches_scalar();
    return 0;
}
