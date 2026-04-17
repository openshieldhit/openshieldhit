#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "beam/osh_beam.h"
#include "beam/osh_beam_model.h"
#include "beam/osh_beam_prepared.h"
#include "common/osh_const.h"
#include "common/osh_rc.h"
#include "particle/osh_particle_pdg.h"
#include "random/osh_rng.h"

#define ASSERT_TRUE(cond)                                                                                              \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            fprintf(stderr, "ASSERT FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__);                                 \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

static void cleanup_manual_wb(struct beam_workspace *wb) {
    if (!wb || !wb->prepared) {
        return;
    }
    free(wb->prepared->cum_wt);
    free(wb->prepared->tm);
    free(wb->prepared);
    wb->prepared = NULL;
}

static void test_single_spot_gaussian_sampling(void) {
    struct beam_workspace wb = {0};
    struct beam_spot spot = {0};
    struct ray_v ray = {0};
    struct osh_rng rng;
    struct osh_rng rng_ref;
    double e_exp;
    double x_exp;
    double xp_exp;
    double y_exp;
    double yp_exp;
    double vnorm;
    enum osh_status rc;

    wb.spots = &spot;
    wb.nspots = 1;
    wb.primary.pdg = OSH_PART_PDG_PROTON;
    wb.primary.z = 1u;
    wb.primary.a = 1u;
    wb.has_primary = 1;
    spot.shape = OSH_BEAM_SHAPE_GAUSSIAN;
    spot.p[0] = 4.0;
    spot.p[1] = 5.0;
    spot.p[2] = 6.0;
    spot.size[0] = 0.2;
    spot.size[1] = 0.3;
    spot.div[0] = 0.01;
    spot.div[1] = 0.02;
    spot.cor[0] = 0.25;
    spot.cor[1] = -0.5;
    spot.t0 = 100.0;
    spot.tsigma = 1.5;

    wb.shared.use_sad = 0;
    wb.shared.theta = 0.0;
    wb.shared.phi = 0.0;
    rc = osh_beam_workspace_prepare(&wb);
    ASSERT_TRUE(rc == OSH_OK);

    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, 42u, 54u);
    osh_rng_init(&rng_ref, OSH_RNG_TYPE_PCG32, 42u, 54u);

    e_exp = osh_rng_gauss(&rng_ref, spot.t0, spot.tsigma);
    if (e_exp < 0.0) {
        e_exp = 0.0;
    }

    x_exp = spot.size[0] * osh_rng_gauss01(&rng_ref);
    xp_exp =
        spot.div[0]
        * (spot.cor[0] * (x_exp / spot.size[0]) + sqrt(1.0 - spot.cor[0] * spot.cor[0]) * osh_rng_gauss01(&rng_ref));
    y_exp = spot.size[1] * osh_rng_gauss01(&rng_ref);
    yp_exp =
        spot.div[1]
        * (spot.cor[1] * (y_exp / spot.size[1]) + sqrt(1.0 - spot.cor[1] * spot.cor[1]) * osh_rng_gauss01(&rng_ref));

    vnorm = sqrt(xp_exp * xp_exp + yp_exp * yp_exp + 1.0);
    xp_exp /= vnorm;
    yp_exp /= vnorm;

    rc = osh_beam_new_primary(&wb, &rng, &ray);

    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(ray.system == OSH_COORD_UNIVERSE);

    ASSERT_TRUE(fabs(ray.p[0] - (spot.p[0] + x_exp)) < 1e-12);
    ASSERT_TRUE(fabs(ray.p[1] - (spot.p[1] + y_exp)) < 1e-12);
    ASSERT_TRUE(fabs(ray.p[2] - spot.p[2]) < 1e-12);
    ASSERT_TRUE(fabs(ray.p[3] - e_exp) < 1e-12);

    ASSERT_TRUE(fabs(ray.v[0] - xp_exp) < 1e-12);
    ASSERT_TRUE(fabs(ray.v[1] - yp_exp) < 1e-12);
    ASSERT_TRUE(fabs(ray.v[2] - (1.0 / vnorm)) < 1e-12);
    ASSERT_TRUE(fabs(ray.v[0] * ray.v[0] + ray.v[1] * ray.v[1] + ray.v[2] * ray.v[2] - 1.0) < 1e-12);
    cleanup_manual_wb(&wb);
}

static void test_single_spot_sad_fanout(void) {
    struct beam_workspace wb = {0};
    struct beam_spot spot = {0};
    struct ray_v ray = {0};
    struct osh_rng rng;
    double vx_exp;
    double vy_exp;
    double vz_exp;
    double x_iso;
    double norm;
    enum osh_status rc;

    wb.spots = &spot;
    wb.nspots = 1;
    wb.primary.pdg = OSH_PART_PDG_PROTON;
    wb.primary.z = 1u;
    wb.primary.a = 1u;
    wb.has_primary = 1;
    spot.shape = OSH_BEAM_SHAPE_PENCIL;
    spot.p[0] = 1.0;
    spot.p[1] = -2.0;
    spot.p[2] = -50.0;
    spot.t0 = 80.0;

    wb.shared.use_sad = 1;
    wb.shared.sad[0] = 100.0;
    wb.shared.sad[1] = 200.0;
    wb.shared.theta = 0.0;
    wb.shared.phi = 0.0;
    rc = osh_beam_workspace_prepare(&wb);
    ASSERT_TRUE(rc == OSH_OK);

    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, 7u, 11u);

    vx_exp = spot.p[0] / (spot.p[2] + wb.shared.sad[0]);
    vy_exp = spot.p[1] / (spot.p[2] + wb.shared.sad[1]);
    vz_exp = 1.0;
    norm = sqrt(vx_exp * vx_exp + vy_exp * vy_exp + vz_exp * vz_exp);
    vx_exp /= norm;
    vy_exp /= norm;
    vz_exp /= norm;
    x_iso = spot.p[0] + (0.0 - spot.p[2]) * (spot.p[0] / (spot.p[2] + wb.shared.sad[0]));

    rc = osh_beam_new_primary(&wb, &rng, &ray);

    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(ray.system == OSH_COORD_UNIVERSE);
    ASSERT_TRUE(fabs(ray.p[0] - spot.p[0]) < 1e-12);
    ASSERT_TRUE(fabs(ray.p[1] - spot.p[1]) < 1e-12);
    ASSERT_TRUE(fabs(ray.p[2] - spot.p[2]) < 1e-12);
    ASSERT_TRUE(fabs(ray.v[0] - vx_exp) < 1e-12);
    ASSERT_TRUE(fabs(ray.v[1] - vy_exp) < 1e-12);
    ASSERT_TRUE(fabs(ray.v[2] - vz_exp) < 1e-12);
    ASSERT_TRUE(fabs(x_iso - 2.0) < 1e-12);
    cleanup_manual_wb(&wb);
}

static void test_single_spot_square_sampling(void) {
    struct beam_workspace wb = {0};
    struct beam_spot spot = {0};
    struct ray_v ray = {0};
    struct osh_rng rng;
    struct osh_rng rng_ref;
    double x_exp;
    double y_exp;
    enum osh_status rc;

    wb.spots = &spot;
    wb.nspots = 1;
    wb.primary.pdg = OSH_PART_PDG_PROTON;
    wb.primary.z = 1u;
    wb.primary.a = 1u;
    wb.has_primary = 1;
    spot.shape = OSH_BEAM_SHAPE_SQUARE;
    spot.p[0] = -3.0;
    spot.p[1] = 2.5;
    spot.p[2] = -40.0;
    spot.size[0] = 1.5;
    spot.size[1] = 0.5;
    spot.t0 = 120.0;

    wb.shared.use_sad = 0;
    wb.shared.theta = 0.0;
    wb.shared.phi = 0.0;
    rc = osh_beam_workspace_prepare(&wb);
    ASSERT_TRUE(rc == OSH_OK);

    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, 123u, 456u);
    osh_rng_init(&rng_ref, OSH_RNG_TYPE_PCG32, 123u, 456u);

    x_exp = (2.0 * osh_rng_double(&rng_ref) - 1.0) * spot.size[0];
    y_exp = (2.0 * osh_rng_double(&rng_ref) - 1.0) * spot.size[1];

    rc = osh_beam_new_primary(&wb, &rng, &ray);

    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(ray.system == OSH_COORD_UNIVERSE);
    ASSERT_TRUE(fabs(ray.p[0] - (spot.p[0] + x_exp)) < 1e-12);
    ASSERT_TRUE(fabs(ray.p[1] - (spot.p[1] + y_exp)) < 1e-12);
    ASSERT_TRUE(fabs(ray.p[2] - spot.p[2]) < 1e-12);
    ASSERT_TRUE(fabs(ray.p[3] - spot.t0) < 1e-12);
    ASSERT_TRUE(fabs(ray.v[0]) < 1e-12);
    ASSERT_TRUE(fabs(ray.v[1]) < 1e-12);
    ASSERT_TRUE(fabs(ray.v[2] - 1.0) < 1e-12);
    ASSERT_TRUE(fabs(x_exp) <= spot.size[0]);
    ASSERT_TRUE(fabs(y_exp) <= spot.size[1]);
    cleanup_manual_wb(&wb);
}

static void test_single_spot_circular_sampling(void) {
    struct beam_workspace wb = {0};
    struct beam_spot spot = {0};
    struct ray_v ray = {0};
    struct osh_rng rng;
    struct osh_rng rng_ref;
    double phi_exp;
    double r2_exp;
    double r_exp;
    double x_exp;
    double y_exp;
    double r_world;
    enum osh_status rc;

    wb.spots = &spot;
    wb.nspots = 1;
    wb.primary.pdg = OSH_PART_PDG_PROTON;
    wb.primary.z = 1u;
    wb.primary.a = 1u;
    wb.has_primary = 1;
    spot.shape = OSH_BEAM_SHAPE_CIRCULAR;
    spot.p[0] = 1.0;
    spot.p[1] = -1.5;
    spot.p[2] = -25.0;
    spot.size[0] = 2.0;
    spot.size[1] = 0.5;
    spot.t0 = 90.0;

    wb.shared.use_sad = 0;
    wb.shared.theta = 0.0;
    wb.shared.phi = 0.0;
    rc = osh_beam_workspace_prepare(&wb);
    ASSERT_TRUE(rc == OSH_OK);

    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, 999u, 1001u);
    osh_rng_init(&rng_ref, OSH_RNG_TYPE_PCG32, 999u, 1001u);

    phi_exp = 2.0 * OSH_M_PI * osh_rng_double(&rng_ref);
    r2_exp = 0.5 * 0.5 + (2.0 * 2.0 - 0.5 * 0.5) * osh_rng_double(&rng_ref);
    r_exp = sqrt(r2_exp);
    x_exp = r_exp * cos(phi_exp);
    y_exp = r_exp * sin(phi_exp);

    rc = osh_beam_new_primary(&wb, &rng, &ray);

    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(ray.system == OSH_COORD_UNIVERSE);
    ASSERT_TRUE(fabs(ray.p[0] - (spot.p[0] + x_exp)) < 1e-12);
    ASSERT_TRUE(fabs(ray.p[1] - (spot.p[1] + y_exp)) < 1e-12);
    ASSERT_TRUE(fabs(ray.p[2] - spot.p[2]) < 1e-12);
    ASSERT_TRUE(fabs(ray.p[3] - spot.t0) < 1e-12);
    ASSERT_TRUE(fabs(ray.v[0]) < 1e-12);
    ASSERT_TRUE(fabs(ray.v[1]) < 1e-12);
    ASSERT_TRUE(fabs(ray.v[2] - 1.0) < 1e-12);

    r_world = sqrt((ray.p[0] - spot.p[0]) * (ray.p[0] - spot.p[0]) + (ray.p[1] - spot.p[1]) * (ray.p[1] - spot.p[1]));
    ASSERT_TRUE(r_world >= 0.5 - 1e-12);
    ASSERT_TRUE(r_world <= 2.0 + 1e-12);
    cleanup_manual_wb(&wb);
}

static void test_weighted_spot_selection(void) {
    struct beam_workspace wb = {0};
    struct beam_spot spots[3] = {0};
    struct ray_v ray = {0};
    struct osh_rng rng;
    struct osh_rng rng_ref;
    double cum_wt[3] = {1.0, 3.0, 10.0};
    double w;
    size_t idx_exp;
    int i;
    enum osh_status rc;

    wb.spots = spots;
    wb.nspots = 3;
    wb.primary.pdg = OSH_PART_PDG_PROTON;
    wb.primary.z = 1u;
    wb.primary.a = 1u;
    wb.has_primary = 1;
    for (i = 0; i < 3; i++) {
        spots[i].shape = OSH_BEAM_SHAPE_PENCIL;
        spots[i].p[0] = 10.0 * (double) i;
        spots[i].p[1] = -5.0 * (double) i;
        spots[i].p[2] = -20.0 - (double) i;
        spots[i].t0 = 70.0 + 10.0 * (double) i;
        spots[i].wt = (i == 0) ? 1.0 : (i == 1) ? 2.0 : 7.0;
    }

    wb.shared.use_sad = 0;
    wb.shared.theta = 0.0;
    wb.shared.phi = 0.0;
    rc = osh_beam_workspace_prepare(&wb);
    ASSERT_TRUE(rc == OSH_OK);

    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, 321u, 654u);
    osh_rng_init(&rng_ref, OSH_RNG_TYPE_PCG32, 321u, 654u);

    for (i = 0; i < 8; i++) {
        w = osh_rng_double(&rng_ref) * wb.prepared->wt_sum;
        if (w <= cum_wt[0]) {
            idx_exp = 0;
        } else if (w <= cum_wt[1]) {
            idx_exp = 1;
        } else {
            idx_exp = 2;
        }

        rc = osh_beam_new_primary(&wb, &rng, &ray);

        ASSERT_TRUE(rc == OSH_OK);
        ASSERT_TRUE(fabs(ray.p[0] - spots[idx_exp].p[0]) < 1e-12);
        ASSERT_TRUE(fabs(ray.p[1] - spots[idx_exp].p[1]) < 1e-12);
        ASSERT_TRUE(fabs(ray.p[2] - spots[idx_exp].p[2]) < 1e-12);
        ASSERT_TRUE(fabs(ray.p[3] - spots[idx_exp].t0) < 1e-12);
    }
    cleanup_manual_wb(&wb);
}

int main(void) {
    test_single_spot_gaussian_sampling();
    test_single_spot_sad_fanout();
    test_single_spot_square_sampling();
    test_single_spot_circular_sampling();
    test_weighted_spot_selection();
    return 0;
}
