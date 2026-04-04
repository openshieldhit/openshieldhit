#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "beam/osh_beam_model.h"
#include "common/osh_rc.h"
#include "common/osh_vect.h"
#include "particle/osh_particle.h"
#include "random/osh_rng.h"

#define ASSERT_TRUE(cond)                                                                                              \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            fprintf(stderr, "ASSERT FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__);                                 \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

static void test_single_spot_gaussian_sampling(void) {
    struct beam_workspace wb = {0};
    struct beam_spot spot = {0};
    struct particle part = {0};
    struct particle *part_out = NULL;
    struct ray_v ray = {0};
    struct osh_rng rng;
    struct osh_rng rng_ref;
    double zdir[3] = {0.0, 0.0, 1.0};
    double e_exp;
    double x_exp;
    double xp_exp;
    double y_exp;
    double yp_exp;
    double vnorm;
    int rc;

    part.pdg = 2212;
    spot.part = &part;
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

    osh_vect_setup_tmatrix_bzalign_affine(spot.p, zdir, spot._tm);

    wb.spots = &spot;
    wb.nspots = 1;
    wb.shared.use_sad = 0;

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

    rc = osh_beam_new_primary(&wb, &rng, &part_out, &ray);

    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(part_out == &part);
    ASSERT_TRUE(ray.system == OSH_COORD_UNIVERSE);

    ASSERT_TRUE(fabs(ray.p[0] - (spot.p[0] + x_exp)) < 1e-12);
    ASSERT_TRUE(fabs(ray.p[1] - (spot.p[1] + y_exp)) < 1e-12);
    ASSERT_TRUE(fabs(ray.p[2] - spot.p[2]) < 1e-12);
    ASSERT_TRUE(fabs(ray.p[3] - e_exp) < 1e-12);

    ASSERT_TRUE(fabs(ray.v[0] - xp_exp) < 1e-12);
    ASSERT_TRUE(fabs(ray.v[1] - yp_exp) < 1e-12);
    ASSERT_TRUE(fabs(ray.v[2] - (1.0 / vnorm)) < 1e-12);
    ASSERT_TRUE(fabs(ray.v[0] * ray.v[0] + ray.v[1] * ray.v[1] + ray.v[2] * ray.v[2] - 1.0) < 1e-12);
}

static void test_single_spot_sad_fanout(void) {
    struct beam_workspace wb = {0};
    struct beam_spot spot = {0};
    struct particle part = {0};
    struct particle *part_out = NULL;
    struct ray_v ray = {0};
    struct osh_rng rng;
    double zdir[3] = {0.0, 0.0, 1.0};
    double vx_exp;
    double vy_exp;
    double vz_exp;
    double x_iso;
    double norm;
    int rc;

    part.pdg = 2212;
    spot.part = &part;
    spot.shape = OSH_BEAM_SHAPE_PENCIL;
    spot.p[0] = 1.0;
    spot.p[1] = -2.0;
    spot.p[2] = -50.0;
    spot.t0 = 80.0;

    osh_vect_setup_tmatrix_bzalign_affine(spot.p, zdir, spot._tm);

    wb.spots = &spot;
    wb.nspots = 1;
    wb.shared.use_sad = 1;
    wb.shared.sad[0] = 100.0;
    wb.shared.sad[1] = 200.0;

    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, 7u, 11u);

    vx_exp = spot.p[0] / (spot.p[2] + wb.shared.sad[0]);
    vy_exp = spot.p[1] / (spot.p[2] + wb.shared.sad[1]);
    vz_exp = 1.0;
    norm = sqrt(vx_exp * vx_exp + vy_exp * vy_exp + vz_exp * vz_exp);
    vx_exp /= norm;
    vy_exp /= norm;
    vz_exp /= norm;
    x_iso = spot.p[0] + (0.0 - spot.p[2]) * (spot.p[0] / (spot.p[2] + wb.shared.sad[0]));

    rc = osh_beam_new_primary(&wb, &rng, &part_out, &ray);

    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(part_out == &part);
    ASSERT_TRUE(ray.system == OSH_COORD_UNIVERSE);
    ASSERT_TRUE(fabs(ray.p[0] - spot.p[0]) < 1e-12);
    ASSERT_TRUE(fabs(ray.p[1] - spot.p[1]) < 1e-12);
    ASSERT_TRUE(fabs(ray.p[2] - spot.p[2]) < 1e-12);
    ASSERT_TRUE(fabs(ray.v[0] - vx_exp) < 1e-12);
    ASSERT_TRUE(fabs(ray.v[1] - vy_exp) < 1e-12);
    ASSERT_TRUE(fabs(ray.v[2] - vz_exp) < 1e-12);
    ASSERT_TRUE(fabs(x_iso - 2.0) < 1e-12);
}

int main(void) {
    test_single_spot_gaussian_sampling();
    test_single_spot_sad_fanout();
    return 0;
}
