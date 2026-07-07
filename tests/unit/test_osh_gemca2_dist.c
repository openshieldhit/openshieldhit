/*
 * Reference-path distance tests for GEMCA (issue #255).
 *
 * osh_gemca_get_distance() is the pointer-linked AST reference implementation
 * used by the standalone tools (SDL viewer, bench). It carries the same quadric
 * distance math as the transport runtime, and had the same bugs:
 *   - the elliptic-cylinder / cone / ellipsoid linear coefficient was missing a
 *     factor of 2 (finding G-1), and
 *   - the cone quadratic did not match the cone membership surface at all.
 *
 * These tests build minimal single-body zones and drive the public entry point.
 * The degenerate-consistency checks (elliptic cylinder with equal radii equals
 * a plain cylinder; ellipsoid with equal radii equals a sphere) are the ones
 * that catch the missing factor of 2.
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_coord.h"
#include "common/osh_ray.h"
#include "gemca/osh_gemca2.h"
#include "gemca/osh_gemca2_defines.h"

#define ASSERT_TRUE(cond)                                                                                              \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            fprintf(stderr, "ASSERT FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__);                                 \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

#define ASSERT_CLOSE(got, want, tol)                                                                                   \
    do {                                                                                                               \
        double _g = (got);                                                                                             \
        double _w = (want);                                                                                            \
        if (fabs(_g - _w) > (tol)) {                                                                                   \
            fprintf(stderr,                                                                                            \
                    "ASSERT_CLOSE FAILED: got %.17g, want %.17g, tol %.3g (%s:%d)\n",                                  \
                    _g,                                                                                                \
                    _w,                                                                                                \
                    (tol),                                                                                             \
                    __FILE__,                                                                                          \
                    __LINE__);                                                                                         \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

enum { ZONE_DISTANCE_MAX_PARAMS = 4 }; /* PLANE has the most surface parameters (A,B,C,D) */

/* Build a single-body zone from one surface, in OSH_COORD_UNIVERSE (identity). */
static double zone_distance(int surf_type, double const *params, int np, struct ray const *r) {
    struct surface sf;
    struct surface *surfs[1];
    struct body b;
    struct zone z;
    struct ray rc;
    double pbuf[ZONE_DISTANCE_MAX_PARAMS]; /* mutable copy so sf.p never aliases const input */
    int i;

    memset(&sf, 0, sizeof(sf));
    memset(&b, 0, sizeof(b));
    memset(&z, 0, sizeof(z));

    ASSERT_TRUE(np >= 0 && np <= ZONE_DISTANCE_MAX_PARAMS);
    for (i = 0; i < np; ++i) {
        pbuf[i] = params[i];
    }

    sf.type = surf_type;
    sf.np = np;
    sf.p = pbuf;

    surfs[0] = &sf;
    b.surfs = surfs;
    b.nsurfs = 1;
    b.coord = (char) OSH_COORD_UNIVERSE;
    b.type = OSH_GEMCA_BODY_NONE;

    z.node.type = _OSH_GEMCA_CGNODE_BODY;
    z.node.body = &b;

    rc = *r;
    rc.system = OSH_COORD_UNIVERSE;
    return osh_gemca_get_distance(&z, &rc);
}

static struct ray make_ray(double px, double py, double pz, double ux, double uy, double uz) {
    struct ray r;
    r.p[0] = px;
    r.p[1] = py;
    r.p[2] = pz;
    r.cp[0] = ux;
    r.cp[1] = uy;
    r.cp[2] = uz;
    r.system = OSH_COORD_UNIVERSE;
    return r;
}

static uint64_t rng_state = 0xcafef00dd15ea5e5ull;

static double next_unit(void) {
    rng_state = rng_state * 6364136223846793005ull + 1442695040888963407ull;
    return (double) (rng_state >> 11) * (1.0 / 9007199254740992.0);
}

static double next_signed(double scale) {
    return (2.0 * next_unit() - 1.0) * scale;
}

static void test_ast_sphere_exit_distance(void) {
    double sphere[1] = {4.0}; /* R = 2 */
    struct ray r = make_ray(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    ASSERT_CLOSE(zone_distance(OSH_GEMCA_SURF_SPHERE, sphere, 1, &r), 2.0, 1e-9);
}

static void test_ast_ellz_matches_cylz(void) {
    double const R2 = 4.0; /* R = 2 */
    double cyl[1] = {R2};
    double ell[2] = {R2, R2};
    size_t i;

    for (i = 0; i < 100u; ++i) {
        /* interior start (radius < R), random direction with an xy component */
        struct ray r = make_ray(
            next_signed(1.2), next_signed(1.2), next_signed(3.0), next_signed(1.0), next_signed(1.0), next_signed(0.5));
        double dc = zone_distance(OSH_GEMCA_SURF_CYLZ, cyl, 1, &r);
        double de = zone_distance(OSH_GEMCA_SURF_ELLZ, ell, 2, &r);
        if (dc >= OSH_GEMCA_INFINITY || de >= OSH_GEMCA_INFINITY) {
            ASSERT_TRUE(dc >= OSH_GEMCA_INFINITY && de >= OSH_GEMCA_INFINITY);
        } else {
            ASSERT_CLOSE(de, dc, 1e-6 * (1.0 + fabs(dc)));
        }
    }
}

static void test_ast_ellipsoid_matches_sphere(void) {
    double const R2 = 9.0; /* R = 3 */
    double sph[1] = {R2};
    double ell[3] = {R2, R2, R2};
    size_t i;

    for (i = 0; i < 100u; ++i) {
        /* interior start (|p| < R), random direction */
        struct ray r = make_ray(
            next_signed(1.5), next_signed(1.5), next_signed(1.5), next_signed(1.0), next_signed(1.0), next_signed(1.0));
        double ds = zone_distance(OSH_GEMCA_SURF_SPHERE, sph, 1, &r);
        double de = zone_distance(OSH_GEMCA_SURF_ELLIPSOID, ell, 3, &r);
        if (ds >= OSH_GEMCA_INFINITY || de >= OSH_GEMCA_INFINITY) {
            ASSERT_TRUE(ds >= OSH_GEMCA_INFINITY && de >= OSH_GEMCA_INFINITY);
        } else {
            ASSERT_CLOSE(de, ds, 1e-6 * (1.0 + fabs(ds)));
        }
    }
}

/*
 * Regression for the cone quadratic. A ray launched from inside a 45-degree
 * cone (x^2 + y^2 - z^2 = 0) up and along the axis never crosses the surface,
 * so the particle never leaves the zone: the distance is infinite. The old
 * cone quadratic reported a spurious finite crossing here.
 */
static void test_ast_cone_ray_never_exits(void) {
    double cone[2] = {0.0, 1.0}; /* apex offset unused, slope^2 = 1 */
    struct ray r = make_ray(0.0, 0.0, 2.0, 0.6, 0.0, 0.8);
    ASSERT_TRUE(zone_distance(OSH_GEMCA_SURF_CONE, cone, 2, &r) >= OSH_GEMCA_INFINITY);
}

/*
 * Exercise the solver's degenerate (a ~ 0) linear branch and its no-real-roots
 * branch, which the "generic quadratic" rays above never reach:
 *   - a purely axial ray in an infinite cylinder gives a = 0 and h = 0;
 *   - a ray parallel to a cone generator gives a = 0 with h != 0 (one root);
 *   - a ray that entirely misses the cone gives a real quadratic with no roots.
 */
static void test_ast_solver_degenerate_branches(void) {
    double cyl[1] = {4.0};       /* R = 2 */
    double cone[2] = {0.0, 1.0}; /* slope^2 = 1 */
    struct ray r;

    /* axial ray: a = 0 and h = 0 -> never crosses the infinite wall */
    r = make_ray(0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
    ASSERT_TRUE(zone_distance(OSH_GEMCA_SURF_CYLZ, cyl, 1, &r) >= OSH_GEMCA_INFINITY);

    /* generator-parallel ray with a single forward crossing at distance sqrt(2) */
    r = make_ray(0.0, 0.0, -2.0, 1.0, 0.0, 1.0);
    ASSERT_CLOSE(zone_distance(OSH_GEMCA_SURF_CONE, cone, 2, &r), sqrt(2.0), 1e-9);

    /* generator-parallel ray whose only crossing is behind the start: never exits */
    r = make_ray(0.0, 0.0, 2.0, 1.0, 0.0, 1.0);
    ASSERT_TRUE(zone_distance(OSH_GEMCA_SURF_CONE, cone, 2, &r) >= OSH_GEMCA_INFINITY);

    /* start outside the cone on a ray that misses it entirely (no real roots);
       the walk never enters the zone, so the path length is zero */
    r = make_ray(5.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    ASSERT_CLOSE(zone_distance(OSH_GEMCA_SURF_CONE, cone, 2, &r), 0.0, 1e-12);
}

int main(void) {
    test_ast_sphere_exit_distance();
    test_ast_ellz_matches_cylz();
    test_ast_ellipsoid_matches_sphere();
    test_ast_cone_ray_never_exits();
    test_ast_solver_degenerate_branches();
    return 0;
}
