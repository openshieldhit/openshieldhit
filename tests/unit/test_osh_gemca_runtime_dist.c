/*
 * Per-surface distance tests for the GEMCA runtime hot path (issue #255).
 *
 * These exercise osh_gemca_runtime_surface_distance() directly for every
 * quadric surface type, covering inside / outside / tangent / grazing / miss
 * rays. The degenerate-consistency tests (elliptic cylinder with equal radii
 * must equal a plain cylinder, ellipsoid with equal radii must equal a sphere)
 * are the checks that would have caught the missing factor-of-2 in the linear
 * coefficient (bug-hunt finding G-1): before the fix the two code paths for the
 * same surface disagreed by ~30 %.
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_coord.h"
#include "common/osh_ray.h"
#include "gemca/osh_gemca2.h" /* OSH_GEMCA_INFINITY, OSH_GEMCA_SMALL */
#include "gemca/osh_gemca2_defines.h"
#include "gemca/runtime/osh_gemca_runtime.h"

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

static struct gemca_rt_surface make_surf(int type, double p0, double p1, double p2) {
    struct gemca_rt_surface sf;
    memset(&sf, 0, sizeof(sf));
    sf.type = type;
    sf.p[0] = p0;
    sf.p[1] = p1;
    sf.p[2] = p2;
    return sf;
}

static struct ray make_ray(double px, double py, double pz, double ux, double uy, double uz) {
    struct ray r;
    double len;
    r.p[0] = px;
    r.p[1] = py;
    r.p[2] = pz;
    len = sqrt(ux * ux + uy * uy + uz * uz);
    r.cp[0] = ux / len;
    r.cp[1] = uy / len;
    r.cp[2] = uz / len;
    r.system = OSH_COORD_UNIVERSE;
    return r;
}

static double surf_dist(struct gemca_rt_surface const *sf, struct ray const *r) {
    return osh_gemca_runtime_surface_distance(sf, r);
}

/* Deterministic tiny PRNG, same style as test_osh_gemca_runtime_surface.c. */
static uint64_t rng_state = 0x243f6a8885a308d3ull;

static double next_unit(void) {
    rng_state = rng_state * 6364136223846793005ull + 1442695040888963407ull;
    return (double) (rng_state >> 11) * (1.0 / 9007199254740992.0);
}

static double next_signed(double scale) {
    return (2.0 * next_unit() - 1.0) * scale;
}

/* ---- Exact analytic distances (hand-computed) ---------------------------- */

static void test_sphere_distances(void) {
    struct gemca_rt_surface sf = make_surf(OSH_GEMCA_SURF_SPHERE, 4.0, 0.0, 0.0); /* R = 2 */
    struct ray r;
    double d;

    /* from centre, outward: hits at x = 2 */
    r = make_ray(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    ASSERT_CLOSE(surf_dist(&sf, &r), 2.0, 1e-12);

    /* from outside (x=3), moving inward: hits near face at x = 2 */
    r = make_ray(3.0, 0.0, 0.0, -1.0, 0.0, 0.0);
    ASSERT_CLOSE(surf_dist(&sf, &r), 1.0, 1e-12);

    /* grazing/tangent: parallel to +z at x = 2, touches at (2,0,0) at t = 5 */
    r = make_ray(2.0, 0.0, -5.0, 0.0, 0.0, 1.0);
    ASSERT_CLOSE(surf_dist(&sf, &r), 5.0, 1e-9);

    /* complete miss: at x=3 moving +y, never reaches R=2 */
    r = make_ray(3.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    ASSERT_TRUE(surf_dist(&sf, &r) >= OSH_GEMCA_INFINITY);

    /* outside moving away: no positive crossing (sentinel 0.0) */
    r = make_ray(3.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    d = surf_dist(&sf, &r);
    ASSERT_TRUE(d <= 0.0 || d >= OSH_GEMCA_INFINITY);
}

static void test_cylz_distances(void) {
    struct gemca_rt_surface sf = make_surf(OSH_GEMCA_SURF_CYLZ, 4.0, 0.0, 0.0); /* R = 2 */
    struct ray r;

    /* on-axis at z=5, radial +x: hits at x = 2 regardless of z */
    r = make_ray(0.0, 0.0, 5.0, 1.0, 0.0, 0.0);
    ASSERT_CLOSE(surf_dist(&sf, &r), 2.0, 1e-12);

    /* parallel to the axis never crosses the infinite wall */
    r = make_ray(0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
    ASSERT_TRUE(surf_dist(&sf, &r) >= OSH_GEMCA_INFINITY);
}

static void test_ellipsoid_distances(void) {
    /* x^2/4 + y^2/9 + z^2/16 = 1; semi-axes 2, 3, 4 */
    struct gemca_rt_surface sf = make_surf(OSH_GEMCA_SURF_ELLIPSOID, 4.0, 9.0, 16.0);
    struct ray r;

    r = make_ray(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    ASSERT_CLOSE(surf_dist(&sf, &r), 2.0, 1e-12); /* hits x semi-axis */
    r = make_ray(0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    ASSERT_CLOSE(surf_dist(&sf, &r), 3.0, 1e-12); /* hits y semi-axis */
    r = make_ray(0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
    ASSERT_CLOSE(surf_dist(&sf, &r), 4.0, 1e-12); /* hits z semi-axis */
}

static void test_cone_distances(void) {
    /* 45-degree cone x^2 + y^2 - z^2 = 0 (slope^2 = 1, apex offset p[0] unused) */
    struct gemca_rt_surface sf = make_surf(OSH_GEMCA_SURF_CONE, 0.0, 1.0, 0.0);
    struct ray r;
    double d;

    /* from the axis at z=2, radial +x: hits where x = z = 2 */
    r = make_ray(0.0, 0.0, 2.0, 1.0, 0.0, 0.0);
    ASSERT_CLOSE(surf_dist(&sf, &r), 2.0, 1e-12);

    /*
     * Regression for the wrong sign/coefficient on the cone's quadratic 'a'
     * term (bug-hunt finding G-1, "additional suspicion"). This up-and-out ray
     * from inside the cone never crosses the surface in the forward direction
     * (both algebraic roots are negative). The old code returned a spurious
     * finite positive hit (~2.95); the fixed, membership-consistent quadratic
     * returns no forward crossing.
     */
    r = make_ray(0.0, 0.0, 2.0, 0.6, 0.0, 0.8);
    d = surf_dist(&sf, &r);
    ASSERT_TRUE(d <= 0.0 || d >= OSH_GEMCA_INFINITY); /* no finite positive crossing */
}

static void test_plane_distances(void) {
    /* PLANEZ: p[0]*z + p[1] = 0 -> plane at z = 2 (p = {1, -2}) */
    struct gemca_rt_surface sf = make_surf(OSH_GEMCA_SURF_PLANEZ, 1.0, -2.0, 0.0);
    struct ray r;

    r = make_ray(0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
    ASSERT_CLOSE(surf_dist(&sf, &r), 2.0, 1e-12);

    r = make_ray(0.0, 0.0, 5.0, 0.0, 0.0, -1.0);
    ASSERT_CLOSE(surf_dist(&sf, &r), 3.0, 1e-12);
}

/* ---- Degenerate consistency: the G-1 catchers ---------------------------- */

/*
 * An elliptic cylinder with equal semi-axes is a circular cylinder, and an
 * ellipsoid with equal semi-axes is a sphere. The specialised and general code
 * paths must return the same distance for the same ray. Before the fix the
 * elliptic paths passed the half linear coefficient to a full-coefficient
 * solver, so they disagreed with the circular/spherical paths for any ray not
 * aimed straight through the centre.
 */
static void test_ellz_matches_cylz(void) {
    double const r2 = 4.0; /* R = 2 */
    struct gemca_rt_surface cyl = make_surf(OSH_GEMCA_SURF_CYLZ, r2, 0.0, 0.0);
    struct gemca_rt_surface ell = make_surf(OSH_GEMCA_SURF_ELLZ, r2, r2, 0.0);
    size_t i;

    for (i = 0; i < 200u; ++i) {
        struct ray r = make_ray(
            next_signed(5.0), next_signed(5.0), next_signed(5.0), next_signed(1.0), next_signed(1.0), next_signed(1.0));
        double dc = surf_dist(&cyl, &r);
        double de = surf_dist(&ell, &r);
        if (dc >= OSH_GEMCA_INFINITY || de >= OSH_GEMCA_INFINITY) {
            ASSERT_TRUE(dc >= OSH_GEMCA_INFINITY && de >= OSH_GEMCA_INFINITY);
        } else {
            ASSERT_CLOSE(de, dc, 1e-9 * (1.0 + fabs(dc)));
        }
    }
}

static void test_ellipsoid_matches_sphere(void) {
    double const r2 = 9.0; /* R = 3 */
    struct gemca_rt_surface sph = make_surf(OSH_GEMCA_SURF_SPHERE, r2, 0.0, 0.0);
    struct gemca_rt_surface ell = make_surf(OSH_GEMCA_SURF_ELLIPSOID, r2, r2, r2);
    size_t i;

    for (i = 0; i < 200u; ++i) {
        struct ray r = make_ray(
            next_signed(6.0), next_signed(6.0), next_signed(6.0), next_signed(1.0), next_signed(1.0), next_signed(1.0));
        double ds = surf_dist(&sph, &r);
        double de = surf_dist(&ell, &r);
        if (ds >= OSH_GEMCA_INFINITY || de >= OSH_GEMCA_INFINITY) {
            ASSERT_TRUE(ds >= OSH_GEMCA_INFINITY && de >= OSH_GEMCA_INFINITY);
        } else {
            ASSERT_CLOSE(de, ds, 1e-9 * (1.0 + fabs(ds)));
        }
    }
}

/* ---- Distance lands on the membership surface ---------------------------- */

/*
 * For a random set of rays that report a finite crossing, advancing the ray by
 * the returned distance must land on the surface, i.e. the membership predicate
 * flips (or sits exactly on the boundary) there. This ties the distance solver
 * to the inside/outside test the transport loop actually uses, and covers the
 * cone whose surface equation is subtle.
 */
static void assert_distance_on_surface(struct gemca_rt_surface const *sf) {
    size_t i;
    size_t hits = 0u;

    for (i = 0; i < 400u; ++i) {
        struct ray r = make_ray(
            next_signed(4.0), next_signed(4.0), next_signed(4.0), next_signed(1.0), next_signed(1.0), next_signed(1.0));
        double d = surf_dist(sf, &r);
        int in_before;
        int in_after;
        struct ray before;
        struct ray after;

        if (d <= 0.0 || d >= OSH_GEMCA_INFINITY) {
            continue;
        }
        hits++;

        /* sample just inside and just past the reported crossing */
        before = r;
        after = r;
        before.p[0] = r.p[0] + r.cp[0] * (d - 1e-6);
        before.p[1] = r.p[1] + r.cp[1] * (d - 1e-6);
        before.p[2] = r.p[2] + r.cp[2] * (d - 1e-6);
        after.p[0] = r.p[0] + r.cp[0] * (d + 1e-6);
        after.p[1] = r.p[1] + r.cp[1] * (d + 1e-6);
        after.p[2] = r.p[2] + r.cp[2] * (d + 1e-6);

        osh_gemca_runtime_check_surface_batch(
            sf, &before.p[0], &before.p[1], &before.p[2], &before.cp[0], &before.cp[1], &before.cp[2], 1u, &in_before);
        osh_gemca_runtime_check_surface_batch(
            sf, &after.p[0], &after.p[1], &after.p[2], &after.cp[0], &after.cp[1], &after.cp[2], 1u, &in_after);

        /* membership must differ across the crossing point */
        ASSERT_TRUE(in_before != in_after);
    }

    ASSERT_TRUE(hits > 20u); /* the ray set must actually exercise crossings */
}

static void test_distance_on_membership_surface(void) {
    struct gemca_rt_surface sph = make_surf(OSH_GEMCA_SURF_SPHERE, 4.0, 0.0, 0.0);
    struct gemca_rt_surface ell = make_surf(OSH_GEMCA_SURF_ELLIPSOID, 4.0, 9.0, 16.0);
    struct gemca_rt_surface cyl = make_surf(OSH_GEMCA_SURF_CYLZ, 4.0, 0.0, 0.0);
    struct gemca_rt_surface elz = make_surf(OSH_GEMCA_SURF_ELLZ, 4.0, 9.0, 0.0);
    struct gemca_rt_surface con = make_surf(OSH_GEMCA_SURF_CONE, 0.0, 0.25, 0.0);

    assert_distance_on_surface(&sph);
    assert_distance_on_surface(&ell);
    assert_distance_on_surface(&cyl);
    assert_distance_on_surface(&elz);
    assert_distance_on_surface(&con);
}

/* ---- Near-surface numerical stability (G-2) ------------------------------ */

/*
 * A particle sitting just off a surface (c ~ 0, the post-nudge state) must
 * still get accurate roots: the tiny near root and the correct far root. The
 * stable sign-split root form returns both without the catastrophic
 * cancellation of the naive (-b +/- sqrt)/2a form.
 */
static void test_near_surface_roots_accurate(void) {
    double const r_big = 5.0;
    double const eps = 1e-10;
    struct gemca_rt_surface sf = make_surf(OSH_GEMCA_SURF_SPHERE, r_big * r_big, 0.0, 0.0);
    struct ray r;
    double d;

    /* just outside the +x face, moving inward: smallest positive root ~ eps */
    r = make_ray(r_big + eps, 0.0, 0.0, -1.0, 0.0, 0.0);
    d = surf_dist(&sf, &r);
    ASSERT_TRUE(d > 0.0);
    ASSERT_CLOSE(d, eps, 1e-3 * eps); /* accurate small root, not garbage */

    /* exactly on the surface, moving outward: no forward crossing */
    r = make_ray(r_big, 0.0, 0.0, 1.0, 0.0, 0.0);
    d = surf_dist(&sf, &r);
    ASSERT_TRUE(d <= 0.0 || d >= OSH_GEMCA_INFINITY);

    /* exactly on the surface, moving inward: far crossing at 2R */
    r = make_ray(r_big, 0.0, 0.0, -1.0, 0.0, 0.0);
    d = surf_dist(&sf, &r);
    ASSERT_CLOSE(d, 2.0 * r_big, 1e-9);
}

int main(void) {
    test_sphere_distances();
    test_cylz_distances();
    test_ellipsoid_distances();
    test_cone_distances();
    test_plane_distances();
    test_ellz_matches_cylz();
    test_ellipsoid_matches_sphere();
    test_distance_on_membership_surface();
    test_near_surface_roots_accurate();
    return 0;
}
