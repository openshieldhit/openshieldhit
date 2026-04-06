/*
 * gemca_bench.c - Headless geometry raycasting benchmark.
 *
 * Loads a geo.dat file, fires N random rays uniformly across the geometry
 * bounding box, and measures throughput of osh_gemca_zone_index() and
 * osh_gemca_dist() calls.  Reports wall-clock time and rays/second so that
 * algorithmic improvements (or regressions) are immediately visible.
 *
 * Usage: gemca_bench <geo.dat> [nrays]
 *   nrays defaults to 1 000 000.
 *
 * Build type matters: run with -DCMAKE_BUILD_TYPE=Release for meaningful
 * numbers.  Debug builds will be ~10-50x slower due to missing inlining.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "common/osh_logger.h"
#include "common/osh_rc.h"
#include "common/osh_vect.h"
#include "gemca/osh_gemca2.h"
#include "random/osh_rng.h"
#include "transport/osh_transport.h"

#define DEFAULT_NRAYS 1000000

/* Returns wall-clock seconds with at least millisecond resolution. */
static double _now(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double) ts.tv_sec + (double) ts.tv_nsec * 1e-9;
}

/*
 * Derive a simple axis-aligned bounding box from the geometry by sampling
 * the transformation matrices stored on bodies.  This is a best-effort
 * estimate; a dedicated AABB on the workspace would be cleaner.
 */
static void _estimate_bbox(struct gemca_workspace const *g, double bbox_min[3], double bbox_max[3]) {
    size_t ib;
    int i;
    double x;

    for (i = 0; i < 3; i++) {
        bbox_min[i] = -1.0;
        bbox_max[i] = 1.0;
    }

    for (ib = 0; ib < g->nbodies; ib++) {
        struct body const *b = g->bodies[ib];
        if (b == NULL || b->na == 0) {
            continue;
        }
        /* first parameter is typically a radius or half-extent */
        x = fabs(b->a[0]);
        if (x < 1e-12) {
            continue;
        }
        for (i = 0; i < 3; i++) {
            /* translation component from the 4th column of each row */
            double centre = b->t[i * 4 + 3];
            if (centre - x < bbox_min[i])
                bbox_min[i] = centre - x;
            if (centre + x > bbox_max[i])
                bbox_max[i] = centre + x;
        }
    }

    /* add 20 % margin */
    for (i = 0; i < 3; i++) {
        double half = (bbox_max[i] - bbox_min[i]) * 0.1;
        bbox_min[i] -= half;
        bbox_max[i] += half;
    }
}

static void _random_ray(struct ray *r, struct osh_rng *rng, double const bbox_min[3], double const bbox_max[3]) {
    int i;
    double cp[3];
    double len;

    for (i = 0; i < 3; i++) {
        r->p[i] = bbox_min[i] + osh_rng_double(rng) * (bbox_max[i] - bbox_min[i]);
    }

    /* uniform direction on the unit sphere (Marsaglia 1972) */
    do {
        cp[0] = 2.0 * osh_rng_double(rng) - 1.0;
        cp[1] = 2.0 * osh_rng_double(rng) - 1.0;
        cp[2] = 2.0 * osh_rng_double(rng) - 1.0;
        len = cp[0] * cp[0] + cp[1] * cp[1] + cp[2] * cp[2];
    } while (len > 1.0 || len < 1e-12);

    len = sqrt(len);
    r->cp[0] = cp[0] / len;
    r->cp[1] = cp[1] / len;
    r->cp[2] = cp[2] / len;
    r->system = OSH_COORD_UNIVERSE;
}

int main(int argc, char *argv[]) {

    struct gemca_workspace g = {0};
    struct osh_rng rng;
    struct ray r;

    double bbox_min[3];
    double bbox_max[3];

    long nrays;
    long i;
    double t0;
    double t1;
    double elapsed;

    /* Setup logger, select OSH_LOG_DEBUG for more information. */
    osh_log_init(OSH_LOG_INFO, OSH_LOG_F_NONE);

    /* --- argument handling ------------------------------------------------ */
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <geo.dat> [nrays]\n", argv[0]);
        fprintf(stderr, "  nrays defaults to %d\n", DEFAULT_NRAYS);
        return EXIT_FAILURE;
    }

    nrays = (argc >= 3) ? atol(argv[2]) : DEFAULT_NRAYS;
    if (nrays <= 0) {
        fprintf(stderr, "nrays must be positive\n");
        return EXIT_FAILURE;
    }

    /* --- load geometry ---------------------------------------------------- */
    printf("Loading geometry: %s\n", argv[1]);
    if (osh_gemca_load(argv[1], &g) != OSH_OK) {
        fprintf(stderr, "osh_gemca_load() failed\n");
        return EXIT_FAILURE;
    }
    printf("  %llu bodies, %llu zones\n", (unsigned long long) g.nbodies, (unsigned long long) g.nzones);

    _estimate_bbox(&g, bbox_min, bbox_max);
    printf("  bbox  x[%.3g, %.3g]  y[%.3g, %.3g]  z[%.3g, %.3g]  (cm)\n",
           bbox_min[0],
           bbox_max[0],
           bbox_min[1],
           bbox_max[1],
           bbox_min[2],
           bbox_max[2]);

    osh_rng_init(&rng, OSH_RNG_TYPE_XOSHIRO256SS, 42, 0);

    /* --- zone-lookup benchmark -------------------------------------------- */
    printf("\nBenchmark: zone lookup  (N = %ld)\n", nrays);
    fflush(stdout);

    /* volatile sink prevents the compiler from optimising away the loop */
    size_t volatile sink = 0;

    t0 = _now();
    for (i = 0; i < nrays; i++) {
        _random_ray(&r, &rng, bbox_min, bbox_max);
        sink = osh_gemca_zone_index(g, r);
    }
    t1 = _now();
    elapsed = t1 - t0;

    printf("  %.3f s   %.2f Mrays/s\n", elapsed, (double) nrays / elapsed * 1e-6);
    (void) sink;

    /* --- distance benchmark ----------------------------------------------- */
    printf("\nBenchmark: distance query  (N = %ld)\n", nrays);
    fflush(stdout);

    double volatile dsink = 0.0;
    size_t zone_idx;

    t0 = _now();
    for (i = 0; i < nrays; i++) {
        _random_ray(&r, &rng, bbox_min, bbox_max);
        zone_idx = osh_gemca_zone_index(g, r);
        dsink = osh_gemca_dist(g.zones[zone_idx], &r);
    }
    t1 = _now();
    elapsed = t1 - t0;

    printf("  %.3f s   %.2f Mrays/s\n", elapsed, (double) nrays / elapsed * 1e-6);
    (void) dsink;

    return 0;
}
