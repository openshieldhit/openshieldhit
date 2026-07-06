/*
 * osh_gpu_zone_parity_main.c — CPU/GPU zone-membership parity + throughput.
 *
 * Loads a geometry, generates deterministic pseudo-random query rays,
 * computes zone indices three ways — CPU scalar (osh_gemca_runtime_get_zone),
 * CPU batch (osh_gemca_runtime_get_zone_batch, AVX2 when available), GPU
 * kernel (osh_gpu_zone_batch) — and requires exact three-way agreement.
 * Because the GPU kernel compiles the identical OSH_HD source lines, the
 * expected mismatch count is zero, not "statistically small".
 *
 * Also reports throughput (rays/s) for all three paths; the GPU number is
 * given twice: end-to-end (alloc + H2D + kernel + D2H) and kernel-only.
 *
 * Usage: osh_gpu_zone_parity <geo.dat> [box_halfwidth_cm] [n_rays]
 */

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "apps/osh/osh_app_osh.h"
#include "common/osh_coord.h"
#include "common/osh_ray.h"
#include "gemca/osh_gemca2.h"
#include "gemca/runtime/osh_gemca_runtime.h"
#include "openshieldhit/geometry.h"

#include "osh_gpu.h"

#define DEFAULT_N_RAYS (1000000u)
#define DEFAULT_HALFWIDTH (30.0)

/* Deterministic LCG (Knuth MMIX) for reproducible query rays; the production
 * RNG is intentionally not used so this harness has no seed-policy coupling. */
static uint64_t g_lcg = 0x123456789abcdef0ULL;

static double lcg_unit(void) {
    g_lcg = (g_lcg * 6364136223846793005ULL) + 1442695040888963407ULL;
    return (double) (g_lcg >> 11) * (1.0 / 9007199254740992.0);
}

static double wall_s(void) {
    struct timespec ts;

    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double) ts.tv_sec + ((double) ts.tv_nsec * 1e-9);
}

int main(int argc, char **argv) {
    struct osh_geometry_workspace *geom = NULL;
    struct osh_gemca_runtime rt = {0};
    struct osh_gpu_gemca_view view;
    double *x, *y, *z, *ux, *uy, *uz;
    size_t *zone_scalar, *zone_batch, *zone_gpu;
    double halfwidth = DEFAULT_HALFWIDTH;
    size_t n = DEFAULT_N_RAYS;
    size_t i;
    size_t bad_batch = 0;
    size_t bad_gpu = 0;
    double t0, t_scalar, t_batch, t_gpu_e2e;
    double kernel_ms = 0.0;
    struct ray r;
    enum osh_status rc;
    int rep;

    if (argc < 2) {
        fprintf(stderr, "usage: %s <geo.dat> [box_halfwidth_cm] [n_rays]\n", argv[0]);
        return 2;
    }
    if (argc > 2) {
        halfwidth = atof(argv[2]);
    }
    if (argc > 3) {
        n = (size_t) atoll(argv[3]);
    }

    rc = osh_geometry_setup_from_path(argv[1], NULL, &geom);
    if (rc != OSH_OK) {
        fprintf(stderr, "failed to load geometry %s (rc=%d)\n", argv[1], rc);
        return 1;
    }
    rc = osh_gemca_compile(geom->prepared, OSH_HU_TABLE_NONE, 0u, NULL, &rt);
    if (rc != OSH_OK) {
        fprintf(stderr, "failed to compile geometry (rc=%d)\n", rc);
        return 1;
    }
    printf("geometry: %s\n", argv[1]);
    printf("  zones %zu, bodies %zu, surfaces %zu, insns_flat %zu, cpu batch dispatch: %s\n", rt.nzones, rt.nbodies,
           rt.nsurfaces, rt.ninsns_flat, osh_gemca_runtime_zone_batch_dispatch_name(&rt));
    printf("  rays %zu in box +-%.1f cm\n", n, halfwidth);

    x = malloc(n * sizeof(double));
    y = malloc(n * sizeof(double));
    z = malloc(n * sizeof(double));
    ux = malloc(n * sizeof(double));
    uy = malloc(n * sizeof(double));
    uz = malloc(n * sizeof(double));
    zone_scalar = malloc(n * sizeof(size_t));
    zone_batch = malloc(n * sizeof(size_t));
    zone_gpu = malloc(n * sizeof(size_t));
    if (!x || !y || !z || !ux || !uy || !uz || !zone_scalar || !zone_batch || !zone_gpu) {
        fprintf(stderr, "out of memory\n");
        return 1;
    }

    for (i = 0; i < n; i++) {
        double cth, sth, phi;

        x[i] = ((2.0 * lcg_unit()) - 1.0) * halfwidth;
        y[i] = ((2.0 * lcg_unit()) - 1.0) * halfwidth;
        z[i] = ((2.0 * lcg_unit()) - 1.0) * halfwidth;
        cth = (2.0 * lcg_unit()) - 1.0;
        sth = sqrt(1.0 - (cth * cth));
        phi = 2.0 * 3.14159265358979323846 * lcg_unit();
        ux[i] = sth * cos(phi);
        uy[i] = sth * sin(phi);
        uz[i] = cth;
    }

    /* CPU scalar reference */
    t0 = wall_s();
    for (i = 0; i < n; i++) {
        r.p[0] = x[i];
        r.p[1] = y[i];
        r.p[2] = z[i];
        r.cp[0] = ux[i];
        r.cp[1] = uy[i];
        r.cp[2] = uz[i];
        r.system = OSH_COORD_UNIVERSE;
        zone_scalar[i] = osh_gemca_runtime_get_zone(&rt, &r);
    }
    t_scalar = wall_s() - t0;

    /* CPU batch (AVX2 when the node supports it) */
    t0 = wall_s();
    osh_gemca_runtime_get_zone_batch(&rt, x, y, z, ux, uy, uz, n, zone_batch);
    t_batch = wall_s() - t0;

    /* GPU: mirror once, one warmup batch, then the timed batch */
    rc = osh_gpu_gemca_view_build(&rt, &view);
    if (rc != OSH_OK) {
        fprintf(stderr, "gpu mirror failed: %s\n", osh_gpu_last_error());
        return 1;
    }
    rc = osh_gpu_zone_batch(&view, x, y, z, ux, uy, uz, n, zone_gpu, &kernel_ms);
    if (rc != OSH_OK) {
        fprintf(stderr, "gpu warmup batch failed: %s\n", osh_gpu_last_error());
        return 1;
    }
    t_gpu_e2e = 0.0;
    for (rep = 0; rep < 3; rep++) {
        double t_rep = wall_s();
        double ms_rep = 0.0;

        rc = osh_gpu_zone_batch(&view, x, y, z, ux, uy, uz, n, zone_gpu, &ms_rep);
        if (rc != OSH_OK) {
            fprintf(stderr, "gpu batch failed: %s\n", osh_gpu_last_error());
            return 1;
        }
        t_rep = wall_s() - t_rep;
        if (rep == 0 || t_rep < t_gpu_e2e) {
            t_gpu_e2e = t_rep;
            kernel_ms = ms_rep;
        }
    }

    /* Three-way exact comparison */
    for (i = 0; i < n; i++) {
        if (zone_batch[i] != zone_scalar[i]) {
            if (bad_batch < 5u) {
                fprintf(stderr, "batch mismatch at ray %zu (%.6f %.6f %.6f): scalar %zd batch %zd\n", i, x[i], y[i],
                        z[i], (ptrdiff_t) zone_scalar[i], (ptrdiff_t) zone_batch[i]);
            }
            bad_batch++;
        }
        if (zone_gpu[i] != zone_scalar[i]) {
            if (bad_gpu < 5u) {
                fprintf(stderr, "gpu mismatch at ray %zu (%.6f %.6f %.6f): scalar %zd gpu %zd\n", i, x[i], y[i], z[i],
                        (ptrdiff_t) zone_scalar[i], (ptrdiff_t) zone_gpu[i]);
            }
            bad_gpu++;
        }
    }

    printf("parity: cpu-batch mismatches %zu / %zu, gpu mismatches %zu / %zu\n", bad_batch, n, bad_gpu, n);
    printf("throughput [Mrays/s]: cpu scalar %.2f | cpu batch %.2f | gpu end-to-end %.2f | gpu kernel-only %.2f\n",
           (double) n / t_scalar / 1e6, (double) n / t_batch / 1e6, (double) n / t_gpu_e2e / 1e6,
           (double) n / (kernel_ms * 1e-3) / 1e6);
    printf("times: cpu scalar %.3f s | cpu batch %.3f s | gpu e2e %.3f s | gpu kernel %.3f ms\n", t_scalar, t_batch,
           t_gpu_e2e, kernel_ms);

    osh_gpu_gemca_view_free(&view);
    osh_gemca_runtime_free(&rt);

    if (bad_gpu != 0u || bad_batch != 0u) {
        printf("FAIL\n");
        return 1;
    }
    printf("PASS\n");
    return 0;
}
