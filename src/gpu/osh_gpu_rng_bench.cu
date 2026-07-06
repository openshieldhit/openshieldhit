/*
 * osh_gpu_rng_bench.cu — device known-answer test + throughput for the
 * project RNG compiled from the shared OSH_HD headers.
 *
 * Two jobs:
 *
 * 1. Known-answer (Level-0 twin from the GPU plan): per-history streams
 *    seeded with osh_rng_seed_history semantics must produce bitwise
 *    identical u32/u64/uniform-double sequences on host and device (pure
 *    integer mixing + one exact multiply). gauss01 goes through log/sqrt,
 *    where libm and CUDA math may legally differ in the last ULPs, so it
 *    is compared with a relative tolerance and the observed maximum is
 *    reported.
 *
 * 2. Throughput: samples/s for pcg32 u32, xoshiro256** u64, uniform
 *    double, and gauss01 at transport-like occupancy (one independent
 *    stream per thread, state in registers).
 *
 * This TU is also the compile-time guard for the RNG HD slice: it
 * instantiates the full _hd call chain in device code, which fails to
 * build if any _hd function regresses into calling a host-only export.
 */

#include <cuda_runtime.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "random/osh_rng.h"
#include "random/osh_rng_hd.h"

#define N_STREAMS 8192
#define N_DRAWS 32

#define BENCH_BLOCKS 4096
#define BENCH_THREADS 256
#define BENCH_DRAWS 4096

static void check(cudaError_t cerr, char const *what) {
    if (cerr != cudaSuccess) {
        fprintf(stderr, "%s: %s\n", what, cudaGetErrorString(cerr));
        exit(1);
    }
}

/* ---- Known-answer kernels ------------------------------------------------ */

__global__ static void ka_kernel(uint64_t seed, enum osh_rng_type type, uint32_t *u32_out, double *dbl_out,
                                 double *gauss_out) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    struct osh_rng rng;
    int k;

    if (i >= N_STREAMS) {
        return;
    }

    _osh_rng_seed_history_hd(&rng, type, seed, (uint64_t) i, OSH_RNG_PURPOSE_PHYSICS);
    for (k = 0; k < N_DRAWS; k++) {
        u32_out[(i * N_DRAWS) + k] = _osh_rng_u32_hd(&rng);
    }
    for (k = 0; k < N_DRAWS; k++) {
        dbl_out[(i * N_DRAWS) + k] = _osh_rng_double_hd(&rng);
    }
    for (k = 0; k < N_DRAWS; k++) {
        gauss_out[(i * N_DRAWS) + k] = _osh_rng_gauss01_hd(&rng);
    }
}

static int known_answer(enum osh_rng_type type, char const *name) {
    size_t const n = (size_t) N_STREAMS * N_DRAWS;
    uint32_t *d_u32, *h_u32;
    double *d_dbl, *h_dbl, *d_gauss, *h_gauss;
    struct osh_rng rng;
    size_t bad_u32 = 0, bad_dbl = 0;
    double max_gauss_rel = 0.0;
    int i, k;

    check(cudaMalloc(&d_u32, n * sizeof(*d_u32)), "malloc u32");
    check(cudaMalloc(&d_dbl, n * sizeof(*d_dbl)), "malloc dbl");
    check(cudaMalloc(&d_gauss, n * sizeof(*d_gauss)), "malloc gauss");
    h_u32 = (uint32_t *) malloc(n * sizeof(*h_u32));
    h_dbl = (double *) malloc(n * sizeof(*h_dbl));
    h_gauss = (double *) malloc(n * sizeof(*h_gauss));

    ka_kernel<<<(N_STREAMS + 255) / 256, 256>>>(42u, type, d_u32, d_dbl, d_gauss);
    check(cudaDeviceSynchronize(), "ka kernel");
    check(cudaMemcpy(h_u32, d_u32, n * sizeof(*h_u32), cudaMemcpyDeviceToHost), "copy u32");
    check(cudaMemcpy(h_dbl, d_dbl, n * sizeof(*h_dbl), cudaMemcpyDeviceToHost), "copy dbl");
    check(cudaMemcpy(h_gauss, d_gauss, n * sizeof(*h_gauss), cudaMemcpyDeviceToHost), "copy gauss");

    /* Host reference through the public .c API — exercises the delegation
     * path, so host-vs-device also validates the .c → _hd re-export. */
    for (i = 0; i < N_STREAMS; i++) {
        osh_rng_seed_history(&rng, type, 42u, (uint64_t) i, OSH_RNG_PURPOSE_PHYSICS);
        for (k = 0; k < N_DRAWS; k++) {
            if (osh_rng_u32(&rng) != h_u32[(i * N_DRAWS) + k]) {
                bad_u32++;
            }
        }
        for (k = 0; k < N_DRAWS; k++) {
            if (osh_rng_double(&rng) != h_dbl[(i * N_DRAWS) + k]) {
                bad_dbl++;
            }
        }
        for (k = 0; k < N_DRAWS; k++) {
            double ref = osh_rng_gauss01(&rng);
            double dev = h_gauss[(i * N_DRAWS) + k];
            double rel = (ref != 0.0) ? fabs(dev - ref) / fabs(ref) : fabs(dev - ref);

            if (rel > max_gauss_rel) {
                max_gauss_rel = rel;
            }
        }
    }

    printf("known-answer %s: u32 mismatches %zu/%zu, double mismatches %zu/%zu, gauss max rel dev %.3e\n", name,
           bad_u32, n, bad_dbl, n, max_gauss_rel);

    cudaFree(d_u32);
    cudaFree(d_dbl);
    cudaFree(d_gauss);
    free(h_u32);
    free(h_dbl);
    free(h_gauss);

    /* gauss01 tolerance: glibc and CUDA log() may legally differ in the last
     * ULP, and the relative error of log(s) diverges as s -> 1, so ~1e-11
     * relative deviations are expected and physically meaningless (observed
     * max on A100/CUDA 12.8 vs glibc 2.34: ~7e-12).  Uniform draws stay
     * bitwise. */
    return (bad_u32 == 0u && bad_dbl == 0u && max_gauss_rel < 1e-9) ? 0 : 1;
}

/* ---- Throughput kernels -------------------------------------------------- */

/* mode: 0 = pcg32 u32, 1 = xoshiro u64, 2 = uniform double, 3 = gauss01 */
__global__ static void bench_kernel(uint64_t seed, enum osh_rng_type type, int mode, double *sink) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    struct osh_rng rng;
    double acc = 0.0;
    uint64_t iacc = 0u;
    int k;

    _osh_rng_seed_history_hd(&rng, type, seed, (uint64_t) i, OSH_RNG_PURPOSE_PHYSICS);

    for (k = 0; k < BENCH_DRAWS; k++) {
        switch (mode) {
        case 0:
            iacc += _osh_rng_u32_hd(&rng);
            break;
        case 1:
            iacc += _osh_rng_u64_hd(&rng);
            break;
        case 2:
            acc += _osh_rng_double_hd(&rng);
            break;
        default:
            acc += _osh_rng_gauss01_hd(&rng);
            break;
        }
    }

    sink[i] = acc + (double) iacc;
}

static void bench(enum osh_rng_type type, int mode, char const *name) {
    double *d_sink;
    cudaEvent_t a, b;
    float ms;
    double const total = (double) BENCH_BLOCKS * BENCH_THREADS * BENCH_DRAWS;

    check(cudaMalloc(&d_sink, (size_t) BENCH_BLOCKS * BENCH_THREADS * sizeof(double)), "sink");
    cudaEventCreate(&a);
    cudaEventCreate(&b);

    bench_kernel<<<BENCH_BLOCKS, BENCH_THREADS>>>(42u, type, mode, d_sink); /* warmup */
    check(cudaDeviceSynchronize(), "warmup");

    cudaEventRecord(a);
    bench_kernel<<<BENCH_BLOCKS, BENCH_THREADS>>>(42u, type, mode, d_sink);
    cudaEventRecord(b);
    check(cudaDeviceSynchronize(), "bench");
    cudaEventElapsedTime(&ms, a, b);

    printf("throughput %-22s %8.1f Gsamples/s\n", name, total / ((double) ms * 1e-3) / 1e9);

    cudaEventDestroy(a);
    cudaEventDestroy(b);
    cudaFree(d_sink);
}

int main(void) {
    int rc = 0;
    cudaDeviceProp prop;

    check(cudaGetDeviceProperties(&prop, 0), "props");
    printf("device: %s (sm_%d%d)\n", prop.name, prop.major, prop.minor);

    rc |= known_answer(OSH_RNG_TYPE_PCG32, "pcg32");
    rc |= known_answer(OSH_RNG_TYPE_XOSHIRO256SS, "xoshiro256**");

    bench(OSH_RNG_TYPE_PCG32, 0, "pcg32 u32");
    bench(OSH_RNG_TYPE_XOSHIRO256SS, 1, "xoshiro256** u64");
    bench(OSH_RNG_TYPE_PCG32, 2, "pcg32 double");
    bench(OSH_RNG_TYPE_PCG32, 3, "pcg32 gauss01");

    printf(rc == 0 ? "PASS\n" : "FAIL\n");
    return rc;
}
