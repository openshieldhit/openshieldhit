/*
 * osh_gpu_atomics_bench.cu — FP64 atomicAdd contention on scoring-like
 * access patterns.
 *
 * The GPU scoring decision (issue #231 D4) is atomicAdd-first with
 * privatization as the measured fallback.  This bench quantifies the two
 * regimes that decide it: a tiny 1D grid where many threads hit the same
 * bins (Bragg-curve depth-dose) and large CT-like grids where collisions
 * are rare.  Bin indices come from the project PCG32 so the pattern is
 * random rather than compiler-foldable.
 */

#include <cuda_runtime.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "random/osh_rng.h"
#include "random/osh_rng_hd.h"

#define BLOCKS 2048
#define THREADS 256
#define OPS_PER_THREAD 1024

static void check(cudaError_t cerr, char const *what) {
    if (cerr != cudaSuccess) {
        fprintf(stderr, "%s: %s\n", what, cudaGetErrorString(cerr));
        exit(1);
    }
}

/* nbins == 0 selects the no-atomic reference: plain store to the thread's
 * own slot, same RNG work, giving the atomic-free ceiling. */
__global__ static void atomic_kernel(double *bins, unsigned int nbins) {
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    struct osh_rng rng;
    unsigned int k;
    unsigned int bin;
    double own = 0.0;

    _osh_rng_seed_history_hd(&rng, OSH_RNG_TYPE_PCG32, 42u, (uint64_t) i, OSH_RNG_PURPOSE_PHYSICS);

    for (k = 0; k < OPS_PER_THREAD; k++) {
        if (nbins == 0u) {
            own += (double) _osh_rng_u32_hd(&rng);
        } else {
            bin = _osh_rng_u32_hd(&rng) % nbins;
            atomicAdd(&bins[bin], 1.0);
        }
    }

    if (nbins == 0u) {
        bins[i] = own;
    }
}

static void run(double *d_bins, unsigned int nbins, char const *name) {
    cudaEvent_t a, b;
    float ms;
    double const total = (double) BLOCKS * THREADS * OPS_PER_THREAD;

    cudaEventCreate(&a);
    cudaEventCreate(&b);

    atomic_kernel<<<BLOCKS, THREADS>>>(d_bins, nbins); /* warmup */
    check(cudaDeviceSynchronize(), name);

    cudaEventRecord(a);
    atomic_kernel<<<BLOCKS, THREADS>>>(d_bins, nbins);
    cudaEventRecord(b);
    check(cudaDeviceSynchronize(), name);
    cudaEventElapsedTime(&ms, a, b);

    printf("%-28s %8.2f Gatomics/s\n", name, total / ((double) ms * 1e-3) / 1e9);

    cudaEventDestroy(a);
    cudaEventDestroy(b);
}

int main(void) {
    double *d_bins;
    size_t const max_bins = 1u << 20;
    cudaDeviceProp prop;

    check(cudaGetDeviceProperties(&prop, 0), "props");
    printf("device: %s (sm_%d%d), %d ops/thread, %d threads\n", prop.name, prop.major, prop.minor, OPS_PER_THREAD,
           BLOCKS * THREADS);

    check(cudaMalloc(&d_bins, max_bins * sizeof(double)), "bins");
    check(cudaMemset(d_bins, 0, max_bins * sizeof(double)), "memset");

    run(d_bins, 0u, "no-atomic reference");
    run(d_bins, 1u, "1 bin (all colliding)");
    run(d_bins, 64u, "64 bins");
    run(d_bins, 256u, "256 bins (1D Bragg-like)");
    run(d_bins, 4096u, "4096 bins");
    run(d_bins, 65536u, "64k bins");
    run(d_bins, 1u << 20, "1M bins (CT-like)");

    cudaFree(d_bins);
    return 0;
}
