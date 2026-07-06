/*
 * osh_gpu_gemca_mirror.cu — GEMCA runtime device mirror and zone kernel.
 *
 * Uploads the runtime's flat arrays (surfaces, bodies, insns_flat,
 * insn_begin) to device memory and runs the zone-membership query as one
 * thread per ray.  The kernel body is _osh_gemca_rt_find_zone_flat_hd()
 * from osh_gemca_runtime_hd.h — the same source lines the CPU executes —
 * so CPU/GPU parity is exact by construction (the arithmetic is identical
 * IEEE-754 double operations in identical order; no fast-math, no FMA
 * contraction differences in the membership predicates' comparisons).
 */

#include <cuda_runtime.h>
#include <stdio.h>
#include <string.h>

#include "common/osh_coord.h"
#include "common/osh_ray.h"
#include "gemca/osh_gemca2_defines.h"
#include "gemca/runtime/osh_gemca_runtime.h"
#include "gemca/runtime/osh_gemca_runtime_hd.h"

#include "osh_gpu.h"

/* ---- Error reporting ----------------------------------------------------- */

static char g_last_err[512] = {0};

extern "C" char const *osh_gpu_last_error(void) {
    return g_last_err[0] ? g_last_err : "no error";
}

static void set_last_error(char const *what, cudaError_t cerr) {
    snprintf(g_last_err, sizeof(g_last_err), "%s: %s", what, cudaGetErrorString(cerr));
}

/* ---- Mirror build/free --------------------------------------------------- */

/* cudaMalloc + cudaMemcpy of one flat array; returns NULL and sets the error
 * message on failure. */
static void *upload_array(void const *host, size_t bytes, char const *what) {
    void *dev = NULL;
    cudaError_t cerr;

    cerr = cudaMalloc(&dev, bytes);
    if (cerr != cudaSuccess) {
        set_last_error(what, cerr);
        return NULL;
    }
    cerr = cudaMemcpy(dev, host, bytes, cudaMemcpyHostToDevice);
    if (cerr != cudaSuccess) {
        set_last_error(what, cerr);
        cudaFree(dev);
        return NULL;
    }
    return dev;
}

extern "C" enum osh_status osh_gpu_gemca_view_build(struct osh_gemca_runtime const *rt,
                                                    struct osh_gpu_gemca_view *view_out) {
    struct osh_gpu_gemca_view view;
    size_t i;
    int dev_count = 0;

    if (!rt || !view_out || !rt->insns_flat || !rt->insn_begin) {
        return OSH_EINVAL;
    }

    if (cudaGetDeviceCount(&dev_count) != cudaSuccess || dev_count <= 0) {
        snprintf(g_last_err, sizeof(g_last_err), "%s", "no CUDA-capable device found");
        return OSH_ESTATE;
    }

    /* struct gemca_rt_body carries a borrowed host pointer (hu) for VOX
     * bodies; a device kernel must never chase it.  Refuse rather than
     * silently mirror a landmine. */
    for (i = 0; i < rt->nbodies; i++) {
        if (rt->bodies[i].type == OSH_GEMCA_BODY_VOX) {
            snprintf(g_last_err, sizeof(g_last_err), "%s", "VOX bodies are not mirrored yet (host hu pointer)");
            return OSH_ENOTSUP;
        }
    }

    memset(&view, 0, sizeof(view));
    view.nsurfaces = rt->nsurfaces;
    view.nbodies = rt->nbodies;
    view.nzones = rt->nzones;
    view.ninsns_flat = rt->ninsns_flat;

    view.surfaces =
        (struct gemca_rt_surface const *) upload_array(rt->surfaces, rt->nsurfaces * sizeof(*rt->surfaces), "surfaces");
    view.bodies = (struct gemca_rt_body const *) upload_array(rt->bodies, rt->nbodies * sizeof(*rt->bodies), "bodies");
    view.insns_flat =
        (struct gemca_rt_insn const *) upload_array(rt->insns_flat, rt->ninsns_flat * sizeof(*rt->insns_flat), "insns_flat");
    view.insn_begin = (int const *) upload_array(rt->insn_begin, (rt->nzones + 1u) * sizeof(*rt->insn_begin), "insn_begin");

    if (!view.surfaces || !view.bodies || !view.insns_flat || !view.insn_begin) {
        osh_gpu_gemca_view_free(&view);
        return OSH_ENOMEM;
    }

    *view_out = view;
    return OSH_OK;
}

extern "C" void osh_gpu_gemca_view_free(struct osh_gpu_gemca_view *view) {
    if (!view) {
        return;
    }
    cudaFree((void *) view->surfaces);
    cudaFree((void *) view->bodies);
    cudaFree((void *) view->insns_flat);
    cudaFree((void *) view->insn_begin);
    memset(view, 0, sizeof(*view));
}

/* ---- Zone-membership kernel ---------------------------------------------- */

__global__ static void zone_batch_kernel(struct osh_gpu_gemca_view view,
                                         double const *x,
                                         double const *y,
                                         double const *z,
                                         double const *ux,
                                         double const *uy,
                                         double const *uz,
                                         size_t n,
                                         size_t *zone_out) {
    size_t i = (size_t) blockIdx.x * blockDim.x + threadIdx.x;
    struct ray r;

    if (i >= n) {
        return;
    }

    r.p[0] = x[i];
    r.p[1] = y[i];
    r.p[2] = z[i];
    r.cp[0] = ux[i];
    r.cp[1] = uy[i];
    r.cp[2] = uz[i];
    r.system = OSH_COORD_UNIVERSE;

    zone_out[i] = _osh_gemca_rt_find_zone_flat_hd(view.surfaces, view.bodies, view.nbodies, view.insns_flat,
                                                  view.insn_begin, view.nzones, &r);
}

extern "C" enum osh_status osh_gpu_zone_batch(struct osh_gpu_gemca_view const *view,
                                              double const *x,
                                              double const *y,
                                              double const *z,
                                              double const *ux,
                                              double const *uy,
                                              double const *uz,
                                              size_t n,
                                              size_t *zone_out,
                                              double *kernel_ms_out) {
    double *d_soa = NULL;
    size_t *d_zone = NULL;
    double const *h_soa[6];
    size_t const nbytes = n * sizeof(double);
    cudaError_t cerr;
    cudaEvent_t ev_a;
    cudaEvent_t ev_b;
    float ms = 0.0f;
    int block = 256;
    int grid;
    size_t k;

    if (!view || !zone_out || n == 0u) {
        return OSH_EINVAL;
    }

    /* One allocation for all six SoA input arrays, laid out back to back. */
    cerr = cudaMalloc((void **) &d_soa, 6u * nbytes);
    if (cerr != cudaSuccess) {
        set_last_error("soa alloc", cerr);
        return OSH_ENOMEM;
    }
    cerr = cudaMalloc((void **) &d_zone, n * sizeof(size_t));
    if (cerr != cudaSuccess) {
        set_last_error("zone alloc", cerr);
        cudaFree(d_soa);
        return OSH_ENOMEM;
    }

    h_soa[0] = x;
    h_soa[1] = y;
    h_soa[2] = z;
    h_soa[3] = ux;
    h_soa[4] = uy;
    h_soa[5] = uz;
    for (k = 0; k < 6u; k++) {
        cerr = cudaMemcpy(d_soa + (k * n), h_soa[k], nbytes, cudaMemcpyHostToDevice);
        if (cerr != cudaSuccess) {
            set_last_error("soa upload", cerr);
            cudaFree(d_soa);
            cudaFree(d_zone);
            return OSH_ESTATE;
        }
    }

    cudaEventCreate(&ev_a);
    cudaEventCreate(&ev_b);

    grid = (int) ((n + (size_t) block - 1u) / (size_t) block);
    cudaEventRecord(ev_a);
    zone_batch_kernel<<<grid, block>>>(*view, d_soa + (0u * n), d_soa + (1u * n), d_soa + (2u * n), d_soa + (3u * n),
                                       d_soa + (4u * n), d_soa + (5u * n), n, d_zone);
    cudaEventRecord(ev_b);

    cerr = cudaDeviceSynchronize();
    if (cerr != cudaSuccess) {
        set_last_error("kernel", cerr);
        cudaEventDestroy(ev_a);
        cudaEventDestroy(ev_b);
        cudaFree(d_soa);
        cudaFree(d_zone);
        return OSH_ESTATE;
    }
    cudaEventElapsedTime(&ms, ev_a, ev_b);
    if (kernel_ms_out) {
        *kernel_ms_out = (double) ms;
    }

    cerr = cudaMemcpy(zone_out, d_zone, n * sizeof(size_t), cudaMemcpyDeviceToHost);

    cudaEventDestroy(ev_a);
    cudaEventDestroy(ev_b);
    cudaFree(d_soa);
    cudaFree(d_zone);

    if (cerr != cudaSuccess) {
        set_last_error("zone download", cerr);
        return OSH_ESTATE;
    }
    return OSH_OK;
}
