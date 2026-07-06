/*
 * osh_gpu_hd_guard.cu — device-buildability guard for the OSH_HD surface.
 *
 * Instantiates every _hd function family (vector ops, RNG, material
 * lookups, atomic physics, GEMCA membership) inside one device kernel, so
 * `OSH_ENABLE_CUDA=ON` fails to *compile* if any OSH_HD body regresses
 * into calling a host-only export or uses a C-only construct.  The
 * archived first port attempt shipped headers that claimed device
 * compilability without ever meeting nvcc; this TU is the guard against
 * that class (rule R2 in docs/dev/gpu_port_plan.md).
 *
 * The kernel is executed once with throwaway inputs; the run only proves
 * the code loads and executes, the numeric twins live in the dedicated
 * parity/known-answer harnesses.
 */

#include <cuda_runtime.h>
#include <stdio.h>

#include "common/osh_ray_hd.h"
#include "common/osh_vect_hd.h"
#include "gemca/runtime/osh_gemca_runtime_hd.h"
#include "material/runtime/osh_material_runtime.h"
#include "physics/atomic/osh_physics_bethe_hd.h"
#include "physics/atomic/osh_physics_scat_highland_hd.h"
#include "physics/atomic/osh_physics_strag_gauss_hd.h"
#include "random/osh_rng_hd.h"

__global__ static void hd_guard_kernel(double *sink) {
    struct osh_rng rng;
    double v[3] = {0.0, 0.0, 1.0};
    double w[3];
    double u1[3];
    double u2[3];
    double acc = 0.0;

    _osh_rng_seed_history_hd(&rng, OSH_RNG_TYPE_PCG32, 1u, 0u, OSH_RNG_PURPOSE_PHYSICS);
    acc += _osh_rng_double_hd(&rng);
    acc += _osh_rng_gauss_hd(&rng, 0.0, 1.0);

    _osh_vect_orthogonal_basis_norm_hd(v, u1, u2);
    acc += _osh_vect_dot_hd(u1, u2);

    acc += _osh_physics_bethe_z_eff_hd(100.0, 1.0, 1.0, 7.42);
    acc += _osh_physics_highland_theta0_hd(100.0, 938.272, 1.0, 0.1, 1.0, 36.08);
    acc += _osh_physics_highland_s_theta_hd(100.0, 938.272, 1.0, 1.0, 36.08, 0.1);
    _osh_physics_highland_scatter_hd(v, w, 0.01, &rng);
    acc += w[0];
    acc += _osh_physics_strag_sigma_hd(1.0, 0.555, 0.1);

    /* Material lookup on a stack-built one-material, two-point table. */
    {
        struct osh_material_runtime tables;
        float sp_tab[2] = {1.0f, 2.0f};
        float range_tab[2] = {0.5f, 1.0f};

        tables.emin = 1.0;
        tables.emax = 100.0;
        tables.log_emin = 0.0;                 /* log(1.0) */
        tables.inv_dlog = 1.0 / log(100.0);    /* 2 grid points */
        tables.nenergy = 2u;
        tables.nprojectiles = 1u;
        tables.mass_stopping_power = sp_tab;
        tables.range_csda = range_tab;

        acc += osh_material_runtime_sp_lookup(&tables, 0u, 0u, 10.0);
        acc += osh_material_runtime_range_lookup(&tables, 0u, 0u, 10.0);
    }

    /* Zone membership on a stack-built one-zone, one-sphere geometry. */
    {
        struct gemca_rt_surface sf;
        struct gemca_rt_body body;
        struct gemca_rt_insn insn;
        int insn_begin[2] = {0, 1};
        struct ray r;

        sf.type = OSH_GEMCA_SURF_SPHERE;
        sf.p[0] = 4.0; /* R^2 */
        sf.p[1] = 0.0;
        sf.p[2] = 0.0;
        sf.p[3] = 0.0;
        body.coord = OSH_COORD_UNIVERSE;
        body.surf_begin = 0u;
        body.nsurfs = 1;
        body.type = OSH_GEMCA_BODY_SPH;
        insn.op = GEMCA_RT_PUSH_BODY;
        insn.operand = 0;
        r.p[0] = 0.5;
        r.p[1] = 0.0;
        r.p[2] = 0.0;
        r.cp[0] = 0.0;
        r.cp[1] = 0.0;
        r.cp[2] = 1.0;
        r.system = OSH_COORD_UNIVERSE;

        acc += (double) _osh_gemca_rt_find_zone_flat_hd(&sf, &body, 1u, &insn, insn_begin, 1u, &r);
    }

    *sink = acc;
}

int main(void) {
    double *d_sink;
    double h_sink = 0.0;
    cudaError_t cerr;

    cudaMalloc(&d_sink, sizeof(double));
    hd_guard_kernel<<<1, 1>>>(d_sink);
    cerr = cudaDeviceSynchronize();
    if (cerr != cudaSuccess) {
        fprintf(stderr, "hd guard kernel failed: %s\n", cudaGetErrorString(cerr));
        return 1;
    }
    cudaMemcpy(&h_sink, d_sink, sizeof(double), cudaMemcpyDeviceToHost);
    cudaFree(d_sink);

    printf("hd guard: ok (sink=%.6f)\n", h_sink);
    return 0;
}
