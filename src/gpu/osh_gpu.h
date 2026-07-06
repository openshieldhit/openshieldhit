#ifndef OSH_GPU_H
#define OSH_GPU_H

/**
 * @file
 * @brief Experimental CUDA backend: device mirrors and query kernels.
 *
 * @details
 * This module is compiled only when OSH_ENABLE_CUDA=ON and nothing in the
 * CPU path depends on it.  It follows the mirror pattern from issue #231:
 * each host runtime gets a small POD "device view" struct holding device
 * pointers + counts, built once per run and passed to kernels by value.
 *
 * Device kernels contain no physics or geometry logic of their own: they
 * call the same OSH_HD static-inline bodies the CPU path executes (see
 * src/gemca/runtime/osh_gemca_runtime_hd.h).
 */

#include <stddef.h>

#include "openshieldhit/status.h"

#ifdef __cplusplus
extern "C" {
#endif

struct osh_gemca_runtime;
struct gemca_rt_surface;
struct gemca_rt_body;
struct gemca_rt_insn;

/**
 * @brief Device-resident mirror of the GEMCA runtime's flat arrays.
 *
 * @details
 * POD struct of device pointers and counts; safe to pass to a kernel by
 * value.  Pointers are owned by the view and released by
 * osh_gpu_gemca_view_free().
 */
struct osh_gpu_gemca_view {
    struct gemca_rt_surface const *surfaces; /**< Device pointer, nsurfaces entries. */
    struct gemca_rt_body const *bodies;      /**< Device pointer, nbodies entries. */
    struct gemca_rt_insn const *insns_flat;  /**< Device pointer, ninsns_flat entries. */
    int const *insn_begin;                   /**< Device pointer, nzones + 1 entries. */
    size_t nsurfaces;                        /**< Number of surfaces. */
    size_t nbodies;                          /**< Number of bodies. */
    size_t nzones;                           /**< Number of zones. */
    size_t ninsns_flat;                      /**< Total flat instructions. */
};

/**
 * @brief Upload the GEMCA runtime's flat arrays to the device.
 *
 * @details
 * Geometries containing VOX (voxel/CT) bodies are refused with OSH_ENOTSUP:
 * struct gemca_rt_body carries a borrowed host pointer to the HU volume and
 * a device kernel must never chase it.  Mirroring CT grids is a separate,
 * explicitly planned step.
 *
 * @param[in]  rt        Compiled gemca runtime (insns_flat populated).
 * @param[out] view_out  Receives device pointers and counts.
 *
 * @returns OSH_OK, or OSH_E* with osh_gpu_last_error() set.
 */
enum osh_status osh_gpu_gemca_view_build(struct osh_gemca_runtime const *rt, struct osh_gpu_gemca_view *view_out);

/** @brief Release all device allocations held by @p view and zero it. */
void osh_gpu_gemca_view_free(struct osh_gpu_gemca_view *view);

/**
 * @brief Zone-membership lookup for a batch of rays on the device.
 *
 * @details
 * Device twin of osh_gemca_runtime_get_zone_batch(): one thread per query
 * ray, each evaluating the identical OSH_HD membership code over the
 * mirrored flat arrays.  Rays outside all zones receive
 * OSH_GEMCA_ZONE_INDEX_INVALID, exactly like the CPU batch API.
 *
 * @param[in]  view          Device mirror from osh_gpu_gemca_view_build().
 * @param[in]  x,y,z         Host SoA positions, length @p n.
 * @param[in]  ux,uy,uz      Host SoA unit directions, length @p n.
 * @param[in]  n             Number of rays.
 * @param[out] zone_out      Host buffer receiving zone indices, length @p n.
 * @param[out] kernel_ms_out Optional; kernel-only time of the last launch
 *                           (CUDA events), excluding H2D/D2H transfers.
 *
 * @returns OSH_OK, or OSH_E* with osh_gpu_last_error() set.
 */
enum osh_status osh_gpu_zone_batch(struct osh_gpu_gemca_view const *view,
                                   double const *x,
                                   double const *y,
                                   double const *z,
                                   double const *ux,
                                   double const *uy,
                                   double const *uz,
                                   size_t n,
                                   size_t *zone_out,
                                   double *kernel_ms_out);

/** @brief Human-readable message for the most recent osh_gpu_* failure. */
char const *osh_gpu_last_error(void);

#ifdef __cplusplus
}
#endif

#endif /* OSH_GPU_H */
