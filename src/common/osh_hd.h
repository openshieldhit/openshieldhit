#ifndef OSH_HD_H
#define OSH_HD_H

/*
 * osh_hd.h — host/device annotation shim for the experimental GPU backend.
 *
 * The macro OSH_HD expands to __host__ __device__ when the translation unit is
 * compiled by nvcc, and to nothing under a plain C compiler.  Annotating a
 * function with OSH_HD makes its body available to both the host CPU path and
 * the CUDA device kernel, satisfying the M2 requirement that the CPU library's
 * API, layout, and behaviour remain byte-identical.
 *
 * Usage in headers that are already `static inline`:
 *   OSH_HD static inline double osh_material_runtime_sp_lookup(...);
 *
 * Usage in .c files whose body must move to a new *_hd.h header:
 *   The .c file includes the *_hd.h header, which contains the
 *   OSH_HD static inline body.  The .c file then re-exports the function with
 *   its unchanged public signature (a one-line non-inline wrapper, or simply
 *   relying on the inline definition being externally visible).
 *
 * Standard <math.h> calls are fine: CUDA supplies device overloads for the
 * C99 math functions used by the physics and geometry helpers.
 */

#if defined(__CUDACC__)
#define OSH_HD __host__ __device__
#else
#define OSH_HD
#endif

/*
 * OSH_RESTRICT — `restrict` is a C99 keyword with no C++ equivalent, so any
 * public header spelling it bare cannot be included from CUDA (or other C++)
 * translation units.  C compilations keep full restrict semantics; C++ ones
 * get the compiler-specific spelling (GCC/Clang/NVCC: __restrict__,
 * MSVC: __restrict).
 */
#if defined(__cplusplus)
#if defined(_MSC_VER)
#define OSH_RESTRICT __restrict
#else
#define OSH_RESTRICT __restrict__
#endif
#else
#define OSH_RESTRICT restrict
#endif

#endif /* OSH_HD_H */
