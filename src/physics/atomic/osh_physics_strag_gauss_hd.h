#ifndef OSH_PHYSICS_STRAG_GAUSS_HD_H
#define OSH_PHYSICS_STRAG_GAUSS_HD_H

/*
 * osh_physics_strag_gauss_hd.h — device-compilable Gaussian (Bohr) σ.
 *
 * The body is marked OSH_HD static inline so it compiles both as a host
 * function (via plain C compilation) and as a device function (via nvcc
 * with __host__ __device__).  The original .c file includes this header
 * and re-exports the function with its unchanged public signature.
 */

#include "common/osh_hd.h"
#include "physics/atomic/osh_physics_strag_gauss.h"

#include <math.h>

#define _BOHR_SQRT_C_HD 0.396128

OSH_HD static inline double _osh_physics_strag_sigma_hd(
    double z_eff, double z_over_a, double thickness_gcm2) {
    if (z_eff <= 0.0 || z_over_a <= 0.0 || thickness_gcm2 <= 0.0) {
        return 0.0;
    }

    return _BOHR_SQRT_C_HD * z_eff * sqrt(z_over_a * thickness_gcm2);
}

#endif /* OSH_PHYSICS_STRAG_GAUSS_HD_H */
