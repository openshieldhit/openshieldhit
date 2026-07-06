#ifndef OSH_PHYSICS_BETHE_HD_H
#define OSH_PHYSICS_BETHE_HD_H

/*
 * osh_physics_bethe_hd.h — device-compilable Hubert effective charge Z_eff.
 *
 * The body is marked OSH_HD static inline so it compiles both as a host
 * function (via plain C compilation) and as a device function (via nvcc
 * with __host__ __device__).  The original .c file includes this header
 * and re-exports the function with its unchanged public signature.
 */

#include "common/osh_hd.h"
#include "physics/atomic/osh_physics_bethe.h"

#include <math.h>

OSH_HD static inline double _osh_physics_bethe_z_eff_hd(
    double t_per_nucleon, double proj_z, double proj_a, double target_z_mean) {
    double dd;
    double u1;
    double u2;
    double u3;
    double u4;
    double gamma_eff;

    if (target_z_mean <= 0.0 || proj_z <= 0.0) {
        return proj_z;
    }

    dd = 1.164 + 0.2319 * exp(-0.004302 * target_z_mean);
    u1 = dd + 1.658 * exp(-0.05170 * proj_z);
    u2 = 8.144 + 0.09876 * log(target_z_mean);
    u3 = 0.3140 + 0.01072 * log(target_z_mean);
    u4 = 0.5218 + 0.02521 * log(target_z_mean);
    gamma_eff = 1.0 - u1 * exp((-u2 * pow(t_per_nucleon, u3)) / pow(proj_z, u4));

    if (gamma_eff < 0.0) {
        gamma_eff = 0.0;
    } else if (gamma_eff > 1.0) {
        gamma_eff = 1.0;
    }

    (void) proj_a;
    return proj_z * gamma_eff;
}

#endif /* OSH_PHYSICS_BETHE_HD_H */
