#ifndef OSH_RNG_PCG32_HD_H
#define OSH_RNG_PCG32_HD_H

/*
 * osh_rng_pcg32_hd.h — device-compilable PCG32 engine body.
 *
 * The body is marked OSH_HD static inline so it compiles both as a host
 * function (via plain C compilation) and as a device function (via nvcc
 * compilation with __host__ __device__).  The original .c file includes
 * this header and re-exports the function with its unchanged public
 * signature.
 */

#include "common/osh_hd.h"
#include "random/osh_rng.h"

OSH_HD static inline uint32_t _osh_rng_pcg32_u32_hd(struct osh_rng *rng) {
    uint64_t oldstate;
    uint32_t xorshifted;
    uint32_t rot;

    oldstate = rng->u.pcg32.state;

    rng->u.pcg32.state = oldstate * 6364136223846793005ULL + (rng->u.pcg32.inc | 1ULL);

    xorshifted = (uint32_t) (((oldstate >> 18u) ^ oldstate) >> 27u);
    rot = (uint32_t) (oldstate >> 59u);

    return (xorshifted >> rot) | (xorshifted << ((32u - rot) & 31u));
}

OSH_HD static inline void _osh_rng_pcg32_init_hd(struct osh_rng *rng, uint64_t seed, uint64_t stream) {
    rng->u.pcg32.state = 0ULL;
    rng->u.pcg32.inc = (stream << 1u) | 1u;

    _osh_rng_pcg32_u32_hd(rng);

    rng->u.pcg32.state += seed;

    _osh_rng_pcg32_u32_hd(rng);
}

#endif /* OSH_RNG_PCG32_HD_H */
