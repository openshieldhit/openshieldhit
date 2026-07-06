#ifndef OSH_RNG_XOSHIRO256SS_HD_H
#define OSH_RNG_XOSHIRO256SS_HD_H

/*
 * osh_rng_xoshiro256ss_hd.h — device-compilable xoshiro256** engine body.
 *
 * Bodies are marked OSH_HD static inline so they compile both as host
 * functions (via plain C compilation) and as device functions (via nvcc
 * with __host__ __device__).  osh_rng_xoshiro256ss.c includes this header
 * and re-exports each function with its unchanged public signature.
 */

#include "common/osh_hd.h"
#include "random/osh_rng.h"

/* rotate left — named rotl64 (no underscore) to avoid colliding with the
 * MSVC compiler intrinsic _rotl64 (C2169: intrinsic cannot be redefined). */
OSH_HD static inline uint64_t _osh_rng_rotl64_hd(uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

/* splitmix64 — see the constant-by-constant rationale in osh_rng_xoshiro256ss.c. */
OSH_HD static inline uint64_t _osh_rng_splitmix64_next_hd(uint64_t *x) {
    uint64_t z;

    *x += 0x9e3779b97f4a7c15ULL;
    z = *x;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

OSH_HD static inline void _osh_rng_xoshiro256ss_init_hd(struct osh_rng *rng, uint64_t seed, uint64_t stream) {
    uint64_t x;

    x = seed ^ (stream * 0x9e3779b97f4a7c15ULL);

    rng->u.xoshiro256ss.s[0] = _osh_rng_splitmix64_next_hd(&x);
    rng->u.xoshiro256ss.s[1] = _osh_rng_splitmix64_next_hd(&x);
    rng->u.xoshiro256ss.s[2] = _osh_rng_splitmix64_next_hd(&x);
    rng->u.xoshiro256ss.s[3] = _osh_rng_splitmix64_next_hd(&x);
}

OSH_HD static inline uint64_t _osh_rng_xoshiro256ss_u64_hd(struct osh_rng *rng) {
    uint64_t *s;
    uint64_t result;
    uint64_t t;

    s = rng->u.xoshiro256ss.s;

    result = _osh_rng_rotl64_hd(s[1] * 5ULL, 7) * 9ULL;
    t = s[1] << 17;

    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];

    s[2] ^= t;
    s[3] = _osh_rng_rotl64_hd(s[3], 45);

    return result;
}

#endif /* OSH_RNG_XOSHIRO256SS_HD_H */
