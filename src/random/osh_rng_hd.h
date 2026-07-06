#ifndef OSH_RNG_HD_H
#define OSH_RNG_HD_H

/*
 * osh_rng_hd.h — device-compilable RNG init, seeding, and distributions.
 *
 * Bodies are marked OSH_HD static inline so they compile both as host
 * functions (via plain C compilation) and as device functions (via nvcc
 * with __host__ __device__).  The original .c file includes this header
 * and re-exports each with its unchanged public signature.
 *
 * Every internal call inside this header MUST target another _hd twin,
 * never a host .c export: nvcc rejects a __host__ __device__ function
 * calling a __host__-only one in device compilation, and that failure
 * only surfaces when a kernel actually instantiates the call chain.  The
 * device build of src/gpu/osh_gpu_rng_bench.cu is the guard that keeps
 * this property true.
 */

#include "common/osh_hd.h"
#include "random/osh_rng.h"
#include "random/osh_rng_pcg32_hd.h"
#include "random/osh_rng_xoshiro256ss_hd.h"

#include <math.h>

/* Body of osh_rng_init() (see osh_rng.c). */
OSH_HD static inline void _osh_rng_init_hd(struct osh_rng *rng, enum osh_rng_type type, uint64_t seed, uint64_t stream) {
    rng->type = type;
    rng->gauss_has_spare = 0;

    switch (type) {
    case OSH_RNG_TYPE_PCG32:
        _osh_rng_pcg32_init_hd(rng, seed, stream);
        break;

    case OSH_RNG_TYPE_XOSHIRO256SS:
        _osh_rng_xoshiro256ss_init_hd(rng, seed, stream);
        break;

    default:
        _osh_rng_pcg32_init_hd(rng, seed, stream);
        rng->type = OSH_RNG_TYPE_PCG32;
        break;
    }
}

OSH_HD static inline uint64_t _rng_mix_stream_hd(uint64_t seed, uint64_t hist_index, uint64_t purpose) {
    uint64_t x = seed ^ (hist_index * 0x9E3779B97F4A7C15ULL) ^ (purpose * 0xD1B54A32D192ED03ULL);

    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

OSH_HD static inline void _osh_rng_seed_history_hd(
    struct osh_rng *rng, enum osh_rng_type type, uint64_t seed, uint64_t hist_index, enum osh_rng_purpose purpose) {
    uint64_t const stream = _rng_mix_stream_hd(seed, hist_index, (uint64_t) purpose);

    _osh_rng_init_hd(rng, type, seed, stream);
}

OSH_HD static inline uint32_t _osh_rng_u32_hd(struct osh_rng *rng) {
    switch (rng->type) {
    case OSH_RNG_TYPE_PCG32:
        return _osh_rng_pcg32_u32_hd(rng);

    case OSH_RNG_TYPE_XOSHIRO256SS:
        return (uint32_t) (_osh_rng_xoshiro256ss_u64_hd(rng) >> 32);

    default:
        return _osh_rng_pcg32_u32_hd(rng);
    }
}

OSH_HD static inline uint64_t _osh_rng_u64_hd(struct osh_rng *rng) {
    uint64_t hi;
    uint64_t lo;

    switch (rng->type) {
    case OSH_RNG_TYPE_PCG32:
        hi = (uint64_t) _osh_rng_pcg32_u32_hd(rng);
        lo = (uint64_t) _osh_rng_pcg32_u32_hd(rng);
        return (hi << 32) | lo;

    case OSH_RNG_TYPE_XOSHIRO256SS:
        return _osh_rng_xoshiro256ss_u64_hd(rng);

    default:
        return ((uint64_t) _osh_rng_u32_hd(rng) << 32) | (uint64_t) _osh_rng_u32_hd(rng);
    }
}

OSH_HD static inline float _osh_rng_float_hd(struct osh_rng *rng) {
    uint32_t r;
    uint32_t mant;

    r = _osh_rng_u32_hd(rng);
    mant = r >> 8;

    return (float) mant * (1.0f / 16777216.0f);
}

OSH_HD static inline double _osh_rng_double_hd(struct osh_rng *rng) {
    uint64_t r;
    uint64_t mant;

    r = _osh_rng_u64_hd(rng);
    mant = r >> 11;

    return (double) mant * (1.0 / 9007199254740992.0);
}

OSH_HD static inline double _osh_rng_gauss01_hd(struct osh_rng *rng) {
    double u;
    double v;
    double s;
    double m;

    if (rng->gauss_has_spare) {
        rng->gauss_has_spare = 0;
        return rng->gauss_spare;
    }

    do {
        u = 2.0 * _osh_rng_double_hd(rng) - 1.0;
        v = 2.0 * _osh_rng_double_hd(rng) - 1.0;
        s = u * u + v * v;
    } while (s <= 0.0 || s >= 1.0);

    m = sqrt(-2.0 * log(s) / s);

    rng->gauss_spare = v * m;
    rng->gauss_has_spare = 1;

    return u * m;
}

OSH_HD static inline double _osh_rng_gauss_hd(struct osh_rng *rng, double mu, double sigma) {
    return mu + (sigma * _osh_rng_gauss01_hd(rng));
}

#endif /* OSH_RNG_HD_H */
