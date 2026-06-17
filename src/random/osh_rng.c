#include "random/osh_rng.h"

#include <math.h>

void osh_rng_init(struct osh_rng *rng, enum osh_rng_type type, uint64_t seed, uint64_t stream) {
    rng->type = type;
    rng->gauss_has_spare = 0;

    switch (type) {
    case OSH_RNG_TYPE_PCG32:
        osh_rng_pcg32_init(rng, seed, stream);
        break;

    case OSH_RNG_TYPE_XOSHIRO256SS:
        osh_rng_xoshiro256ss_init(rng, seed, stream);
        break;

    default:
        osh_rng_pcg32_init(rng, seed, stream);
        rng->type = OSH_RNG_TYPE_PCG32;
        break;
    }
}

/*
 * SplitMix64 finaliser — mixes three 64-bit inputs into a well-distributed
 * stream id.  The golden-ratio and prime odd multipliers de-correlate adjacent
 * history indices and purposes before the avalanche, so consecutive histories
 * map to far-apart, statistically independent streams.  This is the same
 * mixing used to seed xoshiro from (seed, stream); reusing it here keeps the
 * lane-separation properties identical across both engines.
 */
static uint64_t rng_mix_stream(uint64_t seed, uint64_t hist_index, uint64_t purpose) {
    uint64_t x = seed ^ (hist_index * 0x9E3779B97F4A7C15ULL) ^ (purpose * 0xD1B54A32D192ED03ULL);

    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

void osh_rng_seed_history(
    struct osh_rng *rng, enum osh_rng_type type, uint64_t seed, uint64_t hist_index, enum osh_rng_purpose purpose) {
    uint64_t const stream = rng_mix_stream(seed, hist_index, (uint64_t) purpose);

    osh_rng_init(rng, type, seed, stream);
}

void osh_rng_split(struct osh_rng *child, struct osh_rng *parent) {
    /* Two draws from the parent supply independent seed and stream entropy.
     * Deterministic along the parent's lineage, hence reproducible regardless
     * of how the wavefront schedules sibling histories. */
    uint64_t const child_seed = osh_rng_u64(parent);
    uint64_t const child_stream = osh_rng_u64(parent);

    osh_rng_init(child, parent->type, child_seed, child_stream);
}

uint32_t osh_rng_u32(struct osh_rng *rng) {
    switch (rng->type) {
    case OSH_RNG_TYPE_PCG32:
        return osh_rng_pcg32_u32(rng);

    case OSH_RNG_TYPE_XOSHIRO256SS:
        return (uint32_t) (osh_rng_xoshiro256ss_u64(rng) >> 32);

    default:
        return osh_rng_pcg32_u32(rng);
    }
}

uint64_t osh_rng_u64(struct osh_rng *rng) {
    uint64_t hi;
    uint64_t lo;

    switch (rng->type) {
    case OSH_RNG_TYPE_PCG32:
        hi = (uint64_t) osh_rng_pcg32_u32(rng);
        lo = (uint64_t) osh_rng_pcg32_u32(rng);
        return (hi << 32) | lo;

    case OSH_RNG_TYPE_XOSHIRO256SS:
        return osh_rng_xoshiro256ss_u64(rng);

    default:
        return ((uint64_t) osh_rng_u32(rng) << 32) | (uint64_t) osh_rng_u32(rng);
    }
}

float osh_rng_float(struct osh_rng *rng) {
    uint32_t r;
    uint32_t mant;

    /* Use 24 bits for float mantissa */
    r = osh_rng_u32(rng);
    mant = r >> 8; /* top 24 bits */

    return (float) mant * (1.0f / 16777216.0f); /* 2^24 */
}

double osh_rng_double(struct osh_rng *rng) {
    uint64_t r;
    uint64_t mant;

    /* Use 53 bits for double mantissa */
    r = osh_rng_u64(rng);
    mant = r >> 11; /* top 53 bits */

    return (double) mant * (1.0 / 9007199254740992.0); /* 2^53 */
}

double osh_rng_gauss01(struct osh_rng *rng) {
    double u;
    double v;
    double s;
    double m;

    if (rng->gauss_has_spare) {
        rng->gauss_has_spare = 0;
        return rng->gauss_spare;
    }

    /* Rejection sample in the unit disk */
    do {
        u = 2.0 * osh_rng_double(rng) - 1.0;
        v = 2.0 * osh_rng_double(rng) - 1.0;
        s = u * u + v * v;
    } while (s <= 0.0 || s >= 1.0);

    m = sqrt(-2.0 * log(s) / s);

    rng->gauss_spare = v * m;
    rng->gauss_has_spare = 1;

    return u * m;
}

double osh_rng_gauss(struct osh_rng *rng, double mu, double sigma) {
    return mu + (sigma * osh_rng_gauss01(rng));
}

/* vectors */

void osh_rng_double_vec(struct osh_rng *rng, double *restrict x, int n) {
    int i = 0;

    while (i + 3 < n) {
        x[i + 0] = osh_rng_double(rng);
        x[i + 1] = osh_rng_double(rng);
        x[i + 2] = osh_rng_double(rng);
        x[i + 3] = osh_rng_double(rng);
        i += 4;
    }
    while (i < n) {
        x[i] = osh_rng_double(rng);
        i++;
    }
}

void osh_rng_float_vec(struct osh_rng *rng, float *restrict x, int n) {
    int i = 0;

    while (i + 3 < n) {
        x[i + 0] = osh_rng_float(rng);
        x[i + 1] = osh_rng_float(rng);
        x[i + 2] = osh_rng_float(rng);
        x[i + 3] = osh_rng_float(rng);
        i += 4;
    }
    while (i < n) {
        x[i] = osh_rng_float(rng);
        i++;
    }
}

void osh_rng_gauss01_vec(struct osh_rng *rng, double *restrict x, int n) {
    int i = 0;

    while (i + 3 < n) {
        x[i + 0] = osh_rng_gauss01(rng);
        x[i + 1] = osh_rng_gauss01(rng);
        x[i + 2] = osh_rng_gauss01(rng);
        x[i + 3] = osh_rng_gauss01(rng);
        i += 4;
    }
    while (i < n) {
        x[i] = osh_rng_gauss01(rng);
        i++;
    }
}

void osh_rng_gauss_vec(struct osh_rng *rng, double mu, double sigma, double *restrict x, int n) {
    int i = 0;

    while (i + 3 < n) {
        x[i + 0] = mu + sigma * osh_rng_gauss01(rng);
        x[i + 1] = mu + sigma * osh_rng_gauss01(rng);
        x[i + 2] = mu + sigma * osh_rng_gauss01(rng);
        x[i + 3] = mu + sigma * osh_rng_gauss01(rng);
        i += 4;
    }
    while (i < n) {
        x[i] = mu + sigma * osh_rng_gauss01(rng);
        i++;
    }
}

void osh_rng_u32_vec(struct osh_rng *rng, uint32_t *restrict x, int n) {
    int i = 0;

    while (i + 3 < n) {
        x[i + 0] = osh_rng_u32(rng);
        x[i + 1] = osh_rng_u32(rng);
        x[i + 2] = osh_rng_u32(rng);
        x[i + 3] = osh_rng_u32(rng);
        i += 4;
    }
    while (i < n) {
        x[i] = osh_rng_u32(rng);
        i++;
    }
}

int osh_rng_poisson(struct osh_rng *rng, double lambda) {
    double L;
    double p;
    int k;

    if (lambda <= 0.0) {
        return 0;
    }
    L = exp(-lambda);
    p = 1.0;
    k = 0;
    do {
        ++k;
        p *= osh_rng_double(rng);
    } while (p > L && k <= 64);
    return k - 1;
}
