#include "osh_rng.h"
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

uint32_t osh_rng_u32(struct osh_rng *rng) {
    switch (rng->type) {
    case OSH_RNG_TYPE_PCG32:
        return osh_rng_pcg32_u32(rng);

    case OSH_RNG_TYPE_XOSHIRO256SS:
        return (uint32_t)(osh_rng_xoshiro256ss_u64(rng) >> 32);

    default:
        return osh_rng_pcg32_u32(rng);
    }
}

uint64_t osh_rng_u64(struct osh_rng *rng) {
    uint64_t hi;
    uint64_t lo;

    switch (rng->type) {
    case OSH_RNG_TYPE_PCG32:
        hi = (uint64_t)osh_rng_pcg32_u32(rng);
        lo = (uint64_t)osh_rng_pcg32_u32(rng);
        return (hi << 32) | lo;

    case OSH_RNG_TYPE_XOSHIRO256SS:
        return osh_rng_xoshiro256ss_u64(rng);

    default:
        return ((uint64_t)osh_rng_u32(rng) << 32) | (uint64_t)osh_rng_u32(rng);
    }
}

float osh_rng_float(struct osh_rng *rng) {
    uint32_t r;
    uint32_t mant;

    /* Use 24 bits for float mantissa */
    r = osh_rng_u32(rng);
    mant = r >> 8; /* top 24 bits */

    return (float)mant * (1.0f / 16777216.0f); /* 2^24 */
}

double osh_rng_double(struct osh_rng *rng) {
    uint64_t r;
    uint64_t mant;

    /* Use 53 bits for double mantissa */
    r = osh_rng_u64(rng);
    mant = r >> 11; /* top 53 bits */

    return (double)mant * (1.0 / 9007199254740992.0); /* 2^53 */
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

double osh_rng_gauss(struct osh_rng *rng, double mu, double sigma) { return mu + sigma * osh_rng_gauss01(rng); }

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