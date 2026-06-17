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
 * rng_mix_stream() — derive a 64-bit stream id from (seed, history index,
 * purpose).  It hash-combines the three inputs, then runs the SplitMix64
 * finaliser.  Being a pure function of its inputs is the whole point: a
 * history maps to the same stream id no matter when, on which thread, or in
 * which pool capacity it is seeded.
 *
 * Why each constant — these are standard mixing constants with known good
 * diffusion, not arbitrary values:
 *
 *   Step 1 — spread each input across all 64 bits, then XOR-combine.
 *   Each field is multiplied by a large ODD constant.  An odd constant is a
 *   unit in the ring Z/2^64, so the multiply is a bijection (it cannot lose
 *   information) and it scatters the entropy of low-magnitude inputs — a small
 *   history index, a 0/1 purpose — across the full word so they survive the
 *   XOR instead of sitting in the low bits.  The constants only need to be odd
 *   and distinct so the index and purpose axes do not alias:
 *     - 0x9E3779B97F4A7C15 = floor(2^64 / phi), the 64-bit golden-ratio
 *       constant ("Fibonacci hashing").  Consecutive indices land far apart —
 *       the classic Weyl-sequence spacing.  Refs: Knuth, TAOCP Vol. 3 sec. 6.4
 *       (multiplicative hashing); Steele, Lea & Flood, "Fast Splittable
 *       Pseudorandom Number Generators", OOPSLA 2014 (SplitMix's GOLDEN_GAMMA).
 *     - 0xD1B54A32D192ED03 is a second large odd constant applied to the
 *       (tiny) purpose enum so BEAM and PHYSICS get well-separated, full-width
 *       offsets; its only required property is oddness (invertible multiply).
 *
 *   Step 2 — the SplitMix64 finaliser (David Stafford's "Mix13" variant), an
 *   improvement on MurmurHash3's fmix64 with stronger avalanche: flipping one
 *   input bit flips ~half the output bits with low bias.  The multipliers
 *   (0xBF58476D1CE4E5B9, 0x94D049BB133111EB) and shift amounts (30, 27, 31)
 *   are taken verbatim.  Each xorshift-multiply step is invertible, so the
 *   finaliser is a bijection.  Refs: D. Stafford, "Better Bit Mixing —
 *   Improving on MurmurHash3's 64-bit Finalizer" (2011),
 *   http://zimbry.blogspot.com/2011/09/better-bit-mixing-improving-on.html ;
 *   public-domain reference implementation https://prng.di.unimi.it/splitmix64.c
 *   This is the same finaliser used to seed xoshiro256** in
 *   osh_rng_xoshiro256ss_init(), so both engines get identical lane
 *   separation.  Disjointness over the index ranges we actually use is also
 *   checked empirically by test_seed_history_disjoint_ranges.
 */
static uint64_t rng_mix_stream(uint64_t seed, uint64_t hist_index, uint64_t purpose) {
    /* Golden-ratio constant on the index axis; a second odd constant on the
     * purpose axis (see the constant rationale above). */
    uint64_t x = seed ^ (hist_index * 0x9E3779B97F4A7C15ULL) ^ (purpose * 0xD1B54A32D192ED03ULL);

    /* Stafford Mix13 finaliser: xorshift / multiply / xorshift / multiply /
     * xorshift, giving full 64-bit avalanche. */
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
