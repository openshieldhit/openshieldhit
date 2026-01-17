/*
 * xoshiro256** random number generator (64-bit output)
 *
 * Reference: http://prng.di.unimi.it/
 * (xoshiro/xoroshiro family by David Blackman and Sebastiano Vigna)
 *
 * Seeding uses splitmix64 to expand (seed, stream) into 256-bit state.
 */

#include "osh_rng.h"

/* rotate left */
static uint64_t _rotl64(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }

/* splitmix64: good for seeding other generators */
static uint64_t _splitmix64_next(uint64_t *x) {
    uint64_t z;

    *x += 0x9e3779b97f4a7c15ULL;
    z = *x;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

/*
 * Initialize xoshiro256** state from (seed, stream).
 * stream selects an independent lane by perturbing the splitmix input.
 */
void osh_rng_xoshiro256ss_init(struct osh_rng *rng, uint64_t seed, uint64_t stream) {
    uint64_t x;

    /* Mix seed and stream into the splitmix input.
     * This is not cryptographic; it is only to separate lanes.
     */
    x = seed ^ (stream * 0x9e3779b97f4a7c15ULL);

    rng->u.xoshiro256ss.s[0] = _splitmix64_next(&x);
    rng->u.xoshiro256ss.s[1] = _splitmix64_next(&x);
    rng->u.xoshiro256ss.s[2] = _splitmix64_next(&x);
    rng->u.xoshiro256ss.s[3] = _splitmix64_next(&x);

    /* xoshiro256** requires a non-zero state. splitmix64 virtually guarantees
     * this, but if you want belt-and-suspenders, you could check all-zero here.
     */
}

/*
 * Generate next 64-bit random number (xoshiro256**).
 */
uint64_t osh_rng_xoshiro256ss_u64(struct osh_rng *rng) {
    uint64_t *s;
    uint64_t result;
    uint64_t t;

    s = rng->u.xoshiro256ss.s;

    result = _rotl64(s[1] * 5ULL, 7) * 9ULL;
    t = s[1] << 17;

    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];

    s[2] ^= t;
    s[3] = _rotl64(s[3], 45);

    return result;
}
