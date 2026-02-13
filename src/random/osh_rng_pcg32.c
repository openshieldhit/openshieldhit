/*
 * PCG32 random number generator
 *
 * Based on:
 *   M.E. O'Neill, pcg-random.org
 *   Apache License 2.0
 */

#include "random/osh_rng.h"

/*
 * Seed PCG32 state.
 * stream selects an independent sequence (must be distinct across lanes).
 */
void osh_rng_pcg32_init(struct osh_rng *rng, uint64_t seed, uint64_t stream) {
    rng->u.pcg32.state = 0ULL;
    rng->u.pcg32.inc = (stream << 1u) | 1u;

    /* Advance once with inc set */
    osh_rng_pcg32_u32(rng);

    /* Add seed */
    rng->u.pcg32.state += seed;

    /* Advance again to mix seed */
    osh_rng_pcg32_u32(rng); /* test */
}

/*
 * Generate next 32-bit random number (PCG32 XSH RR).
 */
uint32_t osh_rng_pcg32_u32(struct osh_rng *rng) {
    uint64_t oldstate;
    uint32_t xorshifted;
    uint32_t rot;

    oldstate = rng->u.pcg32.state;

    /* Advance internal state */
    rng->u.pcg32.state = oldstate * 6364136223846793005ULL + (rng->u.pcg32.inc | 1ULL);

    /* Output function XSH RR */
    xorshifted = (uint32_t) (((oldstate >> 18u) ^ oldstate) >> 27u);
    rot = (uint32_t) (oldstate >> 59u);

    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}
