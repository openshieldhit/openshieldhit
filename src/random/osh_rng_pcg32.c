/*
 * PCG32 random number generator
 *
 * Based on:
 *   M.E. O'Neill, pcg-random.org
 *   Apache License 2.0
 */

#include "random/osh_rng.h"
#include "random/osh_rng_pcg32_hd.h"

/*
 * Seed PCG32 state.
 * stream selects an independent sequence (must be distinct across lanes).
 */
void osh_rng_pcg32_init(struct osh_rng *rng, uint64_t seed, uint64_t stream) {
    _osh_rng_pcg32_init_hd(rng, seed, stream);
}

/*
 * Generate next 32-bit random number (PCG32 XSH RR).
 */
uint32_t osh_rng_pcg32_u32(struct osh_rng *rng) {
    return _osh_rng_pcg32_u32_hd(rng);
}
