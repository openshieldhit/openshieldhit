/*
 * xoshiro256** random number generator (64-bit output)
 *
 * Reference: http://prng.di.unimi.it/
 * (xoshiro/xoroshiro family by David Blackman and Sebastiano Vigna)
 *
 * Seeding uses splitmix64 to expand (seed, stream) into 256-bit state.
 * Bodies live in osh_rng_xoshiro256ss_hd.h so device kernels compile the
 * same lines; this file re-exports them with unchanged public signatures.
 */

#include "random/osh_rng.h"
#include "random/osh_rng_xoshiro256ss_hd.h"

/*
 * Initialize xoshiro256** state from (seed, stream).
 * stream selects an independent lane by perturbing the splitmix input.
 */
void osh_rng_xoshiro256ss_init(struct osh_rng *rng, uint64_t seed, uint64_t stream) {
    _osh_rng_xoshiro256ss_init_hd(rng, seed, stream);
}

/*
 * Generate next 64-bit random number (xoshiro256**).
 */
uint64_t osh_rng_xoshiro256ss_u64(struct osh_rng *rng) {
    return _osh_rng_xoshiro256ss_u64_hd(rng);
}
