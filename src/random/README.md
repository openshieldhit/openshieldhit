# Random Module (`src/random/`)

OpenShieldHIT's random-number layer: the pseudo-random **engines**, the
**distributions** built on them (uniform, Gaussian, Poisson), and the
**per-history seeding** that makes Monte Carlo transport reproducible and
parallel-ready. Every routine is a pure function of an explicit
`struct osh_rng` passed by pointer — there is **no global RNG state**.

## What is where

| File | Contents |
|---|---|
| `osh_rng.h` | Public API: `struct osh_rng`, engine types, distributions, and the per-history seeding primitives (`osh_rng_seed_history`, `osh_rng_split`). |
| `osh_rng.c` | Engine dispatch, distributions (Marsaglia-polar Gaussian, Knuth Poisson), and the stream-mixing/seeding logic. |
| `osh_rng_pcg32.c` | PCG32 engine (O'Neill) — the default; 16-byte state, built-in stream selector. |
| `osh_rng_xoshiro256ss.c` | xoshiro256\*\* engine (Blackman & Vigna) — 64-bit output, 2²⁵⁶ period. |

## Why it is built this way (in brief)

- **Two small, high-quality engines** (both pass TestU01 BigCrush) because a
  Monte Carlo run needs speed, tiny per-particle state (for SIMD/GPU), and
  splittable streams.
- **`struct osh_rng` is a plain value type** (tagged union, no pointers, 56 B),
  so it embeds in a particle-pool slot and copies to a GPU trivially.
- **Per-history streams** key each particle's randomness to its global index, so
  results do not depend on batch size, thread, or rank — the basis for
  reproducible SIMD/multithread/MPI/GPU transport.

## Full documentation

The detailed treatment — engine motivation and references, the distribution
math, **worked numerical examples** of per-history streams and the bit-mixing
constants, and a **tutorial on adding your own PRNG** — lives in the developer
guide:

➡️ **[`docs/dev/random_numbers.md`](../../docs/dev/random_numbers.md)**

Background and roadmap: [issue #148](https://github.com/openshieldhit/openshieldhit/issues/148)
(per-history RNG streams) and the project
[`TODO.md`](../../TODO.md).
