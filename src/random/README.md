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

## The quick mental model

Parallel workers are assigned **history-index ranges**, not slices of one shared
random-number sequence.  For example, a run with `nstat = 1000` might give
worker 0 histories `[0, 250)`, worker 1 histories `[250, 500)`, and so on.  A
history may consume a few random draws or many random draws during transport,
but that does not affect any other worker because each primary history is seeded
from its own global index:

```text
history_index = rndoffset + hist_lo + worker_local_index
```

That is why the beam fill path takes the first global history index for the
current fill explicitly.  The beam runtime should not ask "how many primaries
have I generated so far?" from a shared counter when workers are filling
disjoint ranges.  The worker already knows the stable history IDs it owns, and
those IDs are enough to recreate the same primary and transport RNG streams
regardless of pool capacity, thread scheduling, or MPI rank.

## Full documentation

The detailed treatment — engine motivation and references, the distribution
math, **worked numerical examples** of per-history streams and the bit-mixing
constants, and a **tutorial on adding your own PRNG** — lives in the developer
guide:

➡️ **[`docs/dev/random_numbers.md`](../../docs/dev/random_numbers.md)**

Background and roadmap: [issue #148](https://github.com/openshieldhit/openshieldhit/issues/148)
(per-history RNG streams) and the project
[`TODO.md`](../../TODO.md).
