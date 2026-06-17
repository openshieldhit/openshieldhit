# Random Module (`src/random/`)

This module is OpenShieldHIT's random-number layer: the pseudo-random
**engines**, the **distributions** built on top of them (uniform, Gaussian,
Poisson), and the **per-history seeding** primitives that make Monte Carlo
transport reproducible and parallel-ready.

It has no dependency on the rest of the simulation — every routine is a pure
function of an explicit `struct osh_rng` passed by pointer. There is **no
global RNG state** anywhere in the library.

```
random/
    osh_rng.h               public API: engines, distributions, seeding
    osh_rng.c               engine dispatch, distributions, mixing/seeding
    osh_rng_pcg32.c         PCG32 engine (O'Neill)
    osh_rng_xoshiro256ss.c  xoshiro256** engine (Blackman & Vigna)
```

---

## 1. Design goals and why these engines

A Monte Carlo dose engine draws on the order of **10³–10⁴ random numbers per
primary** and transports 10⁶–10⁹ primaries, so the generator must be (a) fast,
(b) statistically sound, (c) cheap to hold per particle, and (d) splittable
into independent streams for parallelism. Two engines are provided and selected
at runtime (`enum osh_rng_type`):

| Engine | Output | State | Period | Role |
|---|---|---|---|---|
| **PCG32** (`OSH_RNG_TYPE_PCG32`) | 32-bit | 16 bytes (`state`, `inc`) | 2⁶⁴ per stream | default; tiny state, built-in stream selector |
| **xoshiro256\*\*** (`OSH_RNG_TYPE_XOSHIRO256SS`) | 64-bit | 32 bytes (`s[4]`) | 2²⁵⁶−1 | fast 64-bit output, very long period |

Why these two rather than, say, Mersenne Twister:

- **Speed and footprint.** Both are a handful of arithmetic/shift ops with no
  large state table. MT19937 needs 2.5 KB of state per stream — unworkable when
  every live particle in a GPU/SIMD wavefront carries its own generator. PCG32's
  16-byte state and xoshiro's 32-byte state are designed exactly for that.
- **Statistical quality.** Both pass the full
  [TestU01 BigCrush](http://simul.iro.umontreal.ca/testu01/tu01.html) battery,
  unlike MT19937 (which fails linear-complexity tests) and unlike the legacy
  `rand()` LCGs used in older transport codes.
- **Stream selection for parallelism.** PCG32 exposes an explicit *sequence
  selector* (the odd increment `inc`); xoshiro is seeded through SplitMix64.
  Both give us many provably-distinct sub-sequences from one seed — the
  foundation of the per-history streams in §4.

`struct osh_rng` is a plain stack value (a POD union, no pointers, no heap), so
it can be embedded in a particle pool slot, copied, memcpy'd to a GPU, or held
in an SoA lane array without any ownership concerns.

> **State layout** (`osh_rng.h`): `type` tag + a `union { pcg32; xoshiro256ss; }`
> + a cached Gaussian spare. The whole struct is 56 bytes.

---

## 2. Public API

```c
/* Engine lifecycle */
void     osh_rng_init(struct osh_rng *r, enum osh_rng_type t, uint64_t seed, uint64_t stream);

/* Raw integers */
uint32_t osh_rng_u32(struct osh_rng *r);
uint64_t osh_rng_u64(struct osh_rng *r);

/* Uniform reals in [0, 1) */
float    osh_rng_float(struct osh_rng *r);    /* 24-bit mantissa */
double   osh_rng_double(struct osh_rng *r);   /* 53-bit mantissa */

/* Distributions */
double   osh_rng_gauss01(struct osh_rng *r);                 /* N(0,1)      */
double   osh_rng_gauss(struct osh_rng *r, double mu, double sigma);
int      osh_rng_poisson(struct osh_rng *r, double lambda);

/* Batched fills (same sequence as repeated scalar calls) */
void     osh_rng_double_vec/float_vec/u32_vec/gauss01_vec/gauss_vec(...);

/* Per-history seeding (see §4) */
void     osh_rng_seed_history(struct osh_rng *r, enum osh_rng_type t,
                              uint64_t seed, uint64_t hist_index, enum osh_rng_purpose p);
void     osh_rng_split(struct osh_rng *child, struct osh_rng *parent);
```

### Uniform reals

`osh_rng_double` takes the top 53 bits of a 64-bit draw and scales by 2⁻⁵³,
giving a uniform value in `[0, 1)` with full double-precision mantissa
resolution; `osh_rng_float` uses the top 24 bits and 2⁻²⁴. This "high-bits"
construction avoids the low-bit weakness that the division-by-modulus idiom
suffers on LCG-style generators.

### Gaussian — Marsaglia polar method

`osh_rng_gauss01` uses the **Marsaglia polar method**: rejection-sample a point
in the unit disk, then map it to two independent standard normals. The method
produces two variates per evaluation, so the spare is cached in the struct and
returned on the next call. It is branch-light and avoids the `sin`/`cos` of
Box–Muller. `osh_rng_gauss(mu, sigma)` is the obvious affine rescaling, used for
beam phase-space sampling and Gaussian energy straggling / multiple scattering.

### Poisson — Knuth's method

`osh_rng_poisson(lambda)` uses Knuth's product-of-uniforms algorithm, O(λ) on
average, with a finite loop guard for pathological inputs. Suitable for the
modest λ encountered in nuclear secondary multiplicities; callers needing large
λ should switch to a transformed-rejection sampler.

---

## 3. Seeding model: `seed`, `stream`, and `rndoffset`

`osh_rng_init(r, type, seed, stream)` takes two 64-bit knobs:

- **`seed`** selects the run. The same seed reproduces the run bit-for-bit.
- **`stream`** selects an *independent sub-sequence* within that seed. PCG32
  folds it into the odd increment (`inc = (stream<<1)|1`); xoshiro folds it into
  the SplitMix64 seeding. Two RNGs with the same seed but different streams
  produce statistically independent sequences.

The transport layer exposes these to users as the `RNDSEED` and `RNDOFFSET`
beam-file cards (`--pool-capacity` is unrelated — that is a performance knob).
`RNDOFFSET` is the **global history-index base** (see §4): giving each process
or MPI rank a disjoint offset range yields disjoint, non-overlapping streams, so
independent runs can be merged without correlation.

---

## 4. Per-history RNG streams (issue #148 / this PR)

### The problem this solves

Previously transport used **two run-wide streams** — one for beam sampling, one
for physics — and *every* particle in the wavefront drew from the shared
physics stream in scheduling order. That coupled the random sequence a history
saw to the pool capacity, the compaction order, and how many secondaries its
siblings happened to spawn. Such coupling makes results non-reproducible the
moment you change batch size or add a thread, and it cannot be split across
SIMD lanes / threads / GPU warps / MPI ranks without either serialising the
draws or giving up reproducibility.

### The model

Every particle history owns an **independent stream keyed by its global
index**, carried *with the particle* (one `struct osh_rng` per pool slot) so
that draws follow the particle, not the schedule:

- **Primaries** are seeded from their global history index:
  `stream_id = mix(rndseed, RNDOFFSET + prim_idx, purpose)`, and the slot's RNG
  is initialised on that stream. A separate `BEAM` purpose seeds a transient
  generator that samples the source phase space, so source sampling is identical
  regardless of which physics options (NUCRE/MSCAT/STRAGG) are enabled.
- **Secondaries** (nuclear recoils, abrasion nucleons, Fermi break-up fragments)
  get their stream by **splitting from the parent** with `osh_rng_split`, which
  consumes two draws from the parent to seed the child. The parent's draw
  sequence is deterministic along its own lineage, so the child stream is
  reproducible no matter when the secondary is created.

```
primary 0   stream = mix(seed, base+0, PHYSICS)
primary 1   stream = mix(seed, base+1, PHYSICS)
  │  ... transports, draws, then a nuclear event fires ...
  ├─ secondary 1a   stream = split(parent₁)      ← 1st split (parent draws #k,#k+1)
  └─ secondary 1b   stream = split(parent₁)      ← 2nd split (parent draws #k+2,#k+3)
primary 2   stream = mix(seed, base+2, PHYSICS)
```

Because each stream is a pure function of `(seed, index, purpose)` or of its
parent's lineage, the random inputs to a history are **invariant** under pool
capacity, thread assignment, and rank — which is exactly what makes
reproducible parallelism possible.

### Worked example — primaries (verified output, `seed = 2025`)

First three 64-bit draws of each primary's `PHYSICS` stream:

```
history 0:  0x2AFDBA56F1A9CF64  0x7E5A29C839DC5F5A  0x4CF277989F588CED
history 1:  0x82B0327A775CDE1A  0xC9D750811950770C  0x5E35D2EA3A34B911
history 2:  0x31B9D40392E2FFE6  0x56991ABF09144E1F  0x9E5126F97292991A
```

These values depend only on `(seed, history index)`. Whether history 1 is
transported first, last, on thread 0 or thread 47, in a pool of 1 or 65536, its
sequence is the same — that is the order-independence guarantee, locked by the
`test_seed_history_order_independence` unit test.

### Worked example — a secondary via `split` (verified output)

Primary history 7 takes two physics draws, then a nuclear event spawns a
secondary, which is seeded by splitting the parent:

```
parent (hist 7) draws #1,#2:  0.223872  0.725151
child  first 2 draws:         0.252732  0.955397
parent next draw (after the split advanced it by two): 0.609548
```

Re-running with a different pool capacity, thread, or rank reproduces the child
exactly:

```
child  first 2 draws again:   0.252732  0.955397   (identical)
```

The child is independent of the parent's continuation (`0.252732…` vs
`0.609548…`) yet fully determined by the parent's lineage — see
`test_split_deterministic`.

### What is and isn't guaranteed

Per-history streams make the **random draws** a history consumes independent of
capacity/threads/ranks. Scored output is therefore capacity-independent **up to
floating-point summation order** in the shared scoring accumulators (a bin's
last ULP can differ because histories deposit in a different order). Removing
even that residual requires a deterministic scoring reduction (per-thread
tallies + fixed-order merge), tracked as the follow-up to #148. The
`bench::capacity_invariance` test checks agreement within a tight tolerance
across capacities 1…4096.

---

## 5. The mixing function and the "magic" hex constants

`rng_mix_stream(seed, hist_index, purpose)` (in `osh_rng.c`) turns the three
inputs into a well-distributed 64-bit stream id. It first multiplies each input
by a large **odd** constant and XOR-combines them, then runs the **SplitMix64
finaliser**. None of the constants are arbitrary:

| Constant | Name / origin | Role |
|---|---|---|
| `0x9E3779B97F4A7C15` | ⌊2⁶⁴/φ⌋, the golden-ratio ("Fibonacci hashing") constant | spreads the **history index**; consecutive indices land far apart (Weyl spacing) |
| `0xD1B54A32D192ED03` | a recognised 64-bit mixing multiplier (the SplitMix64-based custom-hash constant from the competitive-programming "blowing up `unordered_map`" idiom) | spreads the tiny **purpose** enum to a full-width offset |
| `0xBF58476D1CE4E5B9`, `0x94D049BB133111EB` + shifts `30, 27, 31` | David Stafford's **"Mix13"** 64-bit finaliser (improves MurmurHash3's `fmix64`) | avalanche: one input bit flips ~half the output bits |

**Why multiply by an odd constant?** An odd number is a *unit* in the ring
ℤ/2⁶⁴ (it has a multiplicative inverse mod 2⁶⁴), so multiplication by it is a
**bijection** — it cannot collapse two inputs together, and it scatters the
entropy of a low-magnitude input (a small history index, a `0`/`1` purpose)
across all 64 bits so it survives the XOR instead of sitting in the low bits.
The two multipliers only need to be **odd and distinct** so the index and
purpose axes do not alias. Note `0xD1B54A32D192ED03` is odd but **not prime** —
primality is irrelevant here; oddness (invertibility) and good bit diffusion are
what matter.

The finaliser constants are the more delicate ones: Stafford searched for the
shift/multiply triple with the best avalanche, and SplitMix64 (Vigna's reference
implementation) adopted them verbatim. The same finaliser seeds xoshiro in
`osh_rng_xoshiro256ss_init`, so both engines inherit identical lane separation.

### Worked example — mixing step by step (verified, `seed = 2025`)

History 0, purpose `PHYSICS` (= 1). The XOR-combine is
`2025 ^ (0 · golden) ^ (1 · 0xD1B54A32D192ED03)`; since the history index is 0
its golden term vanishes and we get `0xD1B54A32D192ED03 ^ 0x7E9 = …EAEA`:

```
  after XOR-combine:   0xD1B54A32D192EAEA
  after 1st mul step:  0xE2D45E5B04E6CED9    (^>>30, ×0xBF58476D1CE4E5B9)
  after 2nd mul step:  0x55756236BEDE89D3    (^>>27, ×0x94D049BB133111EB)
  final stream id:     0x5575623614344DBE    (^>>31)
```

Changing only the **purpose** produces a completely different stream id (full
avalanche from a one-field change), which is why BEAM and PHYSICS never
correlate for the same history:

```
  history 0, BEAM    →  0x0040F0F2B7D74795
  history 0, PHYSICS →  0x5575623614344DBE
```

---

## 6. Status in this PR

- ✅ Per-history seeding (`osh_rng_seed_history`, `osh_rng_split`,
  `struct osh_rng_seeding`, `enum osh_rng_purpose`).
- ✅ One RNG stream per particle-pool slot, carried through compaction;
  `prim_idx` widened to `uint64_t` (the seeding key).
- ✅ `RNDOFFSET` redefined as the global history-index base (disjoint ranges for
  process/MPI splitting).
- ✅ Capacity made a runtime knob; invariance covered by
  `bench::capacity_invariance`.
- ✅ Unit tests for the seeding primitives and the engine/distribution surface.

**Not yet done (future work, beyond this PR):**

- Deterministic scoring reduction for byte-exact capacity invariance.
- An SoA / N-lane RNG draw API for SIMD beam fill.
- A **counter-based** engine (Philox/Threefry) as a third `osh_rng_type` for
  massive GPU parallelism and O(1) skip-ahead — these are stateless `f(key,
  counter)` generators, ideal when millions of GPU threads each need an
  independent stream without carrying state.

See `../../TODO.md` and issue #148 for the broader parallelization roadmap
(SIMD → multithread → MPI → GPU).

---

## 7. Examples

Independent streams for a small batch of histories:

```c
#include "random/osh_rng.h"

struct osh_rng_seeding seeding = {
    .type = OSH_RNG_TYPE_PCG32,
    .seed = 2025u,        /* RNDSEED  */
    .hist_base = 0u,      /* RNDOFFSET */
};

for (uint64_t h = 0; h < nstat; ++h) {
    struct osh_rng rng;
    osh_rng_seed_history(&rng, seeding.type, seeding.seed,
                         seeding.hist_base + h, OSH_RNG_PURPOSE_PHYSICS);
    /* ... transport history h, drawing from rng ... */
}
```

Spawning a reproducible secondary stream:

```c
/* parent is the slot's RNG; child is the new secondary's slot RNG */
osh_rng_split(&pool->rng[child_slot], &pool->rng[parent_slot]);
```

A standalone generator for ad-hoc use:

```c
struct osh_rng r;
osh_rng_init(&r, OSH_RNG_TYPE_XOSHIRO256SS, /*seed*/ 42u, /*stream*/ 0u);
double x = osh_rng_double(&r);          /* uniform [0,1)  */
double g = osh_rng_gauss(&r, 0.0, 1.0); /* standard normal */
```

---

## 8. Testing

- `tests/unit/test_osh_rng.c` — known-answer sequences for both engines, uniform
  ranges, Gaussian known values, vector-helper/scalar equivalence, Poisson, and
  the seeding primitives (order-independence, purpose-independence, disjoint
  ranges, split determinism).
- `tests/bench/run_capacity_invariance.cmake` — end-to-end pool-capacity
  invariance of scored output.

---

## 9. References

- M. E. O'Neill, *PCG: A Family of Simple Fast Space-Efficient Statistically
  Good Algorithms for Random Number Generation*, 2014.
  <https://www.pcg-random.org/>
- D. Blackman, S. Vigna, *Scrambled Linear Pseudorandom Number Generators*, ACM
  Trans. Math. Softw. 47(4), 2021. <https://prng.di.unimi.it/>
- G. L. Steele, D. Lea, C. H. Flood, *Fast Splittable Pseudorandom Number
  Generators*, OOPSLA 2014 (SplitMix).
- D. Stafford, *Better Bit Mixing — Improving on MurmurHash3's 64-bit
  Finalizer*, 2011.
  <http://zimbry.blogspot.com/2011/09/better-bit-mixing-improving-on.html>
- D. E. Knuth, *The Art of Computer Programming*, Vol. 2 (random numbers, Poisson
  sampling) and Vol. 3 §6.4 (multiplicative/Fibonacci hashing).
- G. Marsaglia, T. A. Bray, *A Convenient Method for Generating Normal
  Variables*, SIAM Review 6(3), 1964 (polar method).
- P. L'Ecuyer, R. Simard, *TestU01: A C Library for Empirical Testing of Random
  Number Generators*, ACM TOMS 33(4), 2007.
- J. K. Salmon et al., *Parallel Random Numbers: As Easy as 1, 2, 3*, SC 2011
  (Philox/Threefry counter-based generators — relevant to the GPU roadmap).
