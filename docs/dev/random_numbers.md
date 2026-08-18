# Random numbers in OpenShieldHIT

This page explains the random-number layer (`src/random/`): the pseudo-random
**engines**, the **distributions** built on them, the **per-history seeding**
that makes Monte Carlo transport reproducible and parallel-ready, and — for the
curious or the student — **how to add your own generator**.

It is written for physicists, computer scientists, and engineers. No prior
familiarity with the code base is assumed; every claim is backed by a reference
or a small worked example you can reproduce.

!!! info "Where the code lives"
    Engines and API: `src/random/osh_rng.{h,c}`,
    `src/random/osh_rng_pcg32.c`, `src/random/osh_rng_xoshiro256ss.c`.
    A short orientation README sits next to them in `src/random/README.md`.

---

## 1. Why random numbers matter here

A Monte Carlo dose engine simulates particle histories by sampling random
choices: where a beam particle starts, how far it travels before an
interaction, how it scatters, whether a nuclear reaction occurs. A single
proton history consumes on the order of **10³–10⁴ random draws**, and a run
transports **10⁶–10⁹ histories**. The generator must therefore be:

1. **Fast** — it sits in the innermost loop.
2. **Statistically sound** — correlations show up as systematic dose errors.
3. **Small** — when many histories are in flight at once (SIMD lanes, GPU
   threads), each carries its own generator state.
4. **Splittable** — independent, reproducible sub-streams are the basis of
   parallelism *and* of debugging (a run must be repeatable).

These four requirements drive every choice below.

---

## 2. The engines, and why they were chosen

Two engines are provided and selected at runtime via `enum osh_rng_type`:

| Engine | Output | State | Period | Notes |
|---|---|---|---|---|
| **PCG32** (`OSH_RNG_TYPE_PCG32`, default) | 32-bit | 16 B (`state`, `inc`) | 2⁶⁴ per stream | tiny state, built-in stream selector |
| **xoshiro256\*\*** (`OSH_RNG_TYPE_XOSHIRO256SS`) | 64-bit | 32 B (`s[4]`) | 2²⁵⁶ − 1 | fast 64-bit output, enormous period |

Why these, and not the classic Mersenne Twister (MT19937)?

- **Footprint.** MT19937 needs ~2.5 KB of state. With one generator per live
  particle that is fatal for SIMD/GPU. PCG32 (16 B) and xoshiro (32 B) are
  designed for exactly this regime.
- **Statistical quality.** Both pass the full
  [TestU01 BigCrush](http://simul.iro.umontreal.ca/testu01/tu01.html) battery.
  MT19937 fails the linear-complexity tests, and the legacy `rand()`-style LCGs
  used in older transport codes fail many more.
- **Stream selection.** PCG32 exposes an explicit *sequence selector* (the odd
  increment `inc`); xoshiro is seeded through SplitMix64. Both yield many
  provably distinct sub-sequences from one seed — the foundation of the
  per-history streams in §6.

PCG32 is the default because its 16-byte state is the smallest and its stream
mechanism is the cleanest; xoshiro256\*\* is offered when a native 64-bit output
or a longer period is preferred.

`struct osh_rng` is a plain value type — a tagged `union`, no pointers, no heap
— so it can be embedded in a particle-pool slot, copied, or memcpy'd to a GPU
without ownership concerns. The whole struct is 56 bytes.

---

## 3. The API

```c
/* Engine lifecycle */
void     osh_rng_init(struct osh_rng *r, enum osh_rng_type t, uint64_t seed, uint64_t stream);

/* Raw integers */
uint32_t osh_rng_u32(struct osh_rng *r);
uint64_t osh_rng_u64(struct osh_rng *r);

/* Uniform reals in [0, 1) */
float    osh_rng_float(struct osh_rng *r);   /* 24-bit mantissa */
double   osh_rng_double(struct osh_rng *r);  /* 53-bit mantissa */

/* Distributions */
double   osh_rng_gauss01(struct osh_rng *r);
double   osh_rng_gauss(struct osh_rng *r, double mu, double sigma);
int      osh_rng_poisson(struct osh_rng *r, double lambda);

/* Per-history seeding (see §6) */
void     osh_rng_seed_history(struct osh_rng *r, enum osh_rng_type t,
                              uint64_t seed, uint64_t hist_index, enum osh_rng_purpose p);
void     osh_rng_split(struct osh_rng *child, struct osh_rng const *parent,
                       uint64_t ordinal);
```

**Uniform reals.** `osh_rng_double` takes the top 53 bits of a 64-bit draw and
scales by 2⁻⁵³; `osh_rng_float` uses the top 24 bits and 2⁻²⁴. Using the *high*
bits avoids the low-bit weakness that the divide-by-modulus idiom suffers on
LCG-style generators.

**Gaussian — Marsaglia polar method.** `osh_rng_gauss01` rejection-samples a
point in the unit disk and maps it to two standard normals; the spare is cached
in the struct and returned on the next call. It avoids the `sin`/`cos` of
Box–Muller. `osh_rng_gauss(mu, sigma)` is the affine rescaling used for beam
phase space and Gaussian straggling/scattering.

**Poisson — Knuth's method.** `osh_rng_poisson(lambda)` uses the
product-of-uniforms algorithm, O(λ) on average, with a finite loop guard.
Appropriate for the small λ of nuclear secondary multiplicities.

**Truncated normal — inverse CDF.** `osh_rng_gauss_trunc` draws from
N(μ, σ²) restricted to `[lo, hi]` (used for the `TCUT0` primary-energy window)
with **exactly one** uniform per draw:

```text
x = μ + σ · Φ⁻¹( Φ(α) + u · [Φ(β) − Φ(α)] ),   α = (lo−μ)/σ,  β = (hi−μ)/σ
```

The interval constants `Φ(α)` and the span `Φ(β) − Φ(α)` are precomputed once
per parameter set by `osh_gauss_trunc_prepare` (four `erfc` calls) and reused
per draw; `struct osh_gauss_trunc` is a plain value type holding them.
`Φ⁻¹` is Acklam's rational approximation, ~1.15e-9 relative, with no
Newton/Halley polish — that would pull `erfc` and `exp` into the per-particle
path to refine a value already far below any physical energy resolution.

Nothing loops: there are no draws to discard, so the cost is constant and the
result is exact wherever the window sits. That matters well before the extreme
cases — a window 3 σ off the mean holds 0.14% of the untruncated distribution,
so discarding out-of-window draws would need ~740 of them per accepted sample,
and any bounded retry count would run out and have to fall back to something
that is no longer the requested distribution.

Two numerical details matter and are worth not re-deriving:

- **The span is formed from the two *smaller* CDF values**, choosing the side by
  the sign of `α + β`: `Φ(β) − Φ(α)` when both cuts are low, `Φᶜ(α) − Φᶜ(β)`
  when both are high, where `Φ` and `Φᶜ` each come straight from `erfc` and
  neither is ever computed as `1` minus the other. Subtracting the two large
  values instead cancels away every significant digit when both cuts sit deep in
  the same tail.
- **`p` and `q = 1 − p` are carried separately** into the probit
  (`p_lo + u·span` and `q_lo − u·span`), so each keeps relative accuracy on its
  own side. That is what makes `min(p, q)` exact, folds the upper tail onto the
  lower one, and removes the "mirror the interval into the lower half" trick the
  standard formulation needs as `p → 1`.

A window whose probability mass underflows to zero is flagged `degenerate` and
falls back to a uniform draw over the window — the correct limit for a
vanishingly narrow interval.

The fixed one-deviate cost also buys reproducibility: the stream position after
sampling does not depend on the truncation, so changing `TCUT0` leaves every
subsequent draw in the history — position, direction — untouched, and the two
runs stay comparable under common random numbers.

---

## 4. The seeding model: `seed`, `stream`, `RNDOFFSET`

`osh_rng_init(r, type, seed, stream)` has two 64-bit knobs:

- **`seed`** selects the run — the same seed reproduces it exactly.
- **`stream`** selects an *independent sub-sequence* within that seed. PCG32
  folds it into the odd increment (`inc = (stream << 1) | 1`); xoshiro folds it
  into SplitMix64 seeding. Same seed + different stream ⇒ statistically
  independent sequences.

Users set `seed` through the `RNDSEED` beam-file card; `RNDOFFSET` has no
beam-file card and is set only via the `-N`/`--seedoffset` CLI flag.
`RNDOFFSET` selects an **independent stream family**: `osh_rng_seeding_init()`
hashes `RNDSEED` and `RNDOFFSET` together through the same SplitMix64-style
mixer used for per-history streams — `seed = RNDSEED` unchanged when
`RNDOFFSET` is 0, otherwise `seed = mix(RNDSEED, RNDOFFSET)` — **not** a plain
sum. Two runs sharing `RNDSEED` but using different `RNDOFFSET` values get
statistically independent streams over the *same* history-index range
`[0, nstat)`, so parallel array-job replicas (e.g. the SH12A `generatemc`
convention of consecutive small `-N` values) merge without correlation (issue
#317). A plain sum would instead let two *different* configurations collide —
`(RNDSEED=S, RNDOFFSET=k)` and `(RNDSEED=S+k, RNDOFFSET=0)` would fold to the
same seed and silently produce byte-identical output — so the hash-mix keeps
that down to an ordinary, negligible 64-bit hash coincidence instead of a
guaranteed one. This is a different axis from process/MPI/worker splitting
(see §6), which instead assigns each worker a **disjoint history-index
range** — `RNDOFFSET` does not touch the history index at all.

---

## 5. Reproducibility is a requirement, not a nicety

A dose calculation that cannot be reproduced cannot be validated, debugged, or
trusted in a clinical workflow. "Reproducible" here means: **the same inputs
produce the same result regardless of how the work was scheduled** — pool batch
size, thread count, MPI rank assignment, GPU warp ordering. That is a strong
constraint, and it is what §6 delivers.

---

## 6. Per-history RNG streams

### The problem

A natural but flawed design gives the whole run one or two shared generators and
lets every particle draw from them in whatever order the scheduler runs them.
That couples the random sequence a history sees to the batch size, the order of
compaction, and how many secondaries its neighbours happened to spawn. Such a
design is **not reproducible** the moment you change batch size or add a thread,
and it cannot be split across lanes/threads/ranks without serialising draws or
abandoning reproducibility.

### The ELI5 version

Workers do not reserve chunks of one long random-number sequence.  They reserve
history IDs.

For example, a run with `nstat = 1000` can be split into four worker ranges:

```text
worker 0: histories [0, 250)
worker 1: histories [250, 500)
worker 2: histories [500, 750)
worker 3: histories [750, 1000)
```

History 17 might need ten random draws, while history 18 might need ten thousand
because it scatters more, creates secondaries, or stays in the geometry longer.
That variation is fine because history 18 does not continue where history 17
stopped.  Instead, each primary history starts its own RNG stream from its
global history index:

```text
history_index = worker_range_start + worker_local_index
```

So worker 1's first primary is seeded from history index `250`, regardless of
how many random numbers worker 0 consumed.  The important
parallelism rule is therefore simple: worker ranges must be disjoint and
together cover `[0, nstat)`.  Then each history ID is generated once, and each
history gets the same random stream no matter which worker runs it.

### What "split" and "ordinal" mean

When a nuclear event creates a secondary particle, the secondary also needs
random numbers for its later transport. It must not simply draw from the
parent's stream, because then the parent's future trajectory would depend on
how many children happened to be created, injected, dropped, or transported
first.

The useful mental model is procedural generation. A world seed can derive a
region seed, which can derive a room seed, while the rest of the world remains
unchanged if one room is skipped. Here, a parent particle's current RNG state
derives child streams:

```text
parent particle RNG state
  -> child ordinal 0 RNG state
  -> child ordinal 1 RNG state
  -> child ordinal 2 RNG state
```

An **ordinal** is only the child's position within that one event: ordinal 0 is
the first secondary, ordinal 1 the second, and so on. It is not the value stored
by the child, and it is not a global particle ID. Its job is to give each
sibling a stable address in the event's ordered list of secondaries.

`osh_rng_split(child, parent, ordinal)` therefore does two things:

1. It creates a new RNG state for that child from the parent's **current** RNG
   state plus the child's ordinal.
2. It leaves the real parent RNG untouched.

The implementation hashes the parent's current internal state into a lineage
key, derives the child's `(seed, stream)` pair from that key and the ordinal,
and initializes the child with the normal engine initializer. It reads the
parent's state but never advances it, and — because the key comes from the
parent's raw state words, which the engine permutes before emitting — the child
seed is no longer drawn from, or structurally correlated with, the parent's own
subsequent output (issue #299). The parent continues as if no split had
happened.

That means children do **not** take alternating numbers from the parent's
sequence:

```text
not: parent gets a,c,e,... and child gets b,d,f,...
```

Instead, each child gets its own independent reproducible sequence:

```text
parent continues:       a b c d e ...
child ordinal 0 draws:  x0 x1 x2 ...
child ordinal 1 draws:  y0 y1 y2 ...
child ordinal 2 draws:  z0 z1 z2 ...
```

This is safe across pools. If a recoil proton stays in the ion pool, a neutron
moves to the neutron pool, or a fragment is processed later, the RNG state
travels with that particle. A child that consumes many random numbers cannot
shift its parent or siblings, because it consumes only from its own carried RNG
state.

### The model

Every history owns an **independent stream keyed by its global index**, carried
*with the particle* (one `struct osh_rng` per pool slot) so draws follow the
particle, not the schedule:

- **Primaries** seed from their global history index:
  `stream_id = mix(seed, prim_idx, purpose)`, where `seed` already combines
  `rndseed` and `RNDOFFSET` (§4). A separate `BEAM` purpose seeds a transient
  generator for source sampling, so the source phase space is identical
  regardless of which physics options are enabled.
- **Secondaries** (nuclear recoils, abrasion nucleons, Fermi break-up fragments)
  get their stream by **splitting from the parent** with `osh_rng_split(child,
  parent, ordinal)`. The split seeds the child by hashing the parent's current
  internal state, keyed by the child's `ordinal`; it **does not consume a draw
  from the parent**. The child stream is therefore a pure function of the
  parent's state (itself a pure function of the parent's lineage) and the
  ordinal — reproducible no matter when the secondary is created, and, crucially,
  independent of whether any sibling was injected or dropped (see the note below).

```
primary 0   stream = mix(seed, base+0, PHYSICS)
primary 1   stream = mix(seed, base+1, PHYSICS)
  │  ... transports, draws, then a nuclear event fires ...
  ├─ secondary 1a   stream = split(parent₁, ordinal 0)   ← parent unchanged
  └─ secondary 1b   stream = split(parent₁, ordinal 1)   ← parent unchanged
primary 2   stream = mix(seed, base+2, PHYSICS)
```

!!! note "Splitting is drop-independent (issue #213)"
    Because the split reads but never advances the parent, and each sibling is
    keyed by its ordinal rather than by the *order* it consumes parent draws, a
    secondary that is reordered or silently dropped (e.g. when a pool overflows)
    cannot shift its parent's or its siblings' streams. Reproducibility of the
    whole lineage is thus independent of pool occupancy and the checkpoint batch
    schedule. An earlier consuming `osh_rng_split` coupled the parent's stream to
    how many children were injected, which made scored output depend on the batch
    size when a full primary wavefront overflowed the ion pool.

Each stream is a pure function of `(seed, index, purpose)` or of its parent's
lineage, so the random inputs to a history are **invariant** under pool
capacity, thread, and rank.

### Worked example — primaries (verified output, `seed = 2025`)

First three 64-bit draws of each primary's `PHYSICS` stream:

```
history 0:  0x2AFDBA56F1A9CF64  0x7E5A29C839DC5F5A  0x4CF277989F588CED
history 1:  0x82B0327A775CDE1A  0xC9D750811950770C  0x5E35D2EA3A34B911
history 2:  0x31B9D40392E2FFE6  0x56991ABF09144E1F  0x9E5126F97292991A
```

These depend only on `(seed, index)`. Whether history 1 runs first, last, on
thread 0 or 47, in a pool of 1 or 65536, its sequence is the same — the
order-independence guarantee, locked by `test_seed_history_order_independence`.

### Worked example — a secondary via `split` (verified output, `seed = 2025`)

Primary history 7 takes two physics draws, then a nuclear event spawns a
secondary, seeded by splitting the parent at ordinal 0. The split reads the
parent's state but **never advances it**, so the parent draws the same next
value whether or not the secondary was ever created:

```
parent (hist 7) draws #1,#2:      0.223872  0.725151
child (ordinal 0) first 2 draws:  0.091444  0.313103
parent's next draw, split done:   0.632358
parent's next draw, no split:     0.632358   (identical)
```

Re-run with any pool capacity / thread / rank — or drop the secondary
altogether — and both streams are unchanged:

```
child (ordinal 0) first 2 draws again:  0.091444  0.313103   (identical)
```

The child stream is a pure function of the parent's lineage and its ordinal,
and the parent is independent of whether the child was ever split — see
`test_split_deterministic` and `test_split_nonconsuming_and_drop_independent`.

### What is and isn't guaranteed

Per-history streams make the **random draws** a history consumes independent of
capacity/threads/ranks. Scored output is therefore invariant **up to
floating-point summation order** in the shared scoring accumulators: a bin's
last unit-in-the-last-place can differ because histories deposit in a different
order (~1e-13 relative). Floating-point addition is not associative, so this is
expected; eliminating even that residual requires a deterministic scoring
reduction (per-thread tallies + a fixed-order merge), which is future work. The
`bench::capacity_invariance` test checks agreement within a tight tolerance
across pool capacities 1…4096.

---

## 7. The mixing function and the "magic" constants

`rng_mix_stream(seed, hist_index, purpose)` (in `osh_rng.c`) turns the three
inputs into a well-distributed 64-bit stream id. It multiplies each input by a
large **odd** constant, XOR-combines them, then runs the **SplitMix64
finaliser**. None of the constants are arbitrary:

| Constant | Origin | Role |
|---|---|---|
| `0x9E3779B97F4A7C15` | ⌊2⁶⁴/φ⌋, the golden-ratio ("Fibonacci hashing") constant | spreads the **history index**; consecutive indices land far apart (Weyl spacing) |
| `0xD1B54A32D192ED03` | a recognised 64-bit mixing multiplier (the SplitMix64-based custom-hash constant from the competitive-programming "blowing up `unordered_map`" idiom) | spreads the tiny **purpose** enum to a full-width offset |
| `0xBF58476D1CE4E5B9`, `0x94D049BB133111EB`; shifts `30, 27, 31` | David Stafford's **"Mix13"** finaliser (improves MurmurHash3's `fmix64`) | avalanche: one input bit flips ~half the output bits |

**Why multiply by an odd constant?** An odd number is a *unit* in the ring
ℤ/2⁶⁴ (it has a multiplicative inverse mod 2⁶⁴), so the multiply is a
**bijection**: it cannot map two inputs to the same value, and it scatters the
entropy of a low-magnitude input (a small index, a `0`/`1` purpose) across all
64 bits so it survives the XOR instead of sitting in the low bits. The two
multipliers only need to be **odd and distinct**. Note `0xD1B54A32D192ED03` is
odd but **not prime** — primality is irrelevant; oddness (invertibility) and
good bit diffusion are what matter.

The finaliser constants are the delicate ones: Stafford searched for the
shift/multiply triple with the best avalanche, and SplitMix64 (Vigna's reference
implementation) adopted them. The same finaliser seeds xoshiro in
`osh_rng_xoshiro256ss_init`, so both engines inherit identical lane separation.

### Worked example — mixing step by step (verified, `seed = 2025`)

History 0, purpose `PHYSICS` (= 1). The XOR-combine is
`2025 ^ (0 · golden) ^ (1 · 0xD1B54A32D192ED03)`; the index-0 golden term
vanishes, and `0xD1B54A32D192ED03 ^ 0x7E9 = …EAEA`:

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

## 8. Tutorial: adding your own PRNG

The engine layer is deliberately small so it is easy to drop in a new generator
— for experimentation, teaching, or benchmarking. A generator needs only a
**seed/stream init** and a **raw output** function; everything else
(distributions, per-history seeding, vector helpers) is built on top and comes
for free. Five steps:

**1. Declare the type.** Add an entry to `enum osh_rng_type` in `osh_rng.h`:

```c
enum osh_rng_type {
    OSH_RNG_TYPE_PCG32 = 1,
    OSH_RNG_TYPE_XOSHIRO256SS = 2,
    OSH_RNG_TYPE_MYLCG = 3,        /* new */
};
```

**2. Add state to the union** in `struct osh_rng` (it must be a plain value —
no pointers — to preserve the POD/stack-only contract):

```c
union {
    struct { uint64_t state; uint64_t inc; } pcg32;
    struct { uint64_t s[4]; } xoshiro256ss;
    struct { uint64_t state; } mylcg;   /* new */
} u;
```

**3. Implement the engine** in a new `osh_rng_mylcg.c`. A textbook 64-bit LCG
(constants from Knuth/MMIX) makes a compact teaching example — note it is *not*
recommended for production, it just shows the integration surface:

```c
#include "random/osh_rng.h"

void osh_rng_mylcg_init(struct osh_rng *r, uint64_t seed, uint64_t stream) {
    /* Fold the stream into the seed so different streams diverge. */
    r->u.mylcg.state = seed ^ (stream * 0x9E3779B97F4A7C15ULL);
}

uint64_t osh_rng_mylcg_u64(struct osh_rng *r) {
    r->u.mylcg.state = r->u.mylcg.state * 6364136223846793005ULL + 1442695040888963407ULL;
    return r->u.mylcg.state;
}
```

**4. Wire it into the dispatch** in `osh_rng.c` — three `switch` statements:
`osh_rng_init` (call your init), `osh_rng_u64` (return your output), and
`osh_rng_u32` (e.g. the high 32 bits of your u64). That is all the rest of the
module needs: `osh_rng_double`, `osh_rng_gauss01`, `osh_rng_seed_history`,
`osh_rng_split`, and the vector helpers are all expressed in terms of
`osh_rng_u32`/`osh_rng_u64` and work unchanged.

**5. Build and test.** Add `osh_rng_mylcg.c` to `src/random/CMakeLists.txt`, then
add a known-answer test to `tests/unit/test_osh_rng.c` (seed the engine, print
the first few outputs once, paste them back as the expected values). Run it
through TestU01 if you mean to use it for real work.

!!! warning "Do not ship a toy generator"
    The LCG above is for illustration. A production engine must pass empirical
    test batteries (TestU01) and provide independent streams. If you want
    massive GPU parallelism, a **counter-based** generator (Philox/Threefry,
    [Salmon et al. 2011](https://doi.org/10.1145/2063384.2063405)) is the
    natural next addition: it is stateless — `output = f(key, counter)` — so a
    GPU thread needs no carried state and can skip ahead in O(1).

---

## 9. Current capabilities and roadmap

**Available today**

- Two production engines (PCG32, xoshiro256\*\*), runtime-selectable.
- Uniform, Gaussian, truncated-Gaussian, and Poisson distributions, plus
  batched vector fills.
- Per-history seeding (`osh_rng_seed_history`, `osh_rng_split`) with independent
  `BEAM`/`PHYSICS` sub-streams and lineage-deterministic secondaries.
- One RNG stream carried per particle-pool slot, surviving compaction.
- `RNDOFFSET` hash-combined with the run seed (not added) to select an
  independent stream family across whole runs (parallel array-job replicas)
  without colliding with another run's plain `RNDSEED`; process/MPI/worker
  splitting is the separate, orthogonal disjoint-history-range mechanism.
- Scored-output invariance across pool capacities, up to floating-point
  reduction order (`bench::capacity_invariance`).

**Roadmap** (see the project [TODO](https://github.com/openshieldhit/openshieldhit/blob/main/TODO.md)
and [issue #148](https://github.com/openshieldhit/openshieldhit/issues/148))

- Deterministic scoring reduction for byte-exact capacity invariance.
- An SoA / N-lane RNG draw API for SIMD beam fill.
- A counter-based engine (Philox/Threefry) as a third `osh_rng_type` for GPU.

The broader parallelization sequence is SIMD → multithread → MPI → GPU; the
per-history seeding here is the shared foundation all four build on.

---

## 10. Testing

- `tests/unit/test_osh_rng.c` — known-answer sequences for both engines, uniform
  ranges, Gaussian known values, vector/scalar equivalence, Poisson, and the
  seeding primitives (order-independence, purpose-independence, no aliasing of
  seeded state across a large index range, split determinism).
- `tests/unit/test_osh_rng_gauss_trunc.c` — truncated normal: probit round-trip
  accuracy against `erfc` over ±37 σ, sampled moments against the closed forms
  for five windows (including 8 σ off the mean), a Kolmogorov–Smirnov test
  against the exact truncated CDF, monotonicity and endpoint mapping of the
  quantile transform, fixed one-deviate stream consumption, and the degenerate /
  zero-width / unbounded window paths.
- `tests/bench/run_capacity_invariance.cmake` — end-to-end pool-capacity
  invariance of scored output.

---

## 11. References

- M. E. O'Neill, *PCG: A Family of Simple Fast Space-Efficient Statistically
  Good Algorithms for Random Number Generation*, 2014.
  <https://www.pcg-random.org/>
- D. Blackman, S. Vigna, *Scrambled Linear Pseudorandom Number Generators*, ACM
  TOMS 47(4), 2021. <https://prng.di.unimi.it/>
- G. L. Steele, D. Lea, C. H. Flood, *Fast Splittable Pseudorandom Number
  Generators*, OOPSLA 2014 (SplitMix).
- D. Stafford, *Better Bit Mixing — Improving on MurmurHash3's 64-bit
  Finalizer*, 2011.
  <http://zimbry.blogspot.com/2011/09/better-bit-mixing-improving-on.html>
- D. E. Knuth, *The Art of Computer Programming*, Vol. 2 (random numbers,
  Poisson) and Vol. 3 §6.4 (multiplicative/Fibonacci hashing).
- G. Marsaglia, T. A. Bray, *A Convenient Method for Generating Normal
  Variables*, SIAM Review 6(3), 1964 (polar method).
- P. J. Acklam, *An algorithm for computing the inverse normal cumulative
  distribution function*, unpublished note (2000–2010) — the probit rational
  approximation and its coefficients, stated as relative error < 1.15e-9 over
  the whole range. Acklam's own page has been offline for years; the algorithm
  and its constants are preserved and explained in L. M. Barros, *Acklam's
  Algorithm for the Inverse Normal CDF*, 2017.
  <https://stackedboxes.org/2017/05/01/acklams-normal-quantile-function/>
- P. L'Ecuyer, R. Simard, *TestU01: A C Library for Empirical Testing of Random
  Number Generators*, ACM TOMS 33(4), 2007.
- J. K. Salmon, M. A. Moraes, R. O. Dror, D. E. Shaw, *Parallel Random Numbers:
  As Easy as 1, 2, 3*, SC 2011. <https://doi.org/10.1145/2063384.2063405>
