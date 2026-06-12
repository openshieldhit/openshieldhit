# Algorithmic performance improvements: analysis, prototype, results

*Status: prototype + report, June 2026.  Builds on the profiling baseline of
PR #139 ("Profiling baseline and reusable benchmark harness").  All code
referenced here lives on the branch `claude/algorithmic-performance-improvements-eaik1x`.*

This document is written to be read end-to-end: it explains not just *what*
was changed but *how the problem was found, why the chosen algorithm is
correct, and how the claim "same physics, much faster" was verified*.

---

## 1. Starting point: what the profiling PR already told us

PR #139 committed deterministic callgrind tables for four benchmark cases and
ranked the hotspots.  The headline numbers:

| Rank | Hotspot | Share of instructions |
|------|---------|----------------------|
| 1 | GEMCA `eval_distance` + `get_zone` | 19–31 % on analytic cases, **88 % combined at ~1000 zones** (C8: `eval_distance` 70.9 % + `get_zone` 17.6 %) |
| 2 | libm transcendentals (log/pow/exp/cbrt in Bethe, Tripathi, Molière, straggling) | ~26–27 % on every physics-bearing case |
| 3 | `osh_transport_ion_step` | 8–11 % |
| 4 | Scoring (`score_step` + raytrace) | 9–13 % |
| 5 | RNG | 2.5–4 % |

The phrase to focus on is *"scales catastrophically with zone count"*.  A 26 %
flat cost (libm) can at best give a 1.35× speedup if completely removed.  A
cost that **grows with problem size** has no such ceiling — and complex
geometries are exactly where a Monte Carlo code spends days, not minutes.  So
the first target is asymptotics, not micro-optimization.

### Reproducing the scaling problem

The C8 generator from PR #139 (`tools/bench/gen_gemca_stress.py`) builds an
n×n×n lattice of PMMA spheres carved out of a water cylinder.  Crucially, the
water "matrix" zone is the CSG expression `+target −s1 −s2 … −sN` — one
negative term per sphere.  Baseline wall time (release build, 300 histories,
150 MeV protons, 4-core Xeon @ 2.10 GHz container):

| Zones (spheres + 3) | Baseline wall time |
|---------------------|--------------------|
| 128  (125)  | 3.4 s |
| 346  (343)  | 7.6 s |
| 732  (729)  | 16.3 s |
| 2200 (2197) | 49.4 s |
| 4099 (4096) | 114.2 s |

Almost perfectly linear in zone count, with a ~1 s flat part.  Doubling the
geometry detail doubles the run time even though each particle only ever
interacts with a handful of nearby spheres.  That observation — *most of the
work consults bodies that are provably irrelevant to the query point* — is
the entire optimization in one sentence.

---

## 2. Why the hot path is O(N): reading the runtime

The compiled geometry runtime (`src/gemca/runtime/osh_gemca_runtime.c`)
answers two queries per particle per wavefront iteration:

1. **`osh_gemca_runtime_get_zone(ray)`** — which zone contains this point?
   Implementation: walk `zones[0..N)` in order, evaluate each zone's RPN
   (reverse Polish notation) CSG program, return the first match.  Two nested
   O(N) problems hide here:
   - the *scan* visits every zone until a hit (O(N) zones for a point in the
     matrix zone, which is defined last);
   - evaluating the matrix zone itself executes its full RPN program — with
     2·N+1 instructions and N body membership tests — even though at most one
     or two spheres are anywhere near the point.

2. **`eval_distance(zone, ray)`** — how far to the next boundary of the zone
   I'm in?  Implementation: execute the same RPN program with (distance,
   inside) pairs, combining children per Roth's CSG ray-casting rules.  Every
   `PUSH_BODY` computes the exact ray/quadric intersection distance for its
   body — all N spheres, every call, at ~30–60 flops each (a transform plus a
   quadratic solve).

So one transport step of one particle in the matrix zone costs ~3·N body
evaluations.  The existing `GEMCA_RT_GUARD_BODY` optimization (O(1) rejection
of a zone whose *guard* body misses the point) helps zones that are small,
but cannot help the matrix zone — its guard (the target cylinder) contains
essentially every point.

### The key structural observation

In `eval_distance`, **all three CSG operators combine child distances the
same way**: `minpos(d_left, d_right)` (see the `GEMCA_RT_UNION /
INTERSECT / DIFF` cases).  Leaf distances are strictly positive or +∞ by
construction (`dist_body_rt` only keeps `d > 0`).  A fold of an associative,
commutative `min` over a tree is just the minimum of the leaves:

> *the distance to the next zone boundary equals the plain minimum of the
> boundary distances of the zone's leaf bodies — the CSG tree shape is
> irrelevant for the distance part.*

Only the *inside/outside* part needs the boolean structure.  This decouples
the two halves of `eval_distance` and is what makes spatial acceleration
straightforward to retrofit: a "minimum over leaves" is exactly the query
that grids and BVHs accelerate.

(Why does the step loop stay correct if some leaf distances are candidates
that don't actually exit the zone?  Because `get_distance` is a
*re-evaluating stepper*: it advances to the nearest candidate, re-checks
membership, and keeps going.  Extra candidates cost iterations; *missing* a
real boundary would be a correctness bug.  Every real zone boundary point
lies on the surface of some leaf body, inside that body's bounding box —
that's the invariant the culling below preserves.)

---

## 3. The prototype: a uniform grid with cell-specialized CSG programs

New file `src/gemca/runtime/osh_gemca_runtime_accel.{c,h}` (setup time), hot
paths in `osh_gemca_runtime.c`.  Four layers, each independently testable:

### 3.1 Per-body AABBs (`body_universe_aabb`)

Every body is an intersection of half-spaces in its local frame, so each
surface can only *tighten* an axis-aligned bounding box: spheres/ellipsoids
bound all three axes, a CYLZ bounds x/y, axis planes bound one side of one
axis, a cone bounds x/y once z is bounded by accompanying planes.  The local
box is mapped to the universe frame through the inverse body transform
(transpose of the rotation — checked for orthonormality, with "unbounded" as
the safe fallback).  Boxes are inflated by `1e-6 + 1e-9·|coordinate|` so that
the membership epsilon (`OSH_GEMCA_SMALL` on the implicit function) can never
put an "inside" point outside its body's box.  Unbounded bodies (infinite
planes, unknown transforms) are flagged and **never culled** — conservative
by construction.

### 3.2 Zone AABBs and the grid

Zone boxes come from folding the RPN program over interval arithmetic:
union → box union, intersection → box intersection, difference → left box.
A uniform grid covers the union of all bounded zone boxes (~8 cells per
zone, ≤128 cells per axis).  Zones with unbounded boxes go to a global
"always check" list.

### 3.3 Per-cell *constant-folded* zone programs (`specialize_program`)

This is the part that breaks the second O(N), and the part a plain
"candidate zone list" grid would miss: even with a perfect zone candidate
list, *evaluating the matrix zone's own membership* still walks all N sphere
terms.  Instead, for every (cell, zone) pair the zone's RPN program is
re-compiled at setup time with cell-specific knowledge:

- `PUSH body` where the body's box cannot touch the cell → constant
  **outside** (this is sound for *every* point in the cell, which is why it
  can be done once at setup);
- boolean operators fold: `x ∪ 0 = x`, `x ∩ 0 = 0`, `x − 0 = x`, `0 − x = 0`
  — with rollback of instruction ranges that a hard 0 makes unreachable;
- a guard whose body misses the cell, or a program folding to constant 0,
  drops the zone from that cell's candidate list entirely.

For the 729-sphere matrix zone (1459 instructions), a typical lattice cell's
specialized program is `GUARD target, PUSH target, PUSH s_i, DIFF, PUSH s_j,
DIFF` — about **7 instructions instead of 1459**.  Membership cost now scales
with *cell occupancy*, not body count.

Zone lookup merges the cell's candidate list with the "always" list in
ascending zone order, preserving the linear scan's first-match-wins
semantics exactly (important if zones overlap — the parser does not forbid
it).

### 3.4 Distance: flat leaf minimum + 3D-DDA grid walk

Using the `min`-over-leaves identity from §2:

- unbounded leaf bodies (never cullable) are evaluated once, upfront;
- then an Amanatides–Woo 3D-DDA walks the grid cells pierced by the ray and
  evaluates only the bodies pushed by *this zone's specialized program of the
  visited cell*, each pre-filtered by a ray/slab test against the running
  minimum;
- the walk terminates as soon as the running minimum is closer than the
  current cell's exit (no later cell can produce a closer surface), or when
  the ray exits the grid (the zone is bounded, so its boundary was passed).

Per query this is O(cells visited × bodies per cell) ≈ O(1) in total body
count.  Zones with fewer than `OSH_GEMCA_ACCEL_DIST_MIN_LEAVES` (12) leaf
bodies keep the original RPN walk — for a two-body zone, grid bookkeeping
costs more than it saves (measured, see §5).

### What deliberately did *not* change

- The cold parser/compiler, the voxel (DICOM) path, and the AVX2 batch
  entry point are untouched; voxel zones and tiny geometries take the
  original code paths.
- The accelerator is optional at runtime: it can fail to build (allocation,
  pathological geometry) and everything falls back to the linear scan;
  `OSH_GEMCA_ACCEL=0` disables it for A/B comparison.
- All accel data is immutable after build — important later for
  parallelization (#138): worker threads can share it read-only.

---

## 4. Verification: "faster" is meaningless if the physics moved

Methodology, in increasing order of strictness:

1. **Unit/integration tests**: full `ctest` suite before and after — same 54
   passes, same 3 pre-existing failures (DICOM fixtures needing Git LFS,
   identical on the unmodified tree).  This includes
   `test_osh_gemca_runtime_zone`, which cross-checks the AVX2 batch zone
   lookup (accel-free path) against the scalar lookup (accel path) on random
   rays — an accidental but effective consistency test between the two
   implementations.
2. **Bit-identical outputs**: same binary, same seed, accel on vs off, for
   125 / 729 / 2197 / 4096 spheres at 300 histories and 729 spheres at 5000
   histories.  The BDO files are byte-identical except the embedded
   wall-clock timestamp (verified to be the only difference between *any*
   two runs, per PR #139); with the timestamp masked, SHA-256 digests match
   in all five comparisons.
3. **Semantics argument** (because tests only sample): every culling
   decision is backed by a conservative invariant — inflated boxes contain
   their bodies' epsilon-thickened surfaces; a bounded zone absent from a
   cell's candidates provably cannot contain any point of that cell; the DDA
   stop condition only prunes cells whose nearest possible surface is beyond
   the current minimum.

One honest footnote on bit-identity: `dist_body_rt` computes distances to
*infinite* surface primitives (e.g. the infinite cylinder of an RCC), so the
original code sometimes steps to phantom intersections outside the body.
AABB culling removes those phantom candidates.  The *total* distance is
identical as a real number (the re-evaluating stepper lands on the same
boundary), but the partial sums can in principle round differently for
multi-surface bodies.  In all cases measured here the outputs were
bit-identical (sphere lattices: a sphere's primitive *is* its surface).  If a
future case shows ulp-level differences, that is the place to look first —
it is a candidate-set change, not a physics change.

---

## 5. Results

Same binary, same machine, `OSH_GEMCA_ACCEL=0` (off) vs default (on),
300 histories, single runs (run-to-run noise ~±5 % in this container):

| Spheres | Off | On | Speedup |
|---------|-----|----|---------|
| 125  | 4.3 s | **0.37 s** | 12× |
| 343  | —¹ | **0.44 s** | — |
| 729  | 12.1 s | **0.50 s** | 24× |
| 2197 | 36.4 s | **0.60 s** | **61×** |
| 4096 | 114.2 s | **0.78 s** | **146×** |

¹ off-run not repeated for this size; the original-binary baseline was 7.6 s.

The accelerated curve is the point: **0.37 → 0.78 s while the geometry grows
33×**.  The residual growth is mostly grid build time and the larger beam
field sampling more of the lattice.  Against the pre-change binary (49.4 s at
2197 zones) the end-to-end factor is ~80×; the same-binary A/B is the fairer
number, so both are quoted.

Regression check on a small geometry (`tests/cases/02_sobp`, 4 zones, 20 000
histories, 3 repeats): 7.57–7.78 s (on) vs 7.62–7.69 s (off) — within noise,
*after* adding the small-zone/small-geometry gates.  The first prototype
without gates was 38 % *slower* here, which is the classic acceleration-
structure lesson: an index must pay for its own bookkeeping, so always
benchmark the case it cannot help.

---

## 6. How this was approached (the part you asked to learn from)

1. **Profile first, and read the *shape*, not just the percentages.**  The
   callgrind tables in PR #139 said 88 % geometry at 1000 zones vs 19–31 % at
   4 zones.  A percentage that moves with problem size means the asymptotic
   complexity is wrong; that is worth more than any flat hotspot.
2. **Read the code until the invariant appears.**  The decisive facts were
   (a) all distance operators are `minpos` → distance is a flat min over
   leaves; (b) membership is a pure function of point and program → it can
   be partially evaluated per region of space; (c) the stepper re-checks
   membership → candidate sets may shrink as long as real boundaries stay.
   None of these are visible in the profile — only in the source.
3. **Build the A/B switch before optimizing.**  The `OSH_GEMCA_ACCEL=0`
   escape hatch turned every later question ("did this change physics?",
   "is the small case slower?") into a one-command experiment on a single
   binary.  Comparing across *builds* failed immediately (a longer
   `git describe` string from a newly fetched tag shifted every byte in the
   BDO) — same-binary A/B avoids a whole class of false alarms.
4. **Make the structure conservative, then prove it.**  Every cull is backed
   by "the box contains everything the epsilon logic can accept".  When in
   doubt (non-orthonormal transform, infinite plane), mark unbounded and
   don't cull — correctness degrades to the old behavior, never past it.
5. **Verify at multiple levels** (tests, bit-identity, argument), and
   **benchmark the unfavorable case** — that is how the missing gates were
   found.
6. **Stage the wins.**  AABB slab culling alone gave 4× (still O(N));
   the DDA walk + specialized programs took the same cases to 60–146× and
   flattened the curve.  Landing the first stage made debugging the second
   tractable, because outputs had to stay identical at every step.

---

## 7. What I would do next (ranked)

1. **Productionize this prototype** — see limitations in §8.
2. **Zone continuity in transport** (`src/transport/osh_transport_ion.c:158`):
   every wavefront iteration relocates *every* particle from zone 0 via
   `get_zone_ref_batch`, although a particle that did not cross a boundary is
   still in its zone, and one that did usually entered a known neighbor.
   Carrying the previous `zone_ref` and *verifying* it first (one specialized-
   program evaluation, O(1)) before falling back to the full lookup would
   remove most remaining `get_zone` work in *every* geometry, simple ones
   included.  Caveat: in geometries with overlapping zone definitions,
   verify-first can legitimately return a different (later-indexed) zone than
   first-match-wins; it should be gated on a "zones are disjoint" check or
   validated against the scan in debug builds.
3. **libm transcendentals (~26 % flat)**: the stopping-power / straggling /
   scattering chain calls `log/pow/exp/cbrt` per step per particle.  The
   algorithmic fix is tabulation: precompute per-material tables on a
   log-energy grid at setup (the infrastructure half-exists for `LOADDEDX`
   external tables) and interpolate in the step loop; the wavefront layout
   also admits 4-wide SIMD math as a follow-up.  Expected ceiling ~1.3× on
   physics-bound cases — worthwhile, but after geometry.
4. **Scoring raytrace (9–13 %)**: the mesh scorer already amortizes well;
   the main win is structural (per-thread accumulation + merge) and belongs
   to the parallelization track (#138).  The accelerator is read-only after
   build, so it composes cleanly with threading.
5. **RNG (2.5–4 %)**: leave it; stream-per-slot redesign is a threading
   precondition, not a speed win.

## 8. Prototype limitations (read before merging anything)

- The AVX2 *batch* zone lookup (`osh_gemca_runtime_get_zone_batch`) does not
  use the accelerator (it is not on the transport hot path today — transport
  goes through `get_zone_ref_batch` → scalar `get_zone`).  Unifying them
  would let the batch API benefit too, or alternatively the batch API could
  be retired.
- Gating thresholds (`OSH_GEMCA_ACCEL_DIST_MIN_LEAVES = 12`, min zones 16,
  ~8 cells/zone, ≤128 cells/axis) were set by quick measurement on two cases;
  a sweep over the C1–C8 matrix with the PR #139 harness should tune them.
- Memory is capped (8M specialized instructions) with graceful fallback, but
  a geometry with one huge zone overlapping every cell of a fine grid could
  hit the cap; a per-zone cell-count budget or hierarchical grid would be the
  fix.
- Voxel (DICOM) zones use the original path; CT-dominated cases see no
  change (they are already grid-traversed).
- `cgnode.bb_min/bb_max` in the cold structs are still unused TODOs; if the
  cold AST ever computes boxes, the runtime should consume rather than
  re-derive them.
- Unit tests for the accel internals (AABB derivation per body type, folding
  rules, DDA edge cases: rays along cell faces, zero direction components,
  points on the domain boundary) should be added before this leaves
  prototype status; current coverage is indirect (full suite + bit-identity
  + AVX2-vs-scalar cross-check).

## 9. Reproducing

```bash
cmake --preset release -DOSH_BUILD_EXAMPLES=OFF && cmake --build --preset release -j

# stress case from PR #139's generator (branch: pull/139/head)
python3 tools/bench/gen_gemca_stress.py --zones 1000 --out /tmp/c8

cd /tmp/c8
time build/bin/openshieldhit -n 300 .                    # accelerated
time OSH_GEMCA_ACCEL=0 build/bin/openshieldhit -n 300 .  # original algorithm
# outputs are byte-identical except the embedded timestamp
```
