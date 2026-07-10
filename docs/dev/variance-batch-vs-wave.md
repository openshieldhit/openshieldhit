# Study: making the variance batch size equal the particle-wave size

*Status: design study / analysis (no production behaviour changed). Companion to
[`scoring.md`](scoring.md) §4 (batch-means uncertainty) and issues #169 / #209 / #195 / #213.*

This note answers a concrete design question: **could the Monte-Carlo variance
"batch" be made one particle wave, so the standard-error accumulator is folded
once after every wave instead of once per fixed fraction of the run?** It reports
where the numbers land today, what the change would cost and buy, and how the
picture shifts as the nuclear/secondary models grow.

The short answer: **yes, it is possible with essentially no hot-path change** — it
is a cadence choice the existing checkpoint/merge machinery already supports — but
it swaps today's *fixed batch count* for a *fixed batch size*, and that trade is
only a win in the large-`nstat` regime. The recommendation is **not** to hard-wire
it as the default, but to understand it as one point on a well-defined axis.

---

## 1. Two quantities that are easy to conflate

| Term | What it is | Set by | Default |
|---|---|---|---|
| **Wave** (wavefront width) | Primaries injected per ion-pool fill — "primaries in flight". A cache/parallelism performance knob. | `ion_wavefront_width` = `min(pool_capacity, nstat)` | `OSH_TRANSPORT_POOL_CAPACITY` = **4096** (override with `-p`) |
| **Batch** (variance unit) | The primary range `[b, b+K)` transported between two family-complete checkpoints, folded into the master as one weighted batch-means observation. | `osh_checkpoint_next_batch_size()` | derived: `K = ceil(nstat / 10)` |

Today these are **independent**. The wave is a fill inside
`osh_transport_ion_run_range()`; the batch is a range handed to that function by
the outer loop in `run_master_batched()` (`src/transport/osh_transport.c`). One
batch of `K` primaries therefore contains `ceil(K / width)` waves, and the M2
(Welford) accumulator is folded **once per batch**, not once per wave.

The default count comes from `OSH_SCORING_VARIANCE_DEFAULT_BATCHES = 10`
(`src/scoring/runtime/osh_scoring_runtime.h`), applied in `osh_simulation_run()`
(`src/simulation/osh_simulation.c`): a variance-tracking run with no explicit
cadence and no `--score-replicas` split gets a count cadence of `ceil(nstat/10)`,
so it always folds **exactly 10 batches**.

## 2. Concrete numbers — "after how many primaries do we get a new batch?"

**Today (fixed *count* = 10):** a new batch every `ceil(nstat/10)` primaries.

**Batch = wave (fixed *size* = wave width = 4096):** a new batch every 4096
primaries (fewer for the trailing partial batch and when `nstat < 4096`).

| `nstat` | today: primaries/batch | today: #batches | batch=wave: primaries/batch | batch=wave: #batches |
|---:|---:|---:|---:|---:|
| 1 000 | 100 | 10 | 1 000 | **1 → no error bar** |
| 4 096 | 410 | 10 | 4 096 | **1 → no error bar** |
| 8 192 | 820 | 10 | 4 096 | 2 |
| 10 000 | 1 000 | 10 | 4 096 | 3 |
| 50 000 | 5 000 | 10 | 4 096 | 13 |
| 100 000 | 10 000 | 10 | 4 096 | 25 |
| 1 000 000 | 100 000 | 10 | 4 096 | 245 |
| 10 000 000 | 1 000 000 | 10 | 4 096 | 2 442 |
| 100 000 000 | 10 000 000 | 10 | 4 096 | 24 415 |

The crossover is at `nstat/10 = width`, i.e. **`nstat ≈ 40 960`**:

* **Below** it, today's batch (`nstat/10`) is *smaller* than a wave — the run
  already folds **more** often than once per wave. Batch=wave would make batches
  *bigger*, dropping to `< 2` batches for `nstat ≤ 4096` and losing the error bar
  entirely.
* **Above** it, today's batch spans *many* waves (25 waves/batch at `nstat = 1 M`),
  and batch=wave folds **far more** often — hundreds to tens of thousands of
  batches.

So the proposal is precisely: **replace "always 10 batches" with "always
4096-primary batches".** Fixed count ⇄ fixed size.

## 3. Is it feasible, and how invasive?

**Feasible, and cheap.** A wave boundary is a natural family-quiescence point: the
ion loop only refills when the pool is fully drained (`pool->n == 0`), and the
family scheduler drains neutrons/photons right after the ion range returns. Making
batch = wave only requires the outer batch loop to hand
`osh_transport_ion_run_range()` a range exactly one wave wide — which the existing
checkpoint machinery already does when the count cadence equals the width. The
per-batch private accumulator + `osh_scoring_accumulator_merge()` fold in
`run_master_batched()` is unchanged; **no per-step (§10) hot-path code is touched.**

The minimal implementation is one line in the derived-cadence block of
`osh_simulation_run()`: use `sim->transport_ctx.ion_wavefront_width` as the cadence
instead of `ceil(nstat / OSH_SCORING_VARIANCE_DEFAULT_BATCHES)`.

### Prototype (no code change)

Setting the checkpoint cadence to the wave width via the existing
`--dump-every-primaries` flag realises batch = wave exactly (each range is one
fill). Measured on `tests/cases/01_simple_detect` (1×1×800 mesh, C-12 @ 400 MeV/u,
NUCRE on), Energy + Fluence with `Variance On`:

* **Values are unchanged.** `nstat = 100 000`, default (10 batches) vs batch=wave
  (25 batches): max relative difference in every value column = **0** at the 12-digit
  text precision (`nstat = 10 000`: 1.5 × 10⁻¹³). Batching only reorders the
  floating-point summation — exactly the invariance §4/#213 promise.
* **Only the error estimate moves, and only in its precision.** At `nstat = 100 000`
  the per-bin standard-error ratio (batch=wave / default) had median ≈ 1.15 with a
  spread of ≈ 0.4 across bins — i.e. the two estimates agree within their own
  (substantial) estimation noise. At `nstat = 10 000` batch=wave has only 3 batches
  (2 d.o.f.) and reads systematically low/noisy (median ≈ 0.66). Batch count is a
  **bias-free knob on the *precision* of the error estimate**, not on its
  expectation or on the scored result.
* **Merge overhead was negligible here.** 100 000 primaries: 49.9 s (10 batches) vs
  49.7 s (25 batches) — the merge is `O(npages × nbins)` and this grid is tiny.

## 4. Why it is statistically permissible at all

Batch-means normally demands "large enough" batches because MCMC samples are
autocorrelated. **That constraint does not apply here.** Each history's RNG streams
are a pure function of its global index (`src/transport/osh_transport_ion.c`,
`osh_rng_split` lineage+ordinal keying), so histories are *exactly independent*.
For independent observations the Welford/Schubert–Gertz batch-means estimator is
unbiased for **any** batch size ≥ 1 — a batch of one wave, or even one history,
still yields an unbiased `SE² = M2 / ((B−1)·W)`. The only hard floor is **B ≥ 2**
(one batch = zero degrees of freedom = no error estimate; see scoring.md §4).

This is the crux: batch=wave is *statistically sound*; the question is purely
engineering trade-offs.

## 5. The trade-offs

**In favour of batch = wave (fixed size):**

* **Fine-grained live error.** The M2 accumulator updates after every 4096
  primaries, so a usable error bar exists after ~8192 primaries instead of after
  20 % of the run. Directly valuable for the periodic-dump/preview path (#193):
  each partial dump can carry a meaningful, up-to-date error.
* **DOF grows with the run.** Long runs get hundreds–thousands of batches → a
  much more *stable* error estimate than a fixed 9 d.o.f.
* **Bounded secondary memory (see §6).**

**Against (why it should not simply become the default):**

* **Small runs lose the error bar.** `nstat ≤ 4096` → 1 batch → no estimate. Today's
  fixed-10 guarantees an estimate for any `nstat ≥ 10`.
* **Couples a statistics knob to a performance knob.** The wave width *is* the
  pool-capacity/`-p` cache-tuning dial, explicitly documented as *not* a
  correctness or statistics parameter. Under batch=wave, a user shrinking the pool
  for cache reasons (`-p 1024`) would quietly quadruple the batch count and change
  the error-estimate character. That conflation is a design smell.
* **Merge/checkpoint cost grows with `nstat × nbins`.** Negligible for an 800-bin
  mesh, but a 512³ voxel grid (~1.3 × 10⁸ bins) with variance on would pay
  `B × nbins` element-merges — tens of thousands of full-grid passes at
  `nstat = 10⁸`. Still likely below transport cost, but no longer free.
* **More family-drain cycles.** Every batch boundary drains all secondary families;
  batch=wave multiplies the number of drain cycles (fixed per-cycle overhead) by
  `nstat/width / 10`.

A cleaner framing than "batch = wave" is therefore **"expose the batch size /
count as a policy"** and let batch=wave be one selectable point — the fixed-count
default stays for out-of-the-box runs, while a fixed-size cadence is available when
a live error bar or high DOF is wanted. The cadence plumbing
(`osh_checkpoint_policy`) already supports both (`batch` = explicit size,
`every_primaries` = count); only the *default derivation* would need a policy.

## 6. How this changes with future dev (nuclear models)

The nuclear/secondary models (inelastic + elastic, neutron transport, Fermi
break-up — `src/physics/nuclear/`, `src/transport/osh_transport_neutron.c`) are the
main reason the batch boundary matters, and they push the trade-off in both
directions:

1. **Family-completeness is mandatory, and batch=wave still honours it.** A fold is
   only unbiased once every secondary family the batch's primaries banked has been
   drained into scoring — dropping the neutron/fragment numerator *and* denominator
   biases ratios like `DLET = Σ(LET·dose)/Σdose` (issue #195), it does not merely
   add noise. Both the current scheme and batch=wave fold at a family-complete
   boundary, so batch=wave stays correct; it just makes that boundary the wave
   boundary. The wavefront-drain structure means the wave is *already* a clean
   ion-family quiescence point, so this composes well.

2. **Secondary-pool memory becomes a win.** The neutron pool is currently sized to
   `nstat` (`simulation_neutron_pool_capacity`) because a single final-only batch
   can bank neutrons from every primary. A batch = wave only ever banks *one wave's*
   secondaries before draining, so the neutron pool — and every future secondary
   family's pool (fragments, photons, e⁻) — could shrink from `O(nstat)` to
   `O(width × secondaries_per_primary)`. As more secondary families are added, this
   bounded, wave-local footprint is an increasingly attractive property of
   batch=wave.

3. **Secondary headroom / drops (#213) stay wave-local.** The ion pool already
   reserves one `width` of headroom for a wave's recoils/fragments so nothing is
   dropped in the default config. Batch=wave keeps that per-wave invariant intact.
   Richer nuclear models raise `secondaries_per_primary`, which may require growing
   the headroom *multiple* — orthogonal to the variance question, but it interacts
   with wave sizing (a larger wave banks proportionally more secondaries per drain).

4. **Rare high-variance events make small batches a noisier *estimator*.** Nuclear
   observables (neutron dose, fragment LET) are dominated by rare events. With
   4096-primary batches, many batches contain *zero* such events, so the
   between-batch M2 for those tallies is driven by which wave happened to contain an
   interaction. Still unbiased (independence, §4) — but the error *estimate* is
   heavy-tailed and noisier than with larger batches that average over more events.
   For nuclear-dominated tallies specifically, batch=wave buys a finer live cadence
   at the price of needing more total histories before the error bar itself is
   trustworthy. This is the strongest statistical caution, and it grows with the
   nuclear models rather than shrinking.

## 7. Recommendation

* Batch = wave is **feasible and cheap** to implement (a cadence choice, no hot-path
  change) and **statistically valid** (histories are independent).
* It is a **fixed batch *size*** scheme; today's default is a **fixed batch
  *count*** scheme. Neither dominates: fixed size wins for live error bars, DOF on
  long runs, and secondary-pool memory under nuclear models; fixed count wins for
  small runs, decoupling statistics from the cache knob, and bounded merge cost on
  huge grids.
* Prefer **making the batch cadence an explicit policy** over silently binding it to
  the wave width. Keep the fixed-count default; offer a fixed-size (wave-aligned)
  cadence for the live-preview / high-DOF use cases and for the memory benefits that
  arrive with the growing secondary-family set.
