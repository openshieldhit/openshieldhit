# GPU port plan, revision 2 (measured on Athena)

Status: living document. Supersedes the milestone ladder of `GPU_PORT_TASK.md`
(archived on the `glm/gpu-port-attempt` branch) and revises several decisions
of issue [#231] where on-hardware measurements contradicted or sharpened the
assumptions. Everything in the *Measured facts* section was produced on an
Athena A100 node on 2026-07-06, branch `fable/gpu-port-r2`; commands and full
context live in the commit messages of that branch.

Reading order for a newcomer: issue [#231] for the architectural survey and
rationale (still mostly valid), then this file for what is *proven*, what
changed, and what to do next.

---

## 1. Where the port stands

Done, merged on this branch, all with `ctest` green (126/126) and the c1
dose payload bit-identical to `main`:

| Piece | State |
|---|---|
| GEMCA flat RPN instruction store (`insns_flat[]`/`insn_begin[]`) | done + unit test (`test_gemca_flat_insns`) |
| `OSH_HD` host/device shim (`src/common/osh_hd.h`) + `OSH_RESTRICT` | done |
| Device-compilable slice: RNG (PCG32 + xoshiro256** + seeding/distributions), vector ops, material lookups, Bethe/Highland/Gauss-straggle, ray transform, **full GEMCA zone-membership evaluator chain** | done; entire surface compiled *and executed* as device code (`osh_gpu_hd_guard`) |
| Build system: C flags scoped per-language, `OSH_ENABLE_CUDA` option, `src/gpu/` module | done |
| GPU zone-membership kernel with **exact** CPU parity | done (`osh_gpu_zone_parity`) |
| Device RNG known-answer + throughput bench | done (`osh_gpu_rng_bench`) |
| FP64 atomic scoring-pattern bench | done (`osh_gpu_atomics_bench`) |

Build and run on an Athena GPU node:

```bash
module load CMake/3.31.8 GCC/13.3.0 CUDA/12.8.0
cmake -S . -B build_cuda -DCMAKE_BUILD_TYPE=Release -DOSH_ENABLE_CUDA=ON
cmake --build build_cuda --parallel 16
build_cuda/bin/osh_gpu_hd_guard
build_cuda/bin/osh_gpu_zone_parity examples/01_sdl_viewer/geo_RCC03.dat 110 1000000
build_cuda/bin/osh_gpu_rng_bench
build_cuda/bin/osh_gpu_atomics_bench
```

Verified toolchain: driver 570.195.03 (CUDA ≤ 12.8), nvcc 12.8.61 with GCC
13.3 or 14.3 host compiler, CMake 3.31.8, sm_80. Helios/GH200 remains
unverified assumption territory — first session there starts with [#231] §P0.

---

## 2. Measured facts (A100-SXM4-40GB, CUDA 12.8, this repo)

### 2.1 CPU reference (c1_p100_water_dose, 200k primaries, one EPYC 7742 core)

1212 ± 1 primaries/s over repeats (release, GCC 13.3; the archived
GCC 12.3 baseline was 1143), ≈ 0.97 M steps/s, 799 steps/primary. Phase
split: step 87 %, distance 11 %, zone_ref 1.7 %. Nuclear events and
secondaries in c1: **zero** (relevant for milestone G1 scoping below).

### 2.2 Zone membership (`osh_gpu_zone_parity`, 10⁶ rays/geometry)

| geometry | CPU scalar | CPU AVX2 batch | GPU kernel-only | GPU end-to-end* |
|---|---|---|---|---|
| geo_RCC03 (6 zones, offset spheres/rods) | 14.4 M/s | 24.0 M/s | **4047 M/s** | 143 M/s |
| 08_minimal_cyl (3 zones, DIFF chains) | 17.2 M/s | 46.0 M/s | 5137 M/s | 142 M/s |
| 01_sdl_viewer geo (3 zones) | 19.3 M/s | 51.7 M/s | 5523 M/s | 147 M/s |
| c1 water slab (3 zones) | 19.4 M/s | 53.5 M/s | 5524 M/s | 144 M/s |
| tilted RCC (BZALIGN rotation) | 15.5 M/s | 33.9 M/s | 4960 M/s | 146 M/s |

\* end-to-end = alloc + H2D of 48 MB SoA + kernel + D2H of 8 MB, per call.

**Parity is exact: 0 mismatches in 5×10⁶ rays**, across translated
(BCALIGN), rotated (BZALIGN), CSG-DIFF, and outside-all-zones cases.
The kernel executes the same `OSH_HD` source lines as the CPU
(`_osh_gemca_rt_find_zone_flat_hd`), so this is by construction, not luck.

Two headline consequences:

1. One A100 evaluates membership ~170× faster than one AVX2 CPU core on
   the kernel itself. Geometry will not be the GPU bottleneck.
2. The 28× gap between kernel-only and end-to-end is pure PCIe traffic.
   **The particle pool must live on the device across wavefront
   iterations**; any per-phase or per-batch host round-trip forfeits the
   win. Piecemeal offload (CPU transport calling GPU geometry batches) is
   an anti-goal — at ~143 M/s end-to-end vs ~50 M/s for AVX2 in-cache it
   could not even pay for its own complexity.

### 2.3 Project RNG on device (`osh_gpu_rng_bench`)

Per-thread `struct osh_rng` streams (48 B state, registers), 10⁶ threads:

| draw | throughput |
|---|---|
| PCG32 u32 | 392–500 G/s |
| xoshiro256** u64 | 292–372 G/s |
| uniform double | 245–313 G/s |
| gauss01 (Box–Muller + rejection) | 56–72 G/s |

Known-answer vs host, 8192 histories × 32 draws each, seeded via
`osh_rng_seed_history` semantics: **u32 and uniform-double sequences
bitwise identical**; gauss01 deviates by ≤ 7×10⁻¹² *relative* (glibc vs
CUDA `log()` last-ULP differences, amplified where `log(s) → 0`).

Consequences:

1. The per-history stream contract (`docs/dev/random_numbers.md`)
   transfers to the device **unchanged and bitwise**. A GPU run consumes
   the same per-history random sequences as a CPU run with the same seed.
2. Transport needs ~10³–10⁴ draws/history; even at 10⁹ histories/s the RNG
   budget is < 10 % of a single engine's measured rate. **Philox
   (#231 D3) is demoted from prerequisite to P4 experiment** — its real
   selling point on device would be shrinking the 48 B carried state, a
   register/spill question for `ncu`, not a throughput one.

### 2.4 FP64 atomic scoring patterns (`osh_gpu_atomics_bench`)

| bin pattern | atomicAdd rate |
|---|---|
| no-atomic reference | 750 G/s |
| 1 bin (total collision) | 0.45 G/s |
| 256 bins (1D Bragg depth-dose) | 4.7 G/s |
| 4096 bins | 15.7 G/s |
| 64k–1M bins (CT-like) | 62–69 G/s |

Transport at 100–1000× one CPU core deposits ~0.1–1 G/s. So: **atomicAdd
first (#231 D4 confirmed)**; CT-sized grids never notice; tiny 1D grids
lose at most tens of percent, and per-block shared-memory privatization is
the measured, small-page-only fallback. The `--score-replicas` merge
harness ([#230], now on main) is the CPU dress rehearsal for the replica
variant if determinism ([#168]) is wanted.

---

## 3. Post-mortem of the first attempt (what changed our rules)

The GLM-5.2 session (log archived as `session-ses_0d27.md` on
`glm/gpu-port-attempt`) did competent CPU-side prep work — the flat
instruction store and the HD-header pattern were cherry-picked onto this
branch nearly verbatim — and then failed at the first device milestone in
an instructive way.

**Failure 1: transcribed physics.** Its kernel *re-implemented* surface
evaluation from guessed semantics instead of sharing the CPU source.
Result: four silent divergences —

- the fatal one: it applied a full 4×4 transform to every body, but for
  `OSH_COORD_BCALIGN` bodies the geometry compiler populates **only the
  translation column** of `t[]` (see `osh_gemca2_calc_body.c:172`); the
  rotation block is zeros, so every query point collapsed to a constant
  and almost every ray reported "no zone";
- swapped `DIFF` operand order (right∖left instead of left∖right);
- ELLIPSOID/ELLZ/CONE silently treated as "outside" (CPU default differs);
- on-surface direction disambiguation dropped.

It then spent 2+ hours re-checking opcode tables and stack logic — the
parts it had gotten right — because its only oracle was the end-to-end
zone answer, which cannot localize a data-interpretation bug.

**Failure 2: "device-compilable" code that nvcc had never seen.** The
committed RNG HD headers called host-only `.c` exports from `OSH_HD`
bodies and used a bare C99 `restrict` in a public header. Plain-C builds
(where `OSH_HD` is empty) accept all of that; the first real nvcc build
of the slice fails. The property being claimed was never once exercised.

**Rules derived (binding for all further GPU work):**

R1. **No numeric or geometric logic in `.cu` files.** Kernels are launch
    shims: index math, `struct ray` assembly, one call into an `_hd`
    function. If a kernel needs logic that exists on the CPU, hoist the
    CPU body into an `_hd.h` header first and delegate the `.c` original.

R2. **Every `_hd` slice gets a device TU that instantiates it, in-tree,
    built by `OSH_ENABLE_CUDA=ON`, from the day the header exists.**
    `osh_gpu_rng_bench.cu` plays this role for RNG; the parity harness
    for geometry. A compile-only CI job (no GPU needed) makes it stick —
    see G4 below.

R3. **Parity harnesses precede kernels, and test at the finest available
    granularity.** For geometry that means per-body/per-transform twins
    would have localized the BCALIGN bug in minutes. Build the Level-0
    twin test first, then the kernel that must pass it.

R4. **When CPU and GPU disagree, diff data before logic**: dump the
    mirrored bytes back and compare against the host arrays, then compare
    intermediate values for one failing input. The uploaded-bytes check
    catches interpretation bugs (coord conventions, unpopulated fields)
    that logic-staring cannot.

R5. **Validation gates must match the arithmetic class.** Bitwise gates
    for integer/uniform-RNG and membership predicates (proven attainable);
    ULP-band gates (~10⁻¹¹ relative headroom) for anything through
    `log/exp/pow`; statistical gates ([#231] §8) for full transport.
    #231's blanket "≤ 1 ULP" for Level-0 would spuriously fail gauss01.

---

## 4. Revised decisions (deltas against issue #231)

| # | #231 said | Now |
|---|---|---|
| D1 single-source CUDA C++ backend | keep | **Confirmed by counterexample** (post-mortem above) and by 0-mismatch parity. Strengthened into rules R1–R2. |
| D2 wavefront phase kernels first | keep for the full engine | For the **first end-to-end milestone**, a run-to-completion history-per-thread kernel is simpler and sufficient (see G1): c1 has zero secondaries, so no pool interaction, no compaction, no phase barriers are exercised anyway. Wavefront machinery enters with secondaries (G2). |
| D3 Philox as third engine, early | prerequisite-ish | **Demoted to P4 experiment.** PCG32 per-history streams are bitwise-portable and ~400 G u32/s on device; carried-state size (48 B) is the only open question and it is a profiling question. |
| D4 atomicAdd first, privatize if needed | keep | **Confirmed with numbers** (§2.4). Privatization only ever for small 1D pages. |
| D5 batches = GPU launches | keep | Unchanged; additionally the pool must be device-resident within a batch (§2.2). |
| D6 statistical validation | keep | Refined into the three-class gate ladder (R5). |
| D7 pointer-free mirrors | keep | Implemented for GEMCA; `gemca_rt_body` carries a borrowed host `hu` pointer, so the mirror **refuses VOX geometries** until the CT grid gets its own mirror (G3). |
| D8 FP64 state, FP32 tables | keep | Untouched; A100 FP64 is fine (atomics bench ran entirely FP64 at 60+ G/s). |

New decisions:

| # | Decision |
|---|---|
| D9 | **Families are the GPU dispatch boundary.** The transport scheduler already drains per-family (ion → neutron → photon-stub). Each family gets its own device loop or a CPU fallback, selected per family. No cross-family megakernel. This is what keeps the growing physics roadmap (photons [seam exists, `osh_transport_photon.c` returns `OSH_ENOTSUP`], δ-electrons [STRAGG 3], INC cascade [#221–#225]) plannable: a new family or model lands CPU-first as today, becomes device-eligible when its step call tree is `OSH_HD`-clean, and until then that family drains on the host from the same banks. |
| D10 | **New hot-path physics must be born device-compilable.** A DEVELOPER.md amendment (G4) requires new step-callable code to follow the `_hd.h` pattern from its first PR — retrofitting is what costs sessions. The INC cascade (#221) and the SoA nuclear hot path (#226/#242/#243) should be written against this rule; tabulated σ_R (#242) is *also* the device-friendly choice (table lookup vs branchy parametric evaluation). |
| D11 | **Rare heavy branches become compacted follow-up kernels, not in-step branches.** Nuclear events are ~1 per 10³ steps; on device, service them by collecting flagged slots and running a dedicated kernel. This is the same shape #231 held in reserve for divergence, promoted to the default design for *heavy* rare physics (INC will be far heavier than Tripathi σ + FBU). |

---

## 5. Milestone ladder (revised)

Naming: G-milestones (this plan) to avoid clashing with #231's P-phases
and the archived task's M-milestones. Each is one PR-sized unit with its
gate stated. CPU behaviour must be bit-identical after every step (the
established identity check: c1 dose payload bytes).

### G1 — Tracer bullet: c1 depth-dose entirely on device

Scope fence: protons, CSDA + Highland MSCAT + Gaussian straggling, NUCRE
off (c1 measures zero nuclear events), analytic geometry, one scorer
(depth dose), single GPU.

1. Material-runtime mirror (tables are dense float arrays + scalars —
   the mirror is one struct copy + a handful of `cudaMemcpy`s).
2. Species table: the pool's `struct particle const **species` is a host
   pointer array (#231 G2); add the dense `species_idx[]` + flat property
   table. Coordinate the layout with #226 so nuclear work reuses it.
3. Device pool: the SoA arrays + one `struct osh_rng` per slot,
   `cudaMemcpy`-ed once per batch, resident until the batch drains.
4. Ion-step `_hd` slice completion: the step phase functions in
   `osh_transport_ion_step.c` (setup/vacuum/length/hinge/straggle/commit)
   hoisted per rules R1–R2, with the per-slot status-code seam (#231 G3)
   replacing hot-path diagnostics.
5. **Run-to-completion kernel**: one thread transports one history to
   absorption/escape, depositing via `atomicAdd` into the mirrored
   depth-dose page. No compaction, no phase barriers — c1 has no
   secondaries, and per-history RNG streams make order irrelevant.
6. Gates: Level-0 twins (membership: bitwise; material lookup: bitwise or
   1-ULP; step kinematics: ULP-band), then c1 CPU-vs-GPU with #231 §8
   Gate A/B (per-bin z, R80/peak). Post primaries/s vs the §2.1 baseline.

Expectation to beat before declaring success: ≥ 20× one EPYC core
(#231's P2 target). Given §2.2–2.4 headroom, a healthy result is
plausibly ≥ 100×; if it lands under 20×, profile before proceeding —
something structural is wrong.

### G2 — Full proton/ion physics: secondaries, NUCRE, wavefront

- Nuclear σ + channel dispatch on device (with #226; tabulated σ_R #242
  helps both targets), secondary append via atomic pool cursor
  (drop-independent RNG splits [#213] make this safe), neutron bank
  copy-back to the host drain loop first (D9).
- This is where the **wavefront phase-kernel architecture** (#231 D2)
  lands, with CUB compaction, because live-set divergence and pool
  interactions now exist. Measure run-to-completion vs wavefront on c4
  (NUCRE case) — keep the winner, keep both if within noise (the HD core
  is shared either way; the two drivers are thin).
- Vavilov/Landau straggling (κ-dispatch divergence is a named watch-item),
  Molière tables mirrored (#231 G7).
- Gates: c1–c4, c6–c7 benchmark cases + `tests/cases/` integration suite
  under §8 Gates A–C.

### G3 — Scoring/geometry completeness: CT, LET, filters, checkpoints

- VOX body mirror (HU volume + grid descriptor; lifts the D7 refusal),
  Jacobs traversal on device, per-voxel ρ weighting.
- Fluence/LET/Qeff two-pass pages, filters, differential axes, and the
  voxel-crossing scratch decision (#231 G4: fuse vs fixed buffer —
  measure).
- Checkpoint/batch/wall-budget integration (D5): GPU batches must honour
  `--dump-every`, SIGUSR1, and `--max-time` exactly like CPU ones.
- Gate: c5 (CT case) + full matrix.

### G4 — Guards and productization (can start in parallel with G2)

- **Compile-only CI job**: `OSH_ENABLE_CUDA=ON` build of every `_hd`
  guard TU for sm_80+sm_90 in an `nvidia/cuda:12.x-devel` container — no
  GPU needed, kills the R2 failure class at PR time.
- DEVELOPER.md device-code amendment (D10, `.cu` style, `_hd.h`
  conventions, banned constructs: bare `restrict` in public headers,
  host-only calls inside `OSH_HD`, heap pointers in mirrored structs).
- `--gpu[=index]` CLI with an honest unsupported-feature matrix
  (refuse politely per family/scorer until G2/G3 land).
- On-cluster statistical regression suite as an sbatch script
  (`benchmarks/gpu/`), since GitHub runners have no GPUs.

### G5 — Performance engineering and scale-out

As #231 P4/P5, unchanged in content: pool-capacity sweeps, `ncu`
warp-efficiency audit of the step kernel, launch-overhead/CUDA-Graphs
work on the wavefront drain tail, Philox-vs-PCG32 register study (D3),
FP32 experiments behind §8 gates, then in-process multi-GPU history
sharding (worker-context machinery from #160/#161) with
`osh_scoring_accumulator_merge` at checkpoints, 1→8 GPU scaling study.
GH200/Helios bring-up starts with #231 §P0 verification there.

---

## 6. Relationship to the physics roadmap

The roadmap is growing while the port happens; the plan is shaped so the
two do not fight:

- **Photon transport** ([photon seam], TODO): a stub family today. When
  implemented CPU-first, its step code should be `_hd` from day one
  (D10); the family scheduler seam (D9) then gives it a device loop as a
  bounded follow-up, not a rewrite.
- **δ-electron transport / STRAGG 3 Urban** (TODO): same story; also the
  reason scoring's deposit seam must stay the single choke point on
  device exactly as on CPU (#158 discipline).
- **INC cascade + coalescence + pre-equilibrium** (#221/#223/#224/#225):
  heavy, rare, per-event physics → the D11 compacted-follow-up-kernel
  shape. Its data (nucleon configurations) should be designed SoA with
  #226.
- **Thermal neutrons** (#246, merged) and the neutron family generally:
  stays a host drain in G1–G2 (D9), device candidate only after the ion
  family is proven — JEFF Tier-1 tables are small (35 nuclides), so the
  mirror is easy when its time comes.
- **SoA nuclear hot path** (#226/#242/#243): co-design the species/isotope
  table with G1 step 2, and prefer tabulated σ over parametric on both
  targets.

---

## 7. Numbers to carry in your head

| quantity | value |
|---|---|
| CPU c1 transport | ~1.21 k primaries/s·core, ~0.97 M steps/s·core, 799 steps/primary |
| GPU membership kernel | 4.0–5.5 G rays/s (≈ 170× AVX2 core) |
| GPU PCG32 | ~0.4–0.5 T u32/s; uniform double ~0.25–0.3 T/s; gauss01 ~60–70 G/s |
| FP64 atomicAdd | 0.45 G/s worst-case same-bin; 4.7 G/s @256 bins; 60+ G/s CT-like |
| PCIe round-trip cost | end-to-end membership 143 M rays/s vs 4 G/s kernel-only (28×) |
| A100 pool budget | ~165 B/slot → 10⁷ slots ≈ 1.7 GB of 40 GB (fine) |
| literature (gPMC/FRED class) | 10⁵–10⁶ protons/s per GPU → our G1 target ≥ 20×/core is conservative |

[#161]: https://github.com/openshieldhit/openshieldhit/issues/161
[#168]: https://github.com/openshieldhit/openshieldhit/issues/168
[#213]: https://github.com/openshieldhit/openshieldhit/issues/213
[#221]: https://github.com/openshieldhit/openshieldhit/issues/221
[#226]: https://github.com/openshieldhit/openshieldhit/issues/226
[#230]: https://github.com/openshieldhit/openshieldhit/issues/230
[#231]: https://github.com/openshieldhit/openshieldhit/issues/231
[#242]: https://github.com/openshieldhit/openshieldhit/issues/242
[#243]: https://github.com/openshieldhit/openshieldhit/issues/243
[#246]: https://github.com/openshieldhit/openshieldhit/issues/246
[photon seam]: ../../src/transport/osh_transport_photon.c
