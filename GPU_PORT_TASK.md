# GPU_PORT_TASK.md — autonomous session: first GPU milestones for OpenShieldHIT

You are an autonomous coding agent running on a **Cyfronet Athena GPU node** with **one NVIDIA
A100-SXM4-40GB** allocated to your Slurm job. This file is your complete task specification for a
single session. Where this file conflicts with your own ideas, **this file wins**. Where this file
conflicts with `DEVELOPER.md` or `CLAUDE.md`, **those win** — they are the repository law.

The mission: carry OpenShieldHIT as far as you can, in one session, toward transporting particles
on the GPU — in small, verified, committed steps. The milestones are ordered so that **every one of
them leaves the repository better and produces evidence**, even if you never reach the last one.
A partially completed ladder with honest logs is a success; a broken build with grand claims is a
failure.

---

## 0. Read this first (before touching anything)

Read, in this order, before your first edit:

1. `CLAUDE.md` — the agent operating brief, including the eight hard style rules.
2. `DEVELOPER.md` §1 (declarations), §2 (type/const spelling), §4 (comments), §5 (return
   conventions), §6 (naming), §9 (module layout), §10 (no hot-path allocation), §11 (portability),
   §15 (SoA/offload readiness).
3. `llms.txt` — project map.
4. `src/gemca/runtime/README.md` — contains a **written GPU migration recipe** you will implement
   in M1; treat it as the design document for that milestone.
5. `docs/dev/random_numbers.md` §5–§6 — the reproducibility model (per-history RNG streams). You
   must not break it; you also get to *rely* on it.
6. `src/transport/osh_transport_ion.c` and `src/transport/osh_transport_ion_step.c` — the loop you
   are porting. Read the whole files, not summaries of them.
7. `benchmarks/performance/README.md` — measurement methodology (binding for any number you report).

Background context (do not fetch, summarised here): the maintainers' GPU plan is GitHub issue
openshieldhit/openshieldhit#231. Its strategy, which you follow: **CPU stays canonical; add a CUDA
backend; share physics/geometry code between host and device via `OSH_HD`-annotated headers;
history-per-thread kernels; FP64 transport and scoring with `atomicAdd`; validation is statistical,
never bitwise.**

---

## 1. Non-negotiable session rules

These rules exist because a previous autonomous session on this codebase (the threads port) lost
~40% of its time and shipped avoidable defects in exactly these spots. Each rule is a countermeasure.

- **R1 — Ground before theorising.** If a hypothesis about code fails **twice**, stop reasoning
  from memory. Re-read the actual file (`grep`/open it) before forming a third hypothesis. When
  debugging, quote the offending lines into your log *from the file, not from recall*. The previous
  session spent 30+ minutes asserting code existed that was never written; `grep` found it in 14
  seconds.
- **R2 — House style is not optional.** The most common violation: C99 in-loop declarations.
  `for (size_t i = 0u; ...)` is **banned** (§1: K&R, all declarations at top of block). Before every
  commit run `grep -rn 'for (size_t\|for (int\|for (unsigned' ` on your changed files — the count
  must be zero — and run `./tools/clang-format-all.sh` (if clang-format ≥ 18 is not available via
  modules, note that in the log and match surrounding formatting by hand; never reformat files you
  did not otherwise change). Also: no `typedef struct`, no `//` comments in committed code,
  `double const *p` not `const double *p`, predicates return `int`, operations return
  `enum osh_status`.
- **R3 — No number without its command.** Any measurement you mention anywhere (log, commit
  message, comment) must have been produced **in this session**, and the exact command plus raw
  output must be pasted into `GPU_PORT_PROGRESS.md` first. If two measurements disagree, resolve or
  report the discrepancy — never ship the stale number. Commit messages describe *what changed*;
  numbers live in the progress log.
- **R4 — Keep the tree clean.** All simulation runs happen **outside the repository**:
  `export RUNDIR=${SCRATCH:-$HOME}/osh_gpu_runs && mkdir -p "$RUNDIR"` and always pass
  `--outdir "$RUNDIR/..."` (and copy case directories there before editing them). `git status`
  before each commit must show only files you intend to commit. No `.bdo`, no profiles, no editor
  droppings.
- **R5 — Every CUDA call is checked.** Define one `OSH_CUDA_CHECK(call)` macro (Appendix B) and
  wrap **every** runtime call and every kernel launch (launch check + `cudaGetLastError`, plus a
  `cudaDeviceSynchronize` check in debug builds). Error paths must be defined behaviour: never use
  a resource whose creation failed (the previous session joined a thread that was never created).
- **R6 — Never weaken a test to make it pass.** Existing assertions are immovable. New statistical
  tolerances must be *derived* (the two-seed calibration in §5), not tuned until green. If a gate
  fails, the result is "gate failed, here is the evidence", not a looser gate.
- **R7 — Numeric scripting in Python only** (`python3`), with `LC_ALL=C` exported in your shell.
  No `bc`, no locale-dependent `printf` parsing — decimal-comma locales have broken this before.
- **R8 — Timebox and move on.** Per-milestone timeboxes are listed below. If you exceed one while
  blocked: write the blocker into the log (symptoms, what you tried, best hypothesis), commit what
  is safe to commit (a WIP commit on your branch is fine if it builds; say WIP in the title), and
  move to the next milestone that does not depend on the blocked one. Never delete failing work.
- **R9 — Do not touch:** physics constants/formulas, input parsers, public headers'
  API/ABI semantics, output file formats, or the behaviour of the serial CPU path. After every
  CPU-side refactor (M1, M2), re-run the identity check in §5.1 — CPU output must be
  **byte-identical** before/after your change.
- **R10 — Uncertainty bookkeeping.** If you merge accumulators, follow `docs/dev/scoring.md` §4
  exactly: `weight`/`nbatch` semantics are part of correctness (a prior session silently made
  batch count scale with worker count). For this session's GPU path, one GPU run of one history
  range = **one batch**.
- **R11 — Measure only release builds** (`cmake --preset release`); debug builds are for tests and
  `compute-sanitizer` only. Both must build at every commit.
- **R12 — Git discipline.** First action in the repo:
  `git checkout -b glm/gpu-port-attempt` (branch off the branch this file is on). Commit at every
  green gate, message style `module: summary` (e.g. `gemca: add flat RPN instruction store`). Do
  **not** push, do not open PRs, do not rebase or amend published history — the human collects the
  branch afterwards. End the session with everything committed, including the logs.

---

## 2. Deliverable files (create these early, append as you go)

- `GPU_PORT_PROGRESS.md` — append-only session log. For each action: timestamp, what/why (one
  line), the command(s), the raw output that matters (trimmed to the relevant lines), and the
  conclusion. This file is the primary artifact the humans will read. Honesty here outranks
  progress: a documented dead end is valuable, an undocumented success is not.
- `GPU_PORT_RESULTS.md` — written at the end (M5): environment summary, what was achieved per
  milestone, the validation evidence, the performance table, known defects/simplifications, and a
  prioritised "next session should…" list.

Both live in the repository root on your branch and are committed with regular commits.

---

## 3. M0 — Environment bring-up and CPU baseline (timebox ~45 min)

Nothing may be edited before M0 is green.

1. **Hardware/toolchain inventory** (paste all outputs into the log):
   ```bash
   hostname; uname -m                      # expect x86_64 on Athena
   nvidia-smi                              # expect 1× A100-SXM4-40GB visible
   module avail 2>&1 | grep -i -E 'cuda|cmake|gcc' | head -30   # or: module spider CUDA
   module load CUDA CMake GCC              # adjust to the names you actually found
   nvcc --version; cmake --version; gcc --version
   ```
   If no CUDA module exists, look for `/usr/local/cuda*/bin/nvcc`. If `nvcc` is genuinely
   unavailable, record it, skip M3–M6, and do M1/M2 only — they are pure CPU work.
2. **Build both presets** (out of the box, no options):
   ```bash
   cmake --preset debug   && cmake --build --preset debug   --parallel
   cmake --preset release && cmake --build --preset release --parallel
   ```
3. **Test suite** (note: there are **no CTest presets** — point at the binary dir):
   ```bash
   ctest --test-dir build_debug --output-on-failure
   ```
   Must be green before any edit. Record the pass count.
4. **CPU baseline** on the reference case (release binary, run in `$RUNDIR`):
   ```bash
   build/bin/openshieldhit -n 200000 --profile "$RUNDIR/base_prof.json" \
       --outdir "$RUNDIR/base_out" benchmarks/performance/cases/c1_p100_water_dose
   ```
   Record wall time, `primaries/s` from the run log, and the profile JSON phase breakdown.
5. Create `GPU_PORT_PROGRESS.md` with the inventory; first commit.

**Gate:** builds green, ctest green, baseline recorded, first commit made.

---## 4. The milestone ladder

### M1 — GEMCA flat instruction store (CPU-only; timebox ~60 min)

The geometry runtime is GPU-ready except for one documented gap: each zone's RPN instruction list
sits behind a per-zone heap pointer (`struct gemca_rt_zone.insns`), which a device kernel cannot
follow. `src/gemca/runtime/README.md` §"GPU migration path" specifies the exact fix: add
`insns_flat[]` (all zones' instructions concatenated) and `insn_begin[]` (per-zone offsets) to
`struct osh_gemca_runtime`, filled in `setup_zones`, **additively** — the CPU path keeps using
`zones[j].insns`.

- Implement exactly the README recipe (fields, population, teardown).
- Add a unit test `tests/unit/test_gemca_flat_insns.c` (auto-discovered by glob; see
  `tests/unit/README.md`): for every fixture geometry, assert `insn_begin[]` is monotonically
  increasing, totals match `ninsns` sums, and the flat array's contents equal the per-zone arrays
  element-by-element.
- §5.1 CPU identity check. Format. Commit (`gemca: add flat RPN instruction store for device use`).

**Gate:** new test passes, all old tests pass, CPU output byte-identical.

### M2 — `OSH_HD`: make the M4 call tree device-compilable (CPU-only; timebox ~90 min)

Introduce a tiny header (suggested `src/common/osh_hd.h`):

```c
#ifndef OSH_HD_H
#define OSH_HD_H
#if defined(__CUDACC__)
#define OSH_HD __host__ __device__
#else
#define OSH_HD
#endif
#endif /* OSH_HD_H */
```

Then make the *minimum* set of functions needed by the M4 kernel compilable for device. Mechanics:
functions that are already `static inline` in headers just gain `OSH_HD`; functions living in `.c`
files move their **bodies** to a new `*_hd.h` header (marked `OSH_HD static inline`) that the
original `.c` includes and re-exports with its unchanged public signature — so the CPU library's
API, layout, and behaviour do not move. No global mutable state may be reached from any `OSH_HD`
function; standard `<math.h>` calls are fine (CUDA supplies device overloads).

Minimum set (verify against the real call tree as you read it, R1):

| Area | Functions | Where today |
|---|---|---|
| RNG | PCG32 init/next, `osh_rng_u32/u64/double/float`, `osh_rng_gauss01/gauss`, `osh_rng_seed_history` + the `rng_mix_stream` mixer | `src/random/osh_rng.c`, `osh_rng_pcg32.c` |
| Material | `osh_material_runtime_sp_lookup`, `..._range_lookup`, `..._get_rho` | already `static inline` in `src/material/runtime/osh_material_runtime.h` |
| Geometry | scalar zone-membership evaluation (RPN walk, via M1's flat store) and current-medium boundary distance for the analytic bodies used by c1 (RCC etc.) | `src/gemca/runtime/osh_gemca_runtime.c`, `src/gemca/osh_gemca2_dist.c`, `..._calc_*.c` |
| Physics | effective charge / stopping helpers used by the step, Highland `theta0`, Gaussian straggling sigma | `src/physics/atomic/` (`osh_physics_bethe.c`, `osh_physics_scat_highland.c`, `osh_physics_strag_gauss.c`) |
| Kinematics | anything `osh_transport_ion_step.c` calls for the restricted config | `src/particle/`, `src/common/` |

- The particle pool's per-slot `struct osh_rng` is a plain 56-byte value type — it uploads as-is.
- Do **not** attempt to `OSH_HD` the whole of `osh_transport_ion_step.c` or the scoring module in
  this session; M4 re-expresses the restricted step on device instead (see M4 scope note).
- §5.1 CPU identity check after each sub-refactor; ctest green; format; one commit per module
  touched (`random: make PCG32 engine device-compilable`, etc.).

**Gate:** all tests green, CPU output byte-identical, and a scratch check that
`gcc -fsyntax-only` still accepts every touched header in C11 mode.

### M3 — CUDA skeleton + geometry parity kernel (first device code; timebox ~90 min)

1. **Build system:** root `CMakeLists.txt` gains
   `option(OSH_ENABLE_CUDA "Build the experimental CUDA backend" OFF)`. When ON:
   `enable_language(CUDA)`, `find_package(CUDAToolkit REQUIRED)`, `set(CMAKE_CUDA_ARCHITECTURES 80)`
   (A100), and add subdirectory `src/gpu/`. OFF must leave every existing target byte-identical.
   Configure with: `cmake --preset release -DOSH_ENABLE_CUDA=ON` (presets accept extra `-D`).
2. **Module `src/gpu/`** (respect §9: it may depend on runtime layers; nothing depends on it):
   `osh_gpu.h` (availability query + status enum), `osh_gpu_mirror.cu` (upload the gemca runtime's
   flat arrays — surfaces/bodies/zones/insns_flat/insn_begin — and the material runtime tables into
   device buffers; a small POD "device view" struct of device pointers + counts, passed to kernels
   **by value**), and `osh_gpu_zone_test.cu`.
3. **Parity kernel:** one thread per query point, evaluating zone membership with the same shared
   `OSH_HD` code the CPU uses. Host test program (suggested `examples/` or a `tests/unit` file
   guarded by the CMake option): generate ~10⁵ deterministic pseudo-random points spanning the c1
   geometry (inside water, vacuum gap, blackhole shell, outside), compute zone ids on CPU via
   `osh_gemca_runtime_get_zone_ref_batch` and on GPU, assert **exact** agreement; same for boundary
   distances within 1 ULP-ish (`rel err < 1e-12`).
4. Run it under `compute-sanitizer --tool memcheck` once (debug build). Log; commit
   (`gpu: CUDA skeleton + zone/distance parity kernel`).

**Gate:** parity test passes on the A100; `OSH_ENABLE_CUDA=OFF` build still green everywhere.

### M4 — Tracer bullet: c1 depth-dose entirely on the GPU (the core milestone; timebox ~3 h)

**Scope fence (deliberate, document every simplification in the log):** protons, analytic c1
geometry, `STRAGG 1` (Gaussian), `NUCRE 0` (no secondaries!), and **`MSCAT 1` (Gaussian) instead of
c1's default `MSCAT 2` (Molière)** — Molière needs host-built tables uploaded; that is a stretch
goal, not this milestone. Because CPU is always compared against CPU-with-identical-settings, this
is legitimate: make a **variant case** in `$RUNDIR` (copy `c1_p100_water_dose/`, set `MSCAT 1`
in `beam.dat`; in `detect.dat` set `Filename depth_dose.txt` and add `FileFormat TEXT` inside the
`Output` block for easy numeric comparison).

With `NUCRE 0` there are no secondaries, so **no pool compaction is needed**: use a
**one-thread-one-history megakernel** — each thread transports one primary from birth to death:

```
thread i: seed rng from global history index (osh_rng_seed_history, PHYSICS purpose — identical to CPU)
          read start state from uploaded pool arrays (host-filled)
          loop: zone lookup → boundary distance → step length (DELTAE/DEMIN rules)
                → energy loss + Gaussian straggling → Gaussian MCS deflection
                → deposit → advance; exit on E ≤ cutoff, blackhole, or geometry escape
          write per-thread termination status into status[i]
```

Implementation notes:

- **Reuse the host beam sampler:** fill the SoA pool on the host with
  `osh_beam_runtime_fill_pool_at()` (explicit global base — it exists for exactly this), upload the
  arrays (`x,y,z,ux,uy,uz,e,wt` + the `osh_rng` array), launch in chunks of ≤ 10⁶ histories to
  bound memory. Study how `run_history_range()` in `osh_transport_ion.c` wires seeding first.
- **Step logic:** re-express the restricted subset of `osh_transport_ion_step.c` phases
  (`ion_step_length`, energy/straggling, scatter, commit) as an `OSH_HD` device function, reading
  the CPU file phase by phase (R1: quote the CPU lines you mirror into your log). Vacuum zones
  advance to boundary without energy loss; blackhole kills.
- **Scoring:** c1 scores dose on a 1-D z-mesh (200 bins over z∈[0,10] cm). Deposit exactly like
  the CPU's crossing-split: distribute each step's ΔE over the z-bins it crosses **proportionally
  to path length in each bin**, `atomicAdd` on a device `double dose[200]` (native FP64 atomics on
  sm_80). Download at the end, feed the sums into the existing accumulator so postprocess/save run
  unchanged — or, minimally, write your own TEXT output *and* also run the CPU writer on the
  downloaded sums; prefer reusing the CPU save path.
- **Wiring:** smallest clean entry point wins. Suggested: an env-var switch `OSH_GPU=1` checked in
  `osh_simulation_run()`'s ion-family dispatch (document it), or a `--gpu` CLI flag if the parser
  edit stays small (`src/cli/osh_cli.c`). Unsupported configuration (any physics switch outside the
  scope fence) must fall back to CPU with a clear log line, never silently compute wrong physics.
- Run once under `compute-sanitizer --tool memcheck` and once with `--tool racecheck` (debug),
  small nstat. Then the validation protocol of §5.2.

**Gate:** §5.2 passes; both sanitizer runs clean; commit
(`gpu: one-thread-per-history CSDA megakernel for restricted c1 config`).

### M5 — Honest measurement + profile (timebox ~60 min)

- Throughput table (release, `-n` sized so each run ≥ 20 s; 3 repeats, report median + spread):
  CPU 1-thread vs GPU, in primaries/s, on the M4 variant case. Report the speedup **whatever it
  is** — a slow-but-correct GPU path passes; a fast-but-unvalidated one fails.
- `nsys profile -o "$RUNDIR/osh_m5" build/bin/openshieldhit ...` — record kernel time vs memcpy vs
  host phases; state whether the run is kernel-bound, transfer-bound, or launch-bound.
- One `ncu --set full` (or `--set basic` if `full` exceeds permissions/time) capture of the
  megakernel: note achieved occupancy, registers/thread, warp execution efficiency, DRAM
  throughput. If `ncu` hits `ERR_NVGPUCTRPERM`, record that and skip — do not fight permissions.
- A pool-chunk-size sweep (e.g. 10⁴/10⁵/10⁶ histories per launch) — one table.
- Write `GPU_PORT_RESULTS.md` (see §2). Commit.

**Gate:** results file complete with every number traceable to a logged command.

### M6 — Stretch goals, only with time left (any order, each is optional)

1. Molière MCS on device: upload the tables built by `osh_physics_moliere_init()`; validate c1
   with its original `MSCAT 2` against CPU.
2. Philox4x32-10 as `OSH_RNG_TYPE_PHILOX4X32` (CPU reference + known-answer tests per the tutorial
   in `docs/dev/random_numbers.md` §8) — the designated GPU engine going forward.
3. Wavefront-kernel refactor (separate zone/distance/step kernels + CUB compaction) — the
   architecture the full port needs once secondaries (NUCRE 1) enter; do not start unless ≥ 2 h
   remain.
4. `benchmarks/`: extend `run_bench.py` with a `--gpu` axis and GPU metadata in the machine block.

---

## 5. Verification protocols

### 5.1 CPU identity check (after every CPU-side refactor — M1, M2)

```bash
build/bin/openshieldhit -n 20000 --outdir "$RUNDIR/ident_a" benchmarks/performance/cases/c1_p100_water_dose
# rebuild with your change, then:
build/bin/openshieldhit -n 20000 --outdir "$RUNDIR/ident_b" benchmarks/performance/cases/c1_p100_water_dose
python3 - "$RUNDIR/ident_a/depth_dose.bdo" "$RUNDIR/ident_b/depth_dose.bdo" <<'EOF'
import sys
a, b = (open(p, 'rb').read() for p in sys.argv[1:3])
print('identical' if a == b else 'DIFFER'); sys.exit(a != b)
EOF
```

Byte-identical is the requirement (run both *after* your change from the same binary is a mistake —
run A from the pre-change build you keep in `build_ref/`, or check out the previous commit; note
BDO may embed paths/timestamps — if the naive compare differs on a **no-op rebuild**, fall back to
comparing the ASCII TEXT output of a variant case instead, and record that caveat once).

### 5.2 GPU statistical equivalence (M4 gate) — the two-seed calibration

Bitwise CPU/GPU agreement is impossible (FMA/rounding/atomic order); "close enough" must be
*derived*, not eyeballed:

1. `caseA` = the M4 variant case. `caseB` = copy of `caseA` with a different `RNDSEED` in
   `beam.dat` (e.g. `sed -i 's/^RNDSEED.*/RNDSEED  13579111/' beam.dat`). Do **not** use
   `-N/--seedoffset` for this — offset ranges can overlap history indices between runs.
2. Three runs, same `nstat` (≥ 5×10⁵; scale so the CPU runs stay ≤ 10 min):
   CPU(caseA), CPU(caseB), GPU(caseA).
3. Parse the three TEXT outputs (columns per `docs/user/detect.dat.md`; last column is the value).
   Over bins with value > 1% of the peak, compute per-bin relative differences for
   D_cal = CPU(A) vs CPU(B) — pure statistics — and D_test = CPU(A) vs GPU(A).
4. **Pass iff:** mean(|D_test|) ≤ 1.5 × mean(|D_cal|) **and** max(|D_test|) ≤ 2 × max(|D_cal|)
   **and** the R80 range depth (distal 80% crossing) of GPU is within the CPU(A)–CPU(B) spread
   (or 0.05 cm, whichever is larger) **and** total integrated dose agrees within
   1.5 × the A–B total-dose difference.
5. Paste the comparison script and its output into the log. A systematic shape difference (e.g.
   every bin beyond the peak biased one way) fails the gate even if the means squeak through —
   plot or print the residual-vs-depth pattern and look at it.

---

## 6. Explicit non-goals for this session

Nuclear reactions / secondaries on device; neutron transport; CT/voxel geometries; DICOM;
LET/fluence/differential scorers; multi-GPU; MPI; Windows/MSVC + CUDA; performance tuning beyond
M5's measurements; touching `include/openshieldhit/` public APIs; creating pull requests; pushing.

---

## Appendix A — style traps that have burned agents here (self-check before each commit)

| Trap | Rule |
|---|---|
| `for (size_t i = 0; ...)` | §1 — declare `i` at top of block; K&R |
| `typedef struct {...} foo;` | §2.1 — `struct foo` spelled out everywhere |
| `// comment` in committed code | §4.1 — block `/* */` only |
| `const double *p` | §2.2 — `double const *p` |
| `strcasecmp`, `<unistd.h>`, `<sys/stat.h>`, `<threads.h>`, `mkdtemp`, `getpid` | §11.1 — banned; Windows is a target (host code; device code is Linux-only but keep it clean) |
| Returning 0/1 from an operation, or `enum osh_status` from a predicate | §5 — predicates `int` (1 = yes); operations `enum osh_status` (`OSH_OK` = 0) |
| `malloc`/`free` anywhere under the transport/scoring step | §10 — pre-allocate at setup; device buffers allocated once at mirror time |
| Public function without `osh_` prefix / non-`static` file-local helper | §6.1 |

## Appendix B — CUDA hygiene skeleton

```c
/* src/gpu/osh_gpu_check.h */
#ifndef OSH_GPU_CHECK_H
#define OSH_GPU_CHECK_H
#include <stdio.h>
#include <cuda_runtime.h>

#define OSH_CUDA_CHECK(call)                                                        \
    do {                                                                            \
        cudaError_t err_ = (call);                                                  \
        if (err_ != cudaSuccess) {                                                  \
            fprintf(stderr, "CUDA error %s at %s:%d: %s\n", cudaGetErrorName(err_), \
                    __FILE__, __LINE__, cudaGetErrorString(err_));                  \
            return OSH_ESTATE;                                                      \
        }                                                                           \
    } while (0)
#endif /* OSH_GPU_CHECK_H */
```

Kernel launches: `kernel<<<grid, block>>>(...); OSH_CUDA_CHECK(cudaGetLastError());` and in debug
builds additionally `OSH_CUDA_CHECK(cudaDeviceSynchronize());`. Device views are POD structs of
device pointers + counts passed by value; never pass a struct containing *host* pointers to a
kernel.

## Appendix C — key file map (verify by reading; do not trust from memory)

| What | Where |
|---|---|
| Wavefront loop, seeding, phase timers | `src/transport/osh_transport_ion.c` |
| Per-step physics phases | `src/transport/osh_transport_ion_step.c` |
| SoA pool (+ per-slot RNG) | `src/common/osh_particle_pool.h` |
| Geometry runtime + GPU recipe | `src/gemca/runtime/` (`README.md`!) |
| Material tables + inline lookups | `src/material/runtime/osh_material_runtime.h` |
| RNG engines, per-history seeding | `src/random/osh_rng.c`, `docs/dev/random_numbers.md` |
| Scoring deposit seam & accumulators | `src/scoring/runtime/osh_scoring_step.c`, `osh_scoring_accumulator.h` |
| Simulation driver (dispatch point for `--gpu`) | `src/simulation/osh_simulation.c` |
| CLI parser | `src/cli/osh_cli.c` |
| Reference case | `benchmarks/performance/cases/c1_p100_water_dose/` |
| Unit-test conventions | `tests/unit/README.md` |

Good luck. Small steps, real evidence, honest logs.
