# Parallelization Plan

This document analyses where OpenShieldHIT stands with respect to parallel
execution, surveys the strategies used by other Monte Carlo transport codes,
and lays out a phased roadmap from the current single-threaded wavefront loop
to a dispatcher that drives SIMD lanes, CPU cores, GPUs, and HPC nodes from
one work-decomposition model.

## 1. Where the code stands today

The codebase is unusually well prepared for parallelization — several key
design decisions were made with SIMD and GPU execution in mind:

**Already in place:**

- **Wavefront (BFS) transport loop** — `src/transport/osh_transport_ion.c`
  transports up to `OSH_TRANSPORT_POOL_CAPACITY` (4096) histories
  simultaneously: fill pool → batch zone/distance query → step every live
  slot → compact → repeat.  This is exactly the *event-based / track-slot*
  execution model that GPU transport codes (Celeritas, AdePT, Shift) converged
  on.  The capacity constant is documented as the cache-vs-parallelism knob,
  with `capacity == NSTAT` called out as "natural for GPU offload".
- **SoA particle pool** — `src/common/osh_particle_pool.h` stores phase space
  as flat per-field arrays in a single slab, explicitly so kernels can be
  auto-vectorized without gather/scatter.  The header already anticipates
  "thread-local pools, with atomic scoring tallies".
- **Streamed RNG** — `src/random/osh_rng.h` (PCG32, xoshiro256**) has a
  `stream` parameter for independent sequences, separate beam/physics streams,
  and vector draw APIs (`osh_rng_*_vec`).
- **SIMD runtime-dispatch precedent** —
  `src/gemca/runtime/osh_gemca_runtime_avx2.c` is compiled with
  `-mavx2 -mfma` when available and selected at runtime via
  `__builtin_cpu_supports("avx2")`, with a scalar fallback.
- **No heap allocation on the hot path** (enforced convention), and a
  compile/runtime split that keeps physics tables, geometry, and material
  data read-only during transport — they can be shared by all workers.
- **Transport-family scheduler seam** — `src/transport/osh_transport.c` and
  `src/transport/README.md` already define the orchestration shape: a
  lightweight scheduler on top, family kernels submitting chunks to a shared
  worker pool, thread-local secondary buffers merged at explicit sync points.

**Missing — the actual work items:**

1. **Per-history RNG semantics.**  Today two RNGs (beam, physics) are
   consumed sequentially across the whole run, so results depend on pool
   capacity, execution order, and any future thread count.  This is the
   single biggest blocker for every phase below.
2. **Mergeable scoring.**  `osh_scoring_score_step()` accumulates into shared
   runtime arrays with no partial-result or merge concept (TODO.md already
   lists "result merge API for embarrassingly parallel runs").
3. **A threading layer.**  Nothing concurrent exists yet; DEVELOPER.md bans
   `<threads.h>` (MSVC) and prescribes pthreads/Win32 or OpenMP.
4. **Batched inner-loop phases.**  The per-slot step
   (`osh_transport_ion_step()`) is scalar; TODO.md lists batch stopping-power
   lookup and the AVX2 `get_distance_batch` path as open items.
5. **Any GPU backend.**

## 2. What other Monte Carlo codes do

| Code | Parallel model | Key techniques |
|------|---------------|----------------|
| **Geant4** | Event-level multithreading (since 10.0), task-based with TBB (since 10.7), MPI on top | Master/worker; geometry + physics tables shared read-only; per-worker RNG seeded per event; thread-local scoring merged in reduction step |
| **Celeritas** (ORNL/Fermilab) | GPU (CUDA/HIP), event/track-slot based | Fixed array of track slots in SoA, loop over per-action kernels with masking, optional track sorting to reduce divergence, per-slot RNG state, atomic tally accumulation; integrates as a Geant4 offload library |
| **AdePT** (CERN) | GPU EM transport | Same track-slot shape; per-particle-type queues (e−/e+/γ) to reduce divergence — analogous to OpenShieldHIT's transport families |
| **OpenMC** | OpenMP + MPI history-based; event-based mode for GPU via OpenMP target offload | Per-history RNG streams → results independent of thread count; thread-local tallies + reduction; event-based mode ran at exascale (Frontier/Aurora) |
| **MCNP, PHITS** | Hybrid MPI + OpenMP, history-level | Tally replication per thread/rank, reduction at end |
| **FLUKA** | Independent processes, offline merge | The "zero-risk" baseline: N jobs with different seeds, statistical merge |
| **MCsquare** (UCLouvain, CPU) | OpenMP + hand-written SSE/AVX SIMD | SoA proton batches, tabulated physics, vectorized step phases — the closest published model to OpenShieldHIT's intended CPU path; ~minutes for clinical plans on a workstation |
| **FRED** (Rome/Kraków) | GPU (CUDA/OpenCL) ion transport | Clinical-plan recalculation in minutes on a single GPU; voxel geometry; demonstrates the end-state for the CT use case |
| **gPMC, GPUMCD, MOQUI, RayStation dose engine** | GPU, history-per-thread or event-based | MOQUI uses key-value (hash) scoring to cut tally memory; vendor engines confirm GPU MC is the clinical production standard |

Distilled lessons:

1. **History-level decomposition is universal.**  Every code, from FLUKA jobs
   to Frontier-scale OpenMC, ultimately splits work by primary history.
   Everything else (SIMD, GPU) is an optimization *inside* that decomposition.
2. **Reproducibility comes from per-history RNG streams**, never per-thread
   ones.  Counter/stream-based generators (PCG streams, Philox) let any
   worker compute history *i* identically regardless of scheduling.
3. **Scoring scales by replication + reduction** on CPU (thread-local
   tallies) and by **atomics** on GPU (`atomicAdd` on dose grids is cheap on
   modern hardware; hash-based scoring if memory becomes the limit).
4. **GPU codes converge on exactly the wavefront-pool shape** OpenShieldHIT
   already has: fixed SoA track slots, kernel loop, compaction, family queues
   to limit divergence.
5. **Shared read-only data** (geometry, material tables) must be immutable
   during transport — already guaranteed by the compile/runtime split.

## 3. Target architecture: chunk dispatcher

One abstraction serves every tier:

```
            NSTAT primaries
                 │
        ┌────────┴─────────┐
        │  chunk dispatcher │   chunk = (first_history, n_histories)
        └────────┬─────────┘
   ┌─────────┬───┴─────┬──────────────┐
   serial    CPU pool   GPU device(s)  MPI ranks
   backend   (M threads,  (1 huge       (file- or
   (capacity  1 pool +     pool, kernel  API-level
   ==1 ref)   1 tally set  loop, atomic  merge)
              per thread)  tallies)
   └─────────┴───┬─────┴──────────────┘
        ┌────────┴─────────┐
        │  tally merge (+)  │   associative partial-result merge
        └──────────────────┘
```

- A **chunk** is a contiguous range of global history indices.  All RNG
  streams for history *i* are derived from `(rndseed, rndoffset + f(i))`, so
  a chunk's physics is bit-identical no matter which backend runs it or in
  what order.
- Each backend drains chunks into its own pool(s) and its own partial tally
  set; merging partial tallies is associative, so results are independent of
  scheduling up to floating-point summation order.
- The heterogeneous case (CPU + GPU together) is just dynamic chunk stealing
  from a shared queue, weighted by observed throughput per backend.

## 4. Phased roadmap

Each phase is independently shippable and gated by the equivalence tests in
Phase 0.  Phases 1–3 map directly onto existing TODO.md items.

### Phase 0 — contracts and measurement (prerequisite, small)

1. **Per-history RNG streams.**  Replace the two run-global RNGs with
   per-history derivation: history *i* gets `stream = base + i * K` with
   fixed sub-stream offsets for beam sampling, physics, and (later)
   secondaries.  PCG32 streams already support this; consider a
   counter-based engine (e.g. Philox) later if stream independence needs
   hardening.  Acceptance: results bit-identical for any
   `OSH_TRANSPORT_POOL_CAPACITY`, including capacity 1.
2. **Benchmark harness + perf regression CI.**  A timed run of a small
   analytic case and the case-05 CT case, with steps/s reported, so every
   later phase has a before/after number.
3. **Statistical equivalence checker.**  A tool comparing two BDO outputs
   within MC uncertainty (per-voxel z-test / gamma-style criterion), used as
   a CTest gate for all parallel backends.

### Phase 1 — embarrassingly parallel runs (days of work)

- **Partial-result + merge API** (TODO item): run *N* chunks (in-process,
  sequentially or as separate jobs), save/load partial tallies with history
  counts, merge with correct normalization and (new) per-tally variance for
  uncertainty estimates.
- A small `oshmerge` tool (or `openshieldhit --merge`) for HPC users: SLURM
  array job with `--nstat-offset`/`--chunk` flags, merge afterwards.
- This alone delivers near-linear scaling on clusters with zero concurrency
  risk, and the merge machinery is reused by every later phase.

### Phase 2 — multicore in-process (weeks)

- Worker pool per DEVELOPER.md: **OpenMP** (simplest, supported by MSVC,
  GCC, Clang) or a thin pthreads/Win32 wrapper if finer control is wanted.
  Recommendation: start with OpenMP `parallel for` over chunks; revisit only
  if the scheduler needs work stealing across particle families.
- Per worker: one particle pool, one RNG scratch, one **thread-local tally
  set** (`_Alignas(64)`, per DEVELOPER.md), merged via the Phase 1 API at
  the end of the run.  Memory cost = tally size × threads; if CT-resolution
  grids × 64 threads ever exceed budget, fall back to atomics on a shared
  grid (measure first — thread-local + reduction is usually faster).
- Secondaries from nuclear events stay in the producing worker's pool
  (their RNG streams derive from the parent history index, preserving
  reproducibility).
- Acceptance: per-history trajectories bit-identical to serial; tallies equal
  up to FP summation order; ≥ 0.8 × linear scaling to physical cores on the
  CT benchmark.

### Phase 3 — SIMD inside the wavefront step (weeks, incremental)

Restructure the per-slot scalar loop into **batch phases over the SoA pool**
(TODO: "Batch ion-step phases"):

1. Batch RNG draws for the whole wavefront (`osh_rng_*_vec` already exists).
2. Batch stopping-power / range lookups (TODO: "batch/SIMD lookup helpers").
3. Batch geometry distances (TODO: AVX2 `get_distance_batch`).  For analytic
   GEMCA this has two prerequisites that are worth doing early because they
   also unblock the GPU port: (a) the flattened `insns_flat[]` +
   `insn_begin[]` runtime layout (TODO item) so zone evaluation reads a
   relocatable blob instead of per-zone heap pointers, and (b) **zone-sorted
   batching** — bin live pool slots by current zone ref each wavefront round
   so that SIMD lanes (and later GPU warps) evaluate the same instruction
   stream.  This is the same divergence-reduction technique as Celeritas
   track sorting, and it can be measured on CPU long before any GPU work.
   The voxel Jacobs path is the other high-value target for the CT use case.
4. Vectorized step update: advance positions, apply energy loss, Molière
   MCS, straggling across lanes; handle rare divergent events (nuclear
   interactions, boundary crossings) with a scalar tail pass over a
   compacted index list — the standard masking technique.

Strategy: portable C first (`restrict`, `#pragma omp simd`), then
intrinsics (AVX2/AVX-512/NEON) behind the existing gemca-style runtime
dispatch *only* where profiling shows the compiler failed.  Note that
x86 is not the only CPU target: GH200-based systems (e.g. Cyfronet
Helios GPU partition) pair Hopper GPUs with ARM Grace host CPUs, so the
portable-C path is the primary implementation and AVX2 an optimization,
never the only path.  MCsquare
demonstrates that this combination yields order-of-magnitude gains over
scalar history-based transport on CPUs.

### Phase 4 — GPU backend (months)

Phases 0–3 are deliberate preparation: they force the transport into
kernel-shaped phases over SoA buffers, after which the GPU port is a
translation, not a redesign.

- **Technology choice for a C codebase:**
    - **CUDA first** — best tooling and the dominant install base in both
      clinics and HPC centres; the batch-phase C functions can be compiled
      as `__device__` code with minimal annotation.
    - **HIP port second** — near-mechanical from CUDA; covers AMD HPC
      systems (LUMI, Frontier).
    - **OpenMP target offload** as a portable prototype/fallback — it works
      from plain C (OpenMC proved it at exascale), at some performance cost.
    - SYCL/Kokkos are C++ frameworks and would wrap rather than serve a C
      core; not recommended.
- **Device-resident state:** one large pool (10⁵–10⁶ slots), compiled
  geometry, material tables, HU LUTs, tally grids.  Per-slot RNG state
  (PCG32 is 16 bytes — ideal).  Scoring via `atomicAdd` into device grids;
  copy back and merge through the Phase 1 API, so the GPU is "just another
  backend producing a partial tally".
- **Kernel loop mirrors the wavefront:** fill (on device or host) →
  distance kernel → step kernel → scoring kernel → compact
  (stream-compaction) → refill.  Family queues (the existing scheduler
  seam) and optional energy/material sorting limit divergence.
- **Order of attack:** first an end-to-end **plumbing spike on the simplest
  case** (water box / `00_minimal`: device pool, device RNG, atomic dose
  tally, partial-tally merge — goal is correctness of the pipeline, not
  speed).  Then the **CT/voxel case** — Jacobs traversal over a regular grid
  is GPU-friendly and clinically relevant (FRED, gPMC, MOQUI all do exactly
  this).  Analytic GEMCA follows by uploading the flattened instruction
  blob and reusing the zone-sorted batching from Phase 3 to keep warps
  coherent; if flattening and sorting land early in Phase 3, GEMCA-on-GPU
  moves up nearly for free.

### Phase 5 — dispatcher and HPC integration

- **Heterogeneous scheduling:** CPU pool and GPU(s) pull from one chunk
  queue; chunk size per backend tuned by measured throughput.  Multi-GPU is
  one dispatcher arm per device.
- **Multi-node:** MPI is *optional* — the file-level merge from Phase 1
  composes with SLURM arrays and is operationally simpler.  Add an MPI
  reduction path only if in-job multi-node runs (single output, one job)
  become a real user need.
- CLI surface: `openshieldhit -j N` (threads), `--backend cpu|cuda|hip`,
  `--device i`, with the serial scalar path (`capacity == 1`) kept forever
  as the validation reference.

## 5. Cross-cutting concerns

- **Reproducibility policy:** same `(seed, nstat)` ⇒ identical per-history
  physics on CPU regardless of backend/threads/chunking; tallies identical
  up to FP summation order (document this; offer Kahan/pairwise summation if
  users need tighter determinism).  GPU results match CPU within FP
  tolerance, validated statistically.
- **Testing ladder:** scalar reference (capacity 1) ⇄ wavefront ⇄ threaded
  ⇄ GPU, all on the same RNG streams.  Bitwise trajectory comparison where
  hardware allows, the Phase 0 statistical checker everywhere else.
- **Memory model:** keep the rule "no heap allocation on the hot path";
  workers allocate everything up front.  Watch tally replication × threads
  on CT grids; MOQUI-style sparse scoring is the documented escape hatch.
- **Public API:** extend `osh_simulation_run()` with a run-options struct
  (backend, thread count, chunking) rather than globals, keeping the
  library embeddable (the future DICOM recalculator app will want one
  simulation instance per field across backends).

## 6. Suggested pull-request sequence

Small, independently testable PRs, in dependency order:

1. Per-history RNG streams + capacity-invariance test (Phase 0.1).
2. Benchmark harness + statistical equivalence checker (Phase 0.2–0.3).
3. Scoring partial-result struct, merge API, uncertainty output (Phase 1).
4. Chunked run control + merge tool + HPC docs (Phase 1).
5. OpenMP worker pool with thread-local pools/tallies (Phase 2).
6. Batch-phase refactor of the ion step, scalar (Phase 3, no SIMD yet).
7. SIMD kernels one by one: RNG, dedx lookup, voxel distance, step update.
8. CUDA prototype of the CT benchmark case; merge as partial tally.
9. HIP/OpenMP-target ports; heterogeneous dispatcher.

Every PR ships with the equivalence test as a CTest gate and a before/after
number from the benchmark harness.
