# Parallelisation strategies — discussion deck

Beamer + TikZ slides on parallelisation strategies and memory-access
patterns for OpenShieldHIT. Prepared for a design discussion with Niels
Bassler.

## Contents

- `slides.tex` — source (self-contained; every diagram is TikZ/pgfplots).
- 25 slides. The compiled PDF is intentionally **not** committed (the repo
  `.gitignore` excludes `*.pdf`); build it locally with `pdflatex`.

## What it covers

1. The shape of the problem: `N` independent histories; compute is
   embarrassingly parallel (80–86 % in the transport step), so scoring
   contention, memory layout, and load balance are the hard parts.
2. A short **vocabulary** slide — cold/hot, RNG stream, SoA/AoS, the 64 B
   cache line — for readers new to the terms.
3. **How much parallelism:** the target machines — laptop (~10 cores),
   Cyfronet Ares (48), Athena (128, GPU-first), Helios (192), at ~2–4 GB
   RAM/core — and why memory-per-core, not core count, drives the choice.
4. **Two workloads that parallelise differently:** CSG geometry + mesh/diff
   scoring (compute-bound, replicate freely) vs. CT / patient plan + CT
   scoring (memory-bound, the ×T replication is the binding constraint).
5. The `cold → hot` architecture and per-history RNG streams (the foundation
   that makes any partition valid).
6. Four strategies on a shared → shared-writable axis:
   - **A** Multi-node / MPI (independent runs, merge `.bdo`)
   - **B** Threads + per-thread accumulators (reduction at the end)
   - **C** Threads + atomic scoring — including the *"each thread only touches
     a fragment"* question and why index-partitioning still collides on the
     Bragg peak
   - **D** SIMD / GPU on the SoA pool, with a practical AVX2 / AVX-512 / GPU-warp
     width slide
7. Memory-access diagrams: shared vs. private, the cache hierarchy, and
   host-aware memory budgeting (#152/#153).
8. **The merge as the variance contract** (#169 — batch statistics,
   Welford/West + Schubert–Gertz, unequal batches).
9. **Run-control** (#170 — periodic/on-demand dumps, wall-time budget,
   graceful stop; each dump is a batch boundary) and the **WASM web target**
   (#171 — Web Workers as the MPI shape, streaming partial results).
10. The **reproducibility contract** (#168), the effect of **growing
    cross-section data** for future neutron/photon transport (#154/#111), and
    **licensing & threading-library choices** (MIT core, OpenMP/MPI/GPU).
11. Comparison matrix and a suggested incremental path with open questions.

Diagrams are grounded in the code base, with per-slide source citations
(`src/random/osh_rng.h`, `src/common/osh_particle_pool.h`,
`src/gemca/README.md`, `src/voxel/README.md`, `src/scoring/README.md`,
`docs/dev/architecture.md`, `src/common/osh_sysinfo.*`,
`src/apps/osh/osh_membudget.*`) and the parallelism issues/PR (#162,
#168–#171, #154).

## Building

Needs a TeX distribution with Beamer, TikZ and pgfplots
(`texlive-latex-recommended`, `texlive-pictures`, `texlive-fonts-recommended`,
`lmodern`):

```sh
pdflatex slides.tex   # run twice for the section navigation / page refs
```
