# Parallelisation strategies — discussion deck

Beamer + TikZ slides on parallelisation strategies and memory-access
patterns for OpenShieldHIT. Prepared for a design discussion with Niels
Bassler.

## Contents

- `slides.tex` — source (self-contained; every diagram is TikZ/pgfplots).
- `slides.pdf` — compiled deck (15 slides).

## What it covers

1. The problem: `N` independent histories; compute is embarrassingly
   parallel, scoring and memory layout are the hard parts.
2. The `cold → hot` architecture and why prepared cold state is shareable.
3. Per-history RNG streams as the reproducibility foundation.
4. Four strategies, on a shared → shared-writable axis:
   - **A** Multi-node / MPI (independent runs, merge `.bdo`)
   - **B** Threads + per-thread scoring pages (reduction at the end)
   - **C** Threads + atomic scoring (contention / false sharing)
   - **D** SIMD / GPU on the SoA particle pool
5. Memory-access diagrams: shared vs. private, the cache hierarchy, and
   host-aware memory budgeting (#152/#153).
6. Comparison matrix and a suggested incremental path.

Every diagram is grounded in the code base — sources are cited on each
slide (`src/random/osh_rng.h`, `src/common/osh_particle_pool.h`,
`src/scoring/README.md`, `docs/dev/architecture.md`,
`src/apps/osh/osh_membudget.*`).

## Building

Needs a TeX distribution with Beamer, TikZ and pgfplots
(`texlive-latex-recommended`, `texlive-pictures`, `texlive-fonts-recommended`,
`lmodern`):

```sh
pdflatex slides.tex   # run twice for the section navigation / page refs
```
