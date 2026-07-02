# Architecture

## Design principle: cold → hot

openshieldhit separates every simulation concern into a **cold** (input/configuration)
stage and a **hot** (runtime/transport) stage.

The app layer (`src/apps/osh/`) owns all file I/O: parsing, resolving paths,
and building the four cold workspaces (`osh_beam_workspace`,
`osh_geometry_workspace`, `osh_material_workspace`, `osh_scoring_workspace`).
These are plain structs in the public API (`include/openshieldhit/`).

The library (`src/beam/`, `src/transport/`, `src/scoring/`, …) never reads
files.  It only sees the cold structs.  `osh_beam_workspace_prepare()`,
`osh_geometry_workspace_prepare()`, etc. compile the cold configuration into
internal prepared state, then `osh_simulation_run()` drives the transport loop.

```
App layer (file I/O)
  └─ parse beam.dat  →  osh_beam_workspace (cold)
  └─ parse geo.dat   →  osh_geometry_workspace (cold)
  └─ parse mat.dat   →  osh_material_workspace (cold)
  └─ parse detect.dat→  osh_scoring_workspace (cold)
          │
          ▼
Library (simulation)
  osh_simulation_create()   — links the four cold workspaces
  osh_simulation_run()      — transport loop
  osh_simulation_save()     — write scored output
```

## Checkpoints, batches, and partial results

`osh_simulation_run()` drives transport as an **outer batch loop** around the
**inner family scheduler** (`src/transport/osh_transport.c`).  The family
scheduler transports a range of ion primaries, then drains every secondary family
(neutrons, and in future fragments/electrons) those primaries banked.  The outer
loop slices the run `[0, nstat)` into batches `[b, b+K)` and runs the family
scheduler once per batch.

A **checkpoint** is the boundary between two batches.  It is the single point
where the whole simulation is *quiescent* and *family-complete*: every family
spawned by the completed primaries has been drained into scoring, so a result
observed there is physically complete — not an ion-only fraction.  This matters
because ratio quantities such as `DLET = Σ(LET·dose)/Σdose` are **biased**, not
merely noisier, if the neutron/fragment contribution is missing from one
`completed_nstat`.  A checkpoint is therefore the correct — and only — place to
observe or dump a partial result.  It is also the natural sync point for the
parallel work to come: per-worker accumulator merges, a variance-batch fold, and
optional periodic dumps all hang off this one boundary.

### The one dial

How often the run checkpoints is the only knob, expressed by
`struct osh_checkpoint_policy` (`src/transport/osh_checkpoint_policy.h`):

| Mode | Cadence | Behaviour |
|------|---------|-----------|
| **Final-only** (default) | — | One batch, `K = nstat`. Fastest; byte-for-byte identical to an unbatched run. |
| **Live** | count (`every_primaries`) | Family-complete batches of a fixed primary count. **Deterministic** and order-independent — the reproducible cadence for tests/CI. |
| **Live** | time (`every_s`) | Family-complete batches sized `≈ rate × cadence`. **Self-bounds overhead** independent of core count — the cadence for production/parallel runs. |

A count cadence is reproducible because each history's RNG stream is a pure
function of its global index, so splitting `[0, nstat)` into fixed sub-ranges and
replaying them in any order reproduces the canonical per-history streams (see
[random numbers](random_numbers.md)).  Scored output then matches the final-only
result up to floating-point reduction order.  A time cadence is **not**
reproducible (batch boundaries fall at wall-clock instants), so a time cadence is
never used where determinism is required.

The speed↔visibility trade is explicit: each checkpoint costs a family drain plus
(later) a barrier + merge, so overhead is `≈ C / T` where `T` is the wall time
between checkpoints.  Final-only pays it once; a time cadence pays it a fixed
number of times per wall-hour on any machine; a count cadence pays it *more often
as the machine gets faster*, which is why it is reserved for deterministic tests.

### Completeness labelling

A checkpoint dump is **exact** (all families drained).  A future escape hatch may
dump at the inner ion safe point without draining secondaries — an **approx**
result — which must be stamped honestly so it is never mistaken for a complete
one (`osh_checkpoint_completeness_label()` → `exact` / `families_pending`).  The
graceful wall-time stop (issue #192) routes its early stop through the family
scheduler, so its partial save is always family-exact for the primaries that
finished.

This module lands the batch-aware seam and the quiescence guarantee only.  The
CLI surface, periodic file dumps and the time cadence (#193), variance-batch
folding (#169), and per-worker accumulator merges (#161) grow into the stable
`osh_checkpoint_policy` shape in their own follow-ups.

## Module map

| Directory | Responsibility |
|-----------|---------------|
| `src/apps/osh/` | File parsers, CLI entry point |
| `src/beam/` | Cold beam model, primary sampling |
| `src/transport/` | Particle transport engine |
| `src/physics/` | Cross-sections, stopping power, nuclear models |
| `src/scoring/` | Detector accumulation, output writers |
| `src/geometry/` | Body/zone ray-tracing |
| `src/gemca/` | GEMCA voxel geometry engine |
| `src/material/` | Material database |
| `src/dicom/` | DICOM reader/writer |
| `src/common/` | Shared utilities (vectors, diagnostics, I/O) |
| `include/openshieldhit/` | Public C API headers |

## Public API boundary

Everything under `include/openshieldhit/` is the stable public API.
Internal headers under `src/` are not API — they may change between releases.

See [API reference →](../api/index.html) for the full Doxygen-generated documentation.

## Examples

`examples/` contains small programs that demonstrate the public C API directly,
bypassing the `openshieldhit` application layer.  They include interactive
geometry viewers (SDL2), a BNCT cell demo, and transport benchmarks.  These are
developer tools, not end-user cases, and are not included in release packages.
SDL2 is required for the interactive viewers (`sudo apt-get install libsdl2-dev`).
See [examples/README.md](../../examples/README.md) for build instructions and a
description of each program.
