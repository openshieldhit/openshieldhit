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
developer tools, not end-user cases.  SDL2 is required for the interactive
viewers (`sudo apt-get install libsdl2-dev`).  See
[examples/README.md](../../examples/README.md) for build instructions and a
description of each program.
