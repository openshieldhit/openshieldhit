# OpenShieldHIT

[![CMake CI](https://github.com/openshieldhit/openshieldhit/actions/workflows/test.yml/badge.svg)](https://github.com/openshieldhit/openshieldhit/actions/workflows/test.yml)
[![clang-format](https://github.com/openshieldhit/openshieldhit/actions/workflows/clang-format.yml/badge.svg)](https://github.com/openshieldhit/openshieldhit/actions/workflows/clang-format.yml)
[![clang-tidy](https://github.com/openshieldhit/openshieldhit/actions/workflows/clang-tidy.yml/badge.svg)](https://github.com/openshieldhit/openshieldhit/actions/workflows/clang-tidy.yml)
[![codecov](https://codecov.io/gh/openshieldhit/openshieldhit/branch/main/graph/badge.svg)](https://codecov.io/gh/openshieldhit/openshieldhit)

**Supported platforms:** Linux · macOS · Windows — all three are tested on every commit via CI.


OpenShieldHIT is a modern and lean Monte Carlo particle transport framework written entirely from scratch in C.

## Philosophy

OpenShieldHIT is designed to be:
- Lean
- Fast
- Efficient
- Readable
- Modular
- Easy to extend

It is not intended to become a large, monolithic framework like Geant4. Instead, the goal is to maintain a lightweight and transparent code base where the full transport chain can be understood and modified without navigating millions of lines of abstraction.

## Pedagogical Goal

OpenShieldHIT also serves as a pedagogical project, particularly aimed at
physicists and students who want to understand what happens under the hood of a
Monte Carlo transport code.

The project emphasizes hardware-aware programming with explicit attention to
memory layout, cache behavior, branching patterns, and CPU efficiency.  The
goal is that the full transport chain — from input parsing through geometry
traversal, energy-loss stepping, and scoring — can be traced through a small
number of readable files.

Opaque handles exist at API boundaries to keep the public interface stable, but
they hide only plumbing, not physics.  The stepping logic, geometry compiler,
and scoring accumulation are all open and traceable, without navigating millions
of lines of framework infrastructure.

## Requirements
[SDL2](https://www.libsdl.org/) is required for building examples (optional): `sudo apt-get install libsdl2-dev`

## How to build

```bash
cmake --preset release && cmake --build --preset release   # optimised (-O3)
cmake --preset debug   && cmake --build --preset debug     # debug symbols, -Og
```

Available presets: `debug`, `release`, `relwithdebinfo`, `prof`.
Requires CMake ≥ 3.21.  Binaries land in `build_rel/bin/` (release) or `build/bin/` (debug).

Useful cache variables for development:

```bash
-DOSH_RAYTRACE_ALGORITHM=SIMPLE|SIDDON|JACOBS
-DOSH_VOXEL_LAYOUT=ROW_MAJOR|MORTON8
```

`cmake --install` currently installs the main executable, documentation, and
selected examples. Public library/header installation is still incomplete.

## Run a Minimal Example
```bash
build_rel/bin/openshieldhit tests/cases/00_minimal/
```
Produces a `.bdo` and `.txt` output file with scored dose vs. depth — a Bragg peak for a proton beam in water.

## Run a CT Transport Example

Case 05 transports 126 MeV protons through the DCPT head-phantom CT
(512 × 512 × 177 voxels, Schneider 2000 HU calibration) and scores energy onto
the full CT grid as well as an RTDOSE template file.

Validate parsing and geometry compilation without transporting particles:

```bash
build_rel/bin/openshieldhit --dry-run tests/cases/05_dicom_simple/
```

Run 100 000 primaries with verbose output:

```bash
build_rel/bin/openshieldhit -vvv -n 100000 tests/cases/05_dicom_simple/
```

Output files are written next to the input case:

| File | Content |
|------|---------|
| `NB_ddc.dat` | depth-dose profile along the beam axis (text) |
| `NB_XZ.bdo` | XZ fluence/energy plane (BDO) |
| `NB_ct.bdo` | energy scored onto the full CT voxel grid, ~355 MB (BDO) |
| `NB_rtdose.dcm` | RTDOSE round-trip: scored energy replaces pixel data in the template `.dcm` |

Overlay the RTDOSE output on the CT with the bundled Python tool:

```bash
python3 tools/plot_dicom.py \
    tests/fixtures/dicom/DCPT_headphantom/ \
    tests/cases/05_dicom_simple/NB_rtdose.dcm
```

Add `-o result.png` to save instead of displaying interactively.
Requires `numpy`, `matplotlib`, and `pydicom`.


## The `openshieldhit` application

The main deliverable is the `openshieldhit` (or `openshieldhit.exe` on Windows)
executable built from `src/apps/osh/`.  It reads four plain-text input files from
a working directory and produces scored output:

```bash
openshieldhit path/to/case/          # reads geo.dat, beam.dat, mat.dat, detect.dat
openshieldhit --validate path/to/case/   # parse and validate without running
openshieldhit --help
```

## Public C API

The library API lives under `include/openshieldhit/`.  Input descriptions are
held in plain, caller-owned structs (cold workspaces); the library compiles them
into an opaque simulation handle at run time:

```c
#include "openshieldhit/simulation.h"

// 1. Fill cold workspaces (via parsers or programmatically)
struct osh_beam_workspace     *beam    = /* ... */;
struct osh_geometry_workspace *geo     = /* ... */;
struct osh_material_workspace *mat     = /* ... */;
struct osh_scoring_workspace  *scoring = /* ... */;

// 2. Compile into a simulation (zone→material wiring, runtime tables, etc.)
struct osh_simulation *sim;
osh_simulation_create(beam, geo, mat, scoring, &sim);

// 3. Run transport and save outputs
osh_simulation_run(sim, "/path/to/output/");

// 4. Release — cold workspaces are NOT freed here, caller owns them
osh_simulation_free(sim);
osh_beam_workspace_free(beam);
// ...
```

The runtime representations are private to `src/simulation/` and never visible
to calling code.

## Try out the examples
```bash
build_rel/bin/bnct_sdl examples/02_bnct/geo_cell.dat
```

## Status

The current codebase has a working end-to-end ion transport path with:

- text-input parsing through the `openshieldhit` application
- cold-to-runtime compilation for beam, geometry, materials, and scoring
- analytic GEMCA geometry plus DICOM CT-backed voxel body setup (`DCM`)
- multiple raytrace implementations with a shared grid contract
- dose / fluence scoring with BDO, text, and DICOM RTDOSE output

CT voxel transport is working end-to-end:

- Schneider and Permatassari HU calibration/material registration
- DICOM CT parsed into a voxel body; HU table applied per-voxel at transport time
- compile-time selectable voxel storage layout (`ROW_MAJOR` baseline, `MORTON8` optional)
- `DicomCT` and `DicomRTDOSE` scoring geometries converted to Cartesian meshes at
  the application layer; the transport library remains DICOM-agnostic
- RTDOSE round-trip: read template file, score energy onto its grid, write back `.dcm`
- IEC 61217 gantry and couch rotation: patient body correctly rotated during
  transport for any HFS/HFP/FFS/FFP/HFDL/HFDR/FFDL/FFDR position and angle


The detailed implementation roadmap lives in [TODO.md](TODO.md).

## Repository layout

- `include/openshieldhit/` — public headers
- `src/apps/osh/` — CLI application, file parsing, app-side orchestration
- `src/common/` — shared low-level utilities, math, raytrace, and helpers
- `src/gemca/` — geometry preparation and runtime compilation
- `src/material/`, `src/beam/`, `src/scoring/`, `src/transport/`, `src/simulation/` — core runtime modules


# Disclaimer and License
While conceptually inspired by the application domain of the closed source [SHIELD-HIT12A](https://shieldhit.org), it shares no source code with the original SHIELD-HIT. The architecture, implementation, and code base are completely new and designed according to modern, explicit software engineering principles. Parts of OpenShieldHIT may be reused in SHIELD-HIT12A. OpenShieldHIT is licensed under the MIT License, permitting reuse, modification, and incorporation of its code into other software projects under compatible licensing terms.
