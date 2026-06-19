# OpenShieldHIT

Open-source Monte Carlo particle transport code for shielding calculations and particle therapy physics.


## What it does

- Transports protons, carbon ions, and other heavy charged particles through arbitrary geometries
- Scores dose, LET, fluence, and other quantities onto voxel or mesh detectors
- Reads clinical DICOM RT Plan and CT data directly
- Provides a pure C library API for embedding in larger workflows

## Quick start

Build:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release .
cmake --build build --parallel
```

Run a case (reads `geo.dat`, `beam.dat`, `mat.dat`, and `detect.dat` from the directory):

```bash
build/bin/openshieldhit tests/cases/00_minimal/
```

Write outputs to a separate directory (created automatically if needed):

```bash
build/bin/openshieldhit --outdir /tmp/osh-run tests/cases/00_minimal/
```

Dry-run — parse and compile inputs, skip transport:

```bash
build/bin/openshieldhit --dry-run tests/cases/00_minimal/
```

## Input files

Every simulation case lives in a directory with four plain-text files:

| File | Purpose |
|------|---------|
| [`beam.dat`](user/beam.dat.md) | Particle type, energy, beam geometry, physics switches |
| [`geo.dat`](user/geo.dat.md) | Geometry definition (bodies, zones, or DICOM CT) |
| [`mat.dat`](user/mat.dat.md) | Material compositions |
| [`detect.dat`](user/detect.dat.md) | Scoring geometries and quantities |

## Command line

- [`user/command-line.md`](user/command-line.md) — CLI options, `--outdir`, dry-run, and file overrides

## SH12A compatibility

Conceptually inspired by the same application domain as [SHIELD-HIT12A](https://shieldhit.org) (SH12A) and sharing its input-file format.

openshieldhit reads the same `beam.dat` / `geo.dat` / `mat.dat` / `detect.dat` format as SH12A.
Where the formats differ or where openshieldhit extends SH12A, the relevant reference pages
call out the difference explicitly.


## API reference

The C library API is documented via Doxygen: [API reference →](api/index.html)
