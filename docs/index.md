# OpenShieldHIT

Open-source Monte Carlo particle transport code for shielding calculations and particle therapy physics, built as a clean reimplementation of the SHIELD-HIT12A (SH12A) heritage code.

## What it does

- Transports protons, carbon ions, and other heavy charged particles through arbitrary geometries
- Scores dose, LET, fluence, and other quantities onto voxel or mesh detectors
- Reads clinical DICOM RT Plan and CT data directly
- Provides a pure C library API for embedding in larger workflows

## Quick start

```bash
# build
cmake -B build -DCMAKE_BUILD_TYPE=Release .
cmake --build build --parallel

# run a case (reads geo.dat, beam.dat, mat.dat, detect.dat from the directory)
build/bin/openshieldhit tests/cases/00_minimal/

# dry-run — parse and compile inputs, skip transport
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

## SH12A compatibility

openshieldhit reads the same `beam.dat` / `geo.dat` / `mat.dat` / `detect.dat` format as SH12A.
Where the formats differ or where openshieldhit extends SH12A, the relevant reference pages
call out the difference explicitly.

## API reference

The C library API is documented via Doxygen: [API reference →](api/index.html)
