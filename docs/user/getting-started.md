# Getting started

## Requirements

- CMake ≥ 3.16
- C11 compiler (GCC, Clang, or MSVC)
- Optional: zlib (compressed `.bdz` output files), libSDL2 (interactive geometry viewers in `examples/`, not needed for the main application)
- DICOM: no external library needed — openshieldhit includes a minimal self-contained reader/writer

## Build

```bash
git clone https://github.com/nbassler/openshieldhit.git
cd openshieldhit
cmake --preset release
cmake --build --preset release
```

The `openshieldhit` binary is written to `build/bin/`.

## Install system-wide

```bash
sudo cmake --install build                   # → /usr/local/bin/openshieldhit (default)
sudo cmake --install build --prefix /usr     # → /usr/bin/openshieldhit
cmake --install build --prefix ~/.local      # user-local, no sudo needed
```

## Run a case

```bash
build/bin/openshieldhit path/to/case/
```

The case directory must contain `beam.dat`, `geo.dat`, `mat.dat`, and
`detect.dat`.  Output files are written into the same directory.

```bash
# parse and compile inputs only, no transport
build/bin/openshieldhit --dry-run path/to/case/

# override history count
build/bin/openshieldhit -n 50000 path/to/case/

# verbose logging
build/bin/openshieldhit -v path/to/case/
```

## Try the minimal example

```bash
build/bin/openshieldhit tests/cases/00_minimal/
```

Produces a Bragg-peak depth-dose profile for 200 MeV protons in water.

## Next steps

- [beam.dat reference](beam.dat.md) — particle, energy, beam geometry, physics switches
- [geo.dat reference](geo.dat.md) — geometry bodies and zones
- [detect.dat reference](detect.dat.md) — scoring detectors
