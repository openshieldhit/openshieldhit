# Getting started

## Requirements

- CMake ≥ 3.16
- C11 compiler (GCC, Clang, or MSVC)
- Optional: zlib (compressed `.bdz` output files), libSDL2 (interactive geometry viewers)
- DICOM: no external library needed — openshieldhit includes a minimal self-contained reader/writer

## Build

```bash
git clone https://github.com/nbassler/openshieldhit.git
cd openshieldhit
cmake -B build_rel -DCMAKE_BUILD_TYPE=Release .
cmake --build build_rel --parallel
```

The `openshieldhit` binary is written to `build_rel/bin/`.

## Install system-wide

```bash
sudo cmake --install build_rel                   # → /usr/local/bin/openshieldhit (default)
sudo cmake --install build_rel --prefix /usr     # → /usr/bin/openshieldhit
cmake --install build_rel --prefix ~/.local      # user-local, no sudo needed
```

## Run a case

```bash
build_rel/bin/openshieldhit path/to/case/
```

The case directory must contain `beam.dat`, `geo.dat`, `mat.dat`, and
`detect.dat`.  Output files are written into the same directory.

```bash
# parse and compile inputs only, no transport
build_rel/bin/openshieldhit --dry-run path/to/case/

# override history count
build_rel/bin/openshieldhit -n 50000 path/to/case/

# verbose logging
build_rel/bin/openshieldhit -v path/to/case/
```

## Try the minimal example

```bash
build_rel/bin/openshieldhit tests/cases/00_minimal/
```

Produces a Bragg-peak depth-dose profile for 200 MeV protons in water.

## Next steps

- [beam.dat reference](beam.dat.md) — particle, energy, beam geometry, physics switches
- [geo.dat reference](geo.dat.md) — geometry bodies and zones
- [detect.dat reference](detect.dat.md) — scoring detectors
