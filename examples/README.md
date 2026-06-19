# OpenShieldHIT Examples

Four standalone C programs demonstrating geometry inspection and performance
benchmarking.  All are built automatically by the standard CMake preset.

## Requirements

SDL2 is needed for examples 01 and 02:

```bash
sudo apt-get install libsdl2-dev   # Debian/Ubuntu
brew install sdl2                  # macOS
```

Examples 03 and 04 have no additional dependencies.

## Building

```bash
cmake --preset release && cmake --build --preset release --parallel
```

Binaries land in `build/bin/`.

## 01 — SDL geometry viewer

Interactive GEMCA geometry inspection window.  Several geometry files are
included; the default is a simple RPP arrangement:

```bash
build/bin/gemca_sdl_viewer examples/01_sdl_viewer/geo.dat
build/bin/gemca_sdl_viewer examples/01_sdl_viewer/geo_RCC03.dat
```

## 02 — BNCT cell geometry

SDL visualisation of a BNCT cell geometry:

```bash
build/bin/bnct_sdl examples/02_bnct/geo_cell.dat
```

## 03 — GEMCA benchmark

Geometry evaluation throughput benchmark.  Prints ray-body intersection rates
to stdout; useful for comparing algorithm variants or hardware:

```bash
build/bin/gemca_bench
```

## 04 — Raycast benchmark

Raycasting throughput benchmark over a voxel grid:

```bash
build/bin/gemca_raycast_bench
```
