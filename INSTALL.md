# Build & run — quick reference

Requirements: CMake ≥ 3.21 and a C11 compiler (GCC, Clang, or MSVC).

## Build

Optimised build (`-O3`), binary in `build/bin/`:

```bash
cmake --preset release && cmake --build --preset release --parallel
```

Debug build (`-O0 -g`), binary in `build_debug/bin/`:

```bash
cmake --preset debug && cmake --build --preset debug --parallel
```

Presets: `release`, `debug`, `relwithdebinfo`, `prof`.


## Install

Install into your home directory — no sudo needed (`~/.local/bin/openshieldhit`):

```bash
cmake --install build --prefix ~/.local
```

System-wide (`/usr/local/bin/openshieldhit`):

```bash
sudo cmake --install build
```

## Run

Run the bundled minimal example:

```bash
build/bin/openshieldhit -v tests/cases/00_minimal/
```

Show all command-line options:

```bash
build/bin/openshieldhit --help
```

## Test

```bash
cd build && ctest --output-on-failure -j4
```

## TL;DR
```
$ cmake --preset release && cmake --build --preset release --parallel; cmake --install build --prefix ~/.local
```
