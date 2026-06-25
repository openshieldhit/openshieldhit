# Build & run — quick reference

Copy-paste commands for building OpenShieldHIT from source and running it.
Kept short on purpose: `head INSTALL.md` (or `cat`) should show you everything
you need without scrolling. For the full story see the links at the bottom.

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

## Run

Run a case (reads `geo.dat`, `beam.dat`, `mat.dat`, `detect.dat` from the dir):

```bash
build/bin/openshieldhit path/to/case/
```

Parse and load inputs only, skip transport:

```bash
build/bin/openshieldhit --dry-run path/to/case/
```

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

## Install

Install into your home directory — no sudo needed (`~/.local/bin/openshieldhit`):

```bash
cmake --install build --prefix ~/.local
```

System-wide (`/usr/local/bin/openshieldhit`):

```bash
sudo cmake --install build
```

## More

- [`README.md`](README.md) — project overview and examples
- [`DEVELOPER.md`](DEVELOPER.md) — full developer guide (sanitizers, coverage, profiling)
- [`docs/dev/build.md`](docs/dev/build.md) / [`docs/dev/testing.md`](docs/dev/testing.md) — build & test details
- [`docs/user/getting-started.md`](docs/user/getting-started.md) — step-by-step user guide
- Full docs: <https://openshieldhit.github.io/openshieldhit>
