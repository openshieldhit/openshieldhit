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

OpenShieldHIT also serves as a pedagogical project.

It is particularly aimed at physicists and students who want to understand what happens under the hood of a Monte Carlo transport code.

The project emphasizes hardware-aware programming, with explicit attention to:
- Memory layout
- Cache behavior
- Branching patterns
- CPU efficiency
- Minimal hidden abstractions
The intention is to keep the physics and computational structure inspectable in
the implementation, rather than burying it behind heavy object hierarchies.
The public C API still uses an opaque context handle to preserve a stable
boundary for applications.


## Requirements
[SDL2](https://www.libsdl.org/) is required for building examples (optional): `sudo apt-get install libsdl2-dev`

## How to build

Try:
```bash
cmake -S . -B build && cmake --build build
```

or for debugging:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug && cmake --build build
```

# Run a Minimal Example
```bash
build/bin/openshieldhit tests/cases/00_minimal/
```
Which should produce a `.bdo` and `.txt` output file with scored dose vs. depth, showing a Bragg peak for a proton beam in water.

## Command-line and public API

OpenShieldHIT has a small public C API in `include/openshieldhit/openshieldhit.h`.
The public API is context-based: callers create an opaque
`openshieldhit_context_t`, configure it with `openshieldhit_config_t`, call
`openshieldhit_run()`, and then destroy the context.

The command-line parser in `src/cli/` is an internal front-end. It translates CLI
options into the same public configuration structure, but CLI-specific types are
not part of the public API.

# Try out the examples
```bash
build/bin/bnct_sdl examples/02_bnct/geo_cell.dat
```

## Status

The first end-to-end transport loop is working: proton beams transport through
geometry, accumulate dose in scoring detectors, and produce Bragg-peak depth-dose
distributions.

## TODO
- [x] logger
- [x] vector library
- [x] readline for tag and key parsing
- [x] prng
- [x] geometry parser
- [x] beam parser
- [x] material parser
- [x] detector parser
- [x] raytracer
- [x] straight-line CSDA transport (no scattering)
- [x] stopping power / CSDA range tables
- [x] dose, fluence, LET scoring
- [ ] multiple Coulomb scattering
- [ ] energy straggling (Vavilov/Gaussian)
- [ ] nuclear interactions / secondaries
- [ ] ...


# Disclaimer and License
While conceptually inspired by the application domain of the closed source [SHIELD-HIT12A](https://shieldhit.org), it shares no source code with the original SHIELD-HIT. The architecture, implementation, and code base are completely new and designed according to modern, explicit software engineering principles. Parts of OpenShieldHIT may be reused in SHIELD-HIT12A. OpenShieldHIT is licensed under the MIT License, permitting reuse, modification, and incorporation of its code into other software projects under compatible licensing terms.
