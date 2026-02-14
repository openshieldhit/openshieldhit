# OpenShieldHIT

Lean Monte Carlo Particle Transport in C

OpenShieldHIT is a modern Monte Carlo particle transport framework written entirely from scratch in C.


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
The intention is to expose the physics and the computational structure, rather than obscuring it behind heavy object hierarchies or opaque frameworks.Monte Carlo Particle Transport


## Requirements
SDL2 is required for building examples (optional): `sudo apt-get install libsdl2-dev`

## How to build

Try:
```bash
cmake -S . -B build && cmake --build build
```

or for debugging:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug && cmake --build build
```

# Try out the examples
```bash
build/bin/bnct_sdl examples/02_bnct/geo_cell.dat
```

## TODO
- [x] logger
- [x] vector library
- [x] readline for tag and key parsing
- [x] prng
- [x] geometry parser
- [ ] beam parser
- [ ] material parser
- [ ] detector parser
- [x] raytracer
- [ ] transport
- [x] particle data
- [ ] stopping power / range / optical depths / restricted stopping power
- [ ] ion scattering
- [ ] vavilov straggling
- [ ] ...


# Disclaimer and License
While conceptually inspired by the application domain of the closed source [SHIELD-HIT12A](https://shieldhit.org), it shares no source code with the original SHIELD-HIT. The architecture, implementation, and code base are completely new and designed according to modern, explicit software engineering principles. OpenShieldHIT is licensed under the MIT License, permitting reuse, modification, and incorporation of its code into other software projects under compatible licensing terms.