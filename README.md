# OpenSHIELD-HIT
Monte Carlo Particle Transport

## Requirements
SDL2 is required for building examples (optional): `sudo apt-get install libsdl2-dev`

## How to build
`cmake -S . -B build && cmake --build build`

or for debugging:

`cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug && cmake --build build`

# Try out the examples
```bash
$ build/bin/bnct_sdl examples/02_bnct/geo_cell.dat
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
