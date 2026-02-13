# OpenSHIELD-HIT
Monte Carlo Particle Transport

## How to build
`cmake -S . -B build && cmake --build build`

or for debugging:

`cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug && cmake --build build`


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
