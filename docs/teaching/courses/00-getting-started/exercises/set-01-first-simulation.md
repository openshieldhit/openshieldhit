# Set 01: First Simulation

## Learning goals

- Build or locate the `openshieldhit` executable.
- Run and modify a simple water-tank simulation.
- Configure first mesh and cylindrical scorers.
- Plot depth-dose, fluence, and deposited-energy curves.
- Recognize when a scored result depends on the detector geometry.

## Context

The OpenShieldHIT version uses plain input files and command-line runs instead
of a graphical interface.  The intent is the same: start with a simple water
phantom, change one physical or scoring parameter at a time, and explain what
the plotted output means.

The exercises are suggestions.  You may deviate from them, but keep notes about
what you changed and why.

Use the input-file references when needed:

- [`beam.dat`](../../../../user/beam.dat.md)
- [`geo.dat`](../../../../user/geo.dat.md)
- [`mat.dat`](../../../../user/mat.dat.md)
- [`detect.dat`](../../../../user/detect.dat.md)

## Preparation

These exercises assume that `openshieldhit` is installed on your `PATH`, or that
you know where the locally built executable is.  The examples below use
`build/bin/openshieldhit`; replace that path with `openshieldhit` if your
installation provides the command directly.

Start from a copy of the minimal water-tank case:

```bash
build/bin/openshieldhit --dry-run tests/cases/00_minimal/
build/bin/openshieldhit -v -n 5000 tests/cases/00_minimal/
```

For exploratory work, write output into a scratch directory:

```bash
build/bin/openshieldhit --outdir /tmp/osh-set-01 -n 5000 tests/cases/00_minimal/
```

The default `tests/cases/00_minimal/` case is already close to the first task:
protons in water, a simple water geometry, and a 1 mm depth grid.  Increase the
statistics once the setup is correct.

## Exercise 1: Your First Depth-Dose Curve

Set up a water phantom and a proton beam:

- water tank: 20 x 20 x 20 cm3, starting at `z = 0 cm`, symmetric around the
  z-axis;
- surrounding material: air or vacuum;
- beam: 120 MeV protons;
- transverse size: FWHM = 1.0 cm, no angular divergence;
- beam start: `z = -200 cm`;
- statistics: begin with 5000 histories, then repeat with more histories if the
  curve is too noisy.

Useful conversion:

```text
sigma = FWHM / 2.35482
```

For FWHM = 1.0 cm, set a Gaussian beam sigma of about 0.425 cm.

Configure a mesh scorer with 1 mm z-bins over the water tank.  Score the full
20 x 20 cm2 transverse width, so each z-bin is a slab through the whole phantom.

Example detector shape:

```text
Geometry Mesh
    Name DepthSlabs
    X -10.0  10.0    1
    Y -10.0  10.0    1
    Z   0.0  20.0  200

Output
    Filename depth_slabs.dat
    FileFormat TEXT
    Geo DepthSlabs
    Quantity Dose
    Quantity DoseGy
    Quantity Fluence
    Quantity Energy
```

Tasks:

1. Plot dose versus depth `z`.
2. Plot fluence versus depth `z`.
3. Plot deposited energy versus depth `z`.
4. Check the approximate proton range using PSTAR/NIST or another stopping-power
   reference.  Compare that estimate with the simulated Bragg-peak position.
5. Explain why `Dose`, `Fluence`, and `Energy` have different units and different
   normalisations.

Discussion points:

- A full-width slab scorer should collect nearly the whole primary beam.
- The fluence should be of the order of one primary crossing per scoring area,
  modified by lateral spread, reactions, and transport losses.
- `Energy` is not dose.  Dose divides deposited energy by scored mass.

## Exercise 2: A 2D Dose Map

Keep the setup from Exercise 1 and score a 2D dose map in the plane `x = 0 cm`.
Use 1 mm pixels in `y` and `z`.

One practical way to approximate this is a thin mesh slice:

```text
Geometry Mesh
    Name DoseMapXZ
    X -0.05   0.05    1
    Y -10.0  10.0  200
    Z   0.0  20.0  200
```

Tasks:

1. Plot a colour map of dose in the slice.
2. Mark the Bragg peak.
3. Compare the peak depth with the PSTAR/NIST range estimate.
4. If the peak position differs from the simple range estimate, list possible
   causes: finite step size, energy spread, scoring-bin width, material choice,
   or physics switches.

## Exercise 3: Repeat with Carbon Ions

Repeat Exercise 1 with carbon-12 ions:

```text
PRIMARY  6 12
TMAX0    240.0
```

Tasks:

1. Plot carbon dose versus depth.
2. Compare the carbon Bragg peak with the proton curve.
3. Repeat with nuclear reactions disabled and enabled, and note the qualitative
   change in the entrance and tail regions.

Optional extension:

- Insert a 1 x 1 x 1 cm3 heterogeneity halfway into the beam.  Use an
  instructor-provided cortical-bone material if available; otherwise use an
  available dense material as a placeholder and state that it is only a proxy.

## Exercise 4: Where Did the Bragg Peak Go?

This is a scoring exercise, not primarily a physics-model exercise.

Use a larger water phantom, for example 40 cm depth, and a 230 MeV proton pencil
beam:

```text
PRIMARY    Proton
TMAX0      230.0
BEAMSIGMA  0.0  0.0
```

Score deposited energy versus depth twice using cylindrical detectors:

```text
Geometry Cyl
    Name BroadCyl
    R 0.0   5.0    1
    Z 0.0  40.0  400

Geometry Cyl
    Name NarrowCyl
    R 0.0   0.1    1
    Z 0.0  40.0  400
```

Tasks:

1. Plot energy deposition versus depth for the broad cylinder.
2. Plot energy deposition versus depth for the narrow cylinder.
3. Compare both curves.
4. Use a 2D dose map to explain where the Bragg peak appears to have gone in the
   narrow scorer.

Discussion point:

- A small-radius scorer samples only the central part of the beam.  Scattering,
  range straggling, and lateral transport can move energy deposition outside the
  narrow scoring volume.

## Exercise 5: The Divergence Problem

Assume the surrounding medium is vacuum.  Compare three beams that have the same
spot size at the phantom entrance, but different upstream histories:

1. A parallel broad circular beam with radius 5 cm.
2. A beam starting at `z = -100 cm`, with divergence chosen so that its radius is
   5 cm at `z = 0 cm`.
3. A beam starting at the phantom surface with the same divergence as case 2.

Tasks:

1. For each beam, score fluence versus depth with a broad cylindrical scorer,
   for example `R = 0..10 cm`.
2. Repeat with a narrow cylindrical scorer, for example `R = 0..0.2 cm`.
3. Plot the three broad-scorer curves together.
4. Plot the three narrow-scorer curves together.
5. Explain why the broad and narrow scorer comparisons do not tell the same
   story.

Hint:

- The same entrance field does not imply the same phase space.  Position and
  direction distributions both matter.

## Exercise 6: Lateral Penumbras

Use a square beam and compare lateral dose profiles.

Suggested ion cases:

- 236 MeV protons;
- 465 MeV/u carbon-12 ions.

Use a square transverse beam field, for example 5 x 5 cm2, and plot the lateral
dose profile along `x` at about `z = 30 cm`.

Tasks:

1. Configure a mesh scorer with fine x-binning and one thin z-slice.
2. Plot the proton lateral profile.
3. Repeat for carbon ions and overlay the curves.
4. Compare the penumbra qualitatively.

Note:

- The original classroom exercise also compared 14 MeV photons.  Treat that as
  future or external-comparison material until photon transport is part of the
  supported teaching workflow.

## Exercise 7: Particle Energy Spectra

Use differential scoring to record fluence spectra for 238 MeV protons.

Surface spectrum:

```text
Geometry Mesh
    Name Surface
    X -5.0  5.0  1
    Y -5.0  5.0  1
    Z  0.0  0.1  1

Output
    Filename spectrum_surface.dat
    FileFormat TEXT
    Geo Surface
    Quantity Fluence
    Diff1     0.0  250.0  100
    Diff1Type EKIN
```

Near-peak spectrum:

- move the z-slice close to the Bragg peak;
- choose an energy range and binning that resolve the low-energy part of the
  spectrum.

Tasks:

1. Plot `dPhi/dE` versus kinetic energy near the entrance.
2. Plot `dPhi/dE` close to the Bragg peak.
3. Explain why the spectral shape changes with depth.

## Suggested Deliverables

Submit:

- the modified input files;
- one depth-dose plot;
- one fluence plot;
- one 2D dose map;
- one short explanation of the narrow-scorer Bragg-peak problem;
- one spectrum plot, if Exercise 7 was completed.

## Instructor Notes

This set intentionally mixes tool handling with first conceptual questions.  It
belongs in `00 Getting Started` because the emphasis is on running, plotting,
and changing input files.  Later courses should revisit the same phenomena more
systematically:

- random sampling and uncertainty in Course 01;
- phase space and divergence in Course 02;
- stopping power, range, and scattering in Course 03;
- carbon fragmentation and nuclear tails in Course 04.
