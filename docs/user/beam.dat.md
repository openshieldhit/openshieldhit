# beam.dat reference

`beam.dat` defines the primary particle, beam geometry, physics models, and
simulation statistics.  It is a plain-text keyword file: one keyword per line,
followed by its arguments.  Lines beginning with `#` or `*` are comments.
Keywords are case-insensitive.

```
# example beam.dat
PRIMARY         Proton
TMAX0           200.0       0.0     # 200 MeV, no energy spread
BEAMPOS         0.0  0.0  -50.0
BEAMSAD         200.0  256.0
USECBEAM        sobp.dat
NSTAT           100000      10000
RNDSEED         12345
DELTAE          0.005
DEMIN           0.025
STRAGG          2
MSCAT           2
NUCRE           1
NEUTRLCUT       0.0
```

---

## Particle definition

### PRIMARY

```
PRIMARY  <name>
PRIMARY  <pdg>
PRIMARY  <Z> <A>
```

Defines the primary particle species.

| Form | Example | Notes |
|------|---------|-------|
| Name | `Proton`, `Carbon`, `Neutron` | Case-insensitive |
| PDG code | `2212` (proton), `1000060120` (¹²C) | Standard Monte Carlo PDG numbering |
| Z A | `6 12` | Atomic number and mass number |

### TMAX0

```
TMAX0  <value>  [<spread>]
```

Initial energy or momentum of the primary, and optional energy spread (σ).

| Sign of `value` | Interpretation |
|-----------------|----------------|
| positive | Kinetic energy T₀ [MeV/nucleon for ions, MeV otherwise] |
| negative | Total momentum p₀ [MeV/c], stored as \|value\| |

The same sign convention applies to `spread`.  Both energy and momentum forms are
accepted; the complementary quantity is derived automatically from the particle mass.

```
TMAX0   200.0          # 200 MeV/nucleon, no spread
TMAX0   200.0   2.0    # 200 MeV/nucleon, σ = 2 MeV/nucleon
TMAX0  -1000.0         # p₀ = 1000 MeV/c
```

---

## Beam geometry

### BEAMPOS

```
BEAMPOS  <X>  <Y>  <Z>   [cm]
```

Physical starting position of the beam in UNIVERSE coordinates [cm].
This is the point where primary particles are created.

- For a beam entering along +Z, a negative Z value places the start
  upstream (in front of) the geometry.
- When `USECBEAM` is active, the X and Y here apply to the single-spot
  template only; USECBEAM spot positions override them.  The Z value is
  always used as the beam-entrance plane.

### BEAMSIGMA

```
BEAMSIGMA  <sx>  [<sy>]   [cm, 1σ]
```

Transverse beam spot size (1σ half-width) [cm].  One value gives a round beam;
two values give an elliptical beam.  The sign encodes the spatial profile:

| Value | Shape |
|-------|-------|
| both positive | Gaussian (σ_x, σ_y) |
| both negative | Uniform square (half-widths \|sx\|, \|sy\|) |
| first ≥ 0, second < 0 | Uniform circular disk, radius \|sy\| |
| both zero | Pencil beam (point source) |

```
BEAMSIGMA   0.5          # round Gaussian, σ = 0.5 cm
BEAMSIGMA   0.5   0.3    # elliptical Gaussian
BEAMSIGMA  -1.0  -1.0    # uniform 2×2 cm square
BEAMSIGMA   0.0  -1.0    # uniform circular disk, r = 1 cm
```

### BEAMSAD

```
BEAMSAD  <sadX>  [<sadY>]   [cm]
```

Source-to-Axis Distance (SAD): the distance from the upstream virtual point
source to the isocenter [cm].

SAD is always a **positive** number — it is a physical distance along the beam
axis, not a signed coordinate.  It is independent of the beam direction (a beam
entering at −Z and a beam entering at +Z both have a positive SAD pointing
upstream).

One value gives a symmetric nozzle (same SAD for X and Y).  Two values support
asymmetric nozzles (common in clinical proton/carbon gantries).

The special value `INF` selects parallel beam delivery explicitly (no fan-out
correction, equivalent to SAD = ∞).  Use this to suppress the warning that
appears when `USECBEAM` is present without `BEAMSAD`.

```
BEAMSAD   200.0            # symmetric, SAD = 200 cm (2 m)
BEAMSAD   200.0  256.0     # asymmetric: SAD_x = 200 cm, SAD_y = 256 cm
BEAMSAD   INF              # parallel beam, no fan-out correction
```

!!! note "Required for correct spot-scanning geometry"
    When `USECBEAM` is active and the nozzle has a finite source distance,
    `BEAMSAD` must be set.  Without it, openshieldhit assumes parallel delivery (SAD = ∞)
    and emits a warning.  See [USECBEAM](#usecbeam) for the full coordinate
    convention.

### USECBEAM

```
USECBEAM  <filename>
```

Load a spot-scanning delivery from an external file (one spot per line).
Path is resolved relative to the directory containing `beam.dat`.

#### Spot file format

Whitespace-delimited numeric text, one spot per row.  All rows must have the
same number of columns.  Accepted column layouts:

| Columns | Layout |
|---------|--------|
| 5 | `E  X  Y  FWHM  weight` |
| 6 | `E  X  Y  FWHMx  FWHMy  weight` |
| 7 | `E  dE  X  Y  FWHMx  FWHMy  weight` |
| 9 | `E  dE  X  Y  FWHMx  FWHMy  divX  divY  weight` |
| 11 | `E  dE  X  Y  FWHMx  FWHMy  divX  divY  corX  corY  weight` |

Units:

| Column | Unit |
|--------|------|
| E, dE | GeV/nucleon |
| X, Y | cm — **isocenter coordinates** (see below) |
| FWHM | cm |
| divX, divY | mrad |
| corX, corY | dimensionless (−1 to 1) |
| weight | arbitrary (relative spot weight / MU) |

#### Coordinate convention

**Spot X/Y positions are isocenter coordinates** — the lateral offset of each
spot at the isocenter plane (z = 0 in beam-local coordinates).  This matches
the SH12A legacy format and the DICOM RT Plan standard, where scanning-magnet
setpoints are always expressed at isocenter.

openshieldhit back-projects these to the beam-entrance plane using `BEAMSAD` and `BEAMPOS`
before storing them in the internal cold struct.  The `osh_beam_spot.p[]` field
always holds physical beam-start coordinates.

Back-projection formula (per axis):

```
x_beam_start = x_iso × (SAD + z_start) / SAD
```

where `z_start` = `BEAMPOS` Z (negative = upstream of isocenter).

**Example:** x_iso = 5 cm, z_start = −50 cm, SAD = 200 cm  
→ factor = 150/200 = 0.75 → x_beam_start = 3.75 cm  
Without back-projection the SAD fan-out would make the field ~33% too wide.

!!! warning "Parallel beam without BEAMSAD"
    If `BEAMSAD` is absent, openshieldhit uses the spot positions as-is (no
    back-projection) and emits a warning.  This is physically correct only
    when no scanning-magnet fan-out is present (broad research fields,
    parallel pencil beams).  Set `BEAMSAD INF` to suppress the warning
    when parallel delivery is intentional.

---

## Physics switches

### STRAGG

```
STRAGG  <mode>
```

Energy straggling model.

| Mode | Model |
|------|-------|
| 0 | Off |
| 1 | Gaussian (Bohr straggling) |
| 2 | Vavilov (statistically correct for thin absorbers; **recommended**) |

Vavilov converges to Gaussian for thick absorbers, so mode 2 is the safe
default.  Gaussian is faster and acceptable when the absorber is thick
compared to the particle range.

### MSCAT

```
MSCAT  <mode>
```

Multiple Coulomb scattering model.

| Mode | Model |
|------|-------|
| 0 | Off |
| 1 | Gaussian (Rossi–Greisen approximation, fast) |
| 2 | Molière (more accurate; **recommended** for clinical use) |

### NUCRE

```
NUCRE  <0|1|2|3>
```

Nuclear reactions switch.  This controls **hadronic** (nuclear) interactions
only.  Coulomb (atomic) elastic scattering is always handled separately by the
multiple Coulomb scattering model selected with [`MSCAT`](#mscat) and is
unaffected by this switch.

| Value | Effect |
|-------|--------|
| 0 | Off — electromagnetic transport only, no hadronic interactions |
| 1 | All nuclear reactions — Tripathi nuclear inelastic absorption and pp nuclear elastic scattering with secondary proton production |
| 2 | Elastic only — pp nuclear elastic scattering only; no inelastic absorption |
| 3 | Inelastic only — Tripathi nuclear inelastic absorption; diagnostic mode |

Mode 1 is the most complete current setting for proton transport in tissue:
inelastic reactions on oxygen and other nuclei remove primaries from the beam
(Tripathi cross section), while pp nuclear elastic scattering deflects primaries
and produces recoil proton secondaries scored at generation ≥ 1.

Mode 2 is useful for isolating the pp elastic contribution, for example when
cross-validating secondary proton spectra against analytic estimates.

Mode 3 is useful for isolating inelastic primary attenuation without pp elastic
energy transfer.

Disabling nuclear reactions entirely (`NUCRE 0`) is useful for pure
range/dose validation where hadronic secondaries would complicate comparison
with analytic models.

### NEUTRLCUT

```
NEUTRLCUT  <e>   [MeV]
```

Lower energy cutoff for neutron transport [MeV].  Neutrons below this
threshold are killed and their energy deposited locally.  `0.0` disables the
cutoff (all neutrons transported).

### DELTAE

```
DELTAE  <f>
```

Maximum relative energy loss per transport step (dimensionless fraction,
e.g. `0.005` = 0.5%).  Smaller values improve accuracy at the cost of
more steps and longer run time.

### DEMIN

```
DEMIN  <e>   [MeV/nucleon]
```

Minimum kinetic energy step size [MeV/nucleon].  Prevents the step-size
algorithm from taking extremely small steps near the Bragg peak.

---

## Statistics and reproducibility

### NSTAT

```
NSTAT  <n>  [<step>]
```

Total number of primary particle histories, and optional save interval.

- `n` — total histories to simulate (required)
- `step` — write intermediate results to disk every `step` histories; `−1`
  disables intermediate saves (write only at the end)

### RNDSEED

```
RNDSEED  <seed>
```

Random number generator seed (integer).  Fixing the seed makes a run
bit-for-bit reproducible.  Using different seeds produces statistically
independent runs — the standard approach for Monte Carlo uncertainty estimation.

---

## Beam mode modifiers (advanced)

| Key | Status | Notes |
|-----|--------|-------|
| `BEAMDIV` | supported | Angular divergence [mrad] and optional focus position |
| `BEAMDIR` | supported | Beam direction angles θ, φ [rad] |
| `USEPARLEV` | not yet implemented | Parallel lever beam optics file |
| `USEBMOD` | not yet implemented | Ripple filter / beam modifier |
| `BMODMC` | stub | Logs a warning, ignored |
| `BMODTRANS` | deprecated | Logs a warning, ignored |
| `EXTSPEC` | supported | External energy spectrum file |
| `MAKELN` | supported | Write primary phase-space to file (0/1) |
| `TCUT0` | supported | Kinetic energy cutoff for primary particle |
| `APCORR` | supported | AP correction factor |
