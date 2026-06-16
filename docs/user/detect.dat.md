# detect.dat reference

`detect.dat` defines scoring geometries, quantities, optional filters/settings,
and output formats.

## Geometry types

### Mesh

Cartesian `(X, Y, Z)` scoring grid.

```text
Geometry Mesh
    Name MyMesh
    X  -5.0   5.0   10    # lo  hi  nbins
    Y  -5.0   5.0   10
    Z   0.0  20.0  200
```

Each axis line gives the lower bound, upper bound, and number of equally spaced
bins. Axis order in the file does not matter.

### Cyl

Cylindrical `(R, Z)` scoring grid, rotationally symmetric around its local
`Z` axis.

```text
Geometry Cyl
    Name MyCyl
    R   0.0   5.0    5    # lo  hi  nbins
    Z   0.0  20.0  200
```

Rules and semantics:

- `R` is the radial coordinate in cm, from inner radius to outer radius.
- `Z` is the axial coordinate in cm.
- Both `R` and `Z` must be present; their declaration order does not matter.
- `Geometry Cyl` is an `R/Z` detector only. There is currently no explicit
  `phi` binning; the detector always represents the full azimuth.
- Without `Rotation`, the cylinder axis is aligned with the universe `Z` axis
  and centered on `(x, y) = (0, 0)`.
- `Rotation <theta_deg> <phi_deg>` applies a pure rotation from universe space
  into the detector's local frame. It does not apply any translation. Axis
  bounds remain local coordinates.

Voxel volumes are computed exactly per radial bin:

`V(ir) = pi * (r1^2 - r0^2) * dz`

This exact annular volume is used for `Fluence`, `Dose`, and related scorers, so
large radial bins are normalized correctly.

## Output

Example:

```text
Output
    Filename NB_cyl.dat
    FileFormat TEXT
    Geo MyCyl
    Quantity Energy
    Quantity Fluence
```

For `Geometry Cyl`:

- `TEXT`/`ASCII` output is supported.
- `BDO2019` output is supported.
- ASCII rows are written in local detector coordinates with columns
  `Z R <QUANTITIES...>`.
- Flat voxel order is `idx = ir + nr * iz`.
- `BDO2019` stores the geometry as legacy `CYL` metadata with an implicit
  full-azimuth span (`phi = 0..360`, one bin) for compatibility.

## Differential scoring

A `Quantity` line can be followed by `Diff1`/`Diff1Type` (and optionally `Diff2`/`Diff2Type`)
to produce a differential (spectral) scorer.  The accumulator is expanded to
`geo_nbins × diff1_nbins` (single) or `geo_nbins × diff1_nbins × diff2_nbins` (double)
bins, matching the SH12A BDO data layout.

```text
Output
    Filename spectra.dat
    FileFormat TEXT
    Geo MyMesh
    Quantity Fluence
    Diff1     0.1  200.0  100  LOG   # lo hi nbins [LOG]
    Diff1Type EKIN                    # ekin | let | qeff | enuc | eamu
    Quantity Dose                     # plain dose — no differential axis
```

`Diff1` and `Diff2` each apply to the most recently parsed `Quantity` line.

| Sub-keyword | Required | Description |
|-------------|----------|-------------|
| `Diff1 lo hi nbins [LOG]` | no | Activates single-differential mode.  `LOG` selects logarithmic binning (requires `lo > 0`). |
| `Diff1Type <kind>` | no | Physical quantity for the first axis (default: `EKIN`). |
| `Diff2 lo hi nbins [LOG]` | no | Activates double-differential mode (requires `Diff1`). |
| `Diff2Type <kind>` | no | Physical quantity for the second axis (default: `EKIN`). |

Supported axis type keywords (same for Diff1Type and Diff2Type):

| Keyword | Synonyms | Axis quantity |
|---------|----------|---------------|
| `EKIN` | `E` | Kinetic energy [MeV] at step midpoint |
| `ENUC` | | Kinetic energy per nucleon [MeV/u] |
| `EAMU` | | Kinetic energy per atomic mass unit [MeV/u] |
| `LET` | `DEDX` | Electronic stopping power in transport medium [MeV/cm] |
| `QEFF` | `ZEFF2BETA2` | (z_eff/β)² |

Differential scoring is supported for `Fluence`, `Dose`, `DoseGy`, and `Energy`.
Averaged quantities (`DLET`, `TLET`, `DQEFF`, `TQEFF`) cannot carry a differential axis.

In `TEXT` output each diff-axis bin centre is written as an extra column between the
spatial coordinates and the quantity values.  For double-differential, Diff1 bins are the
outer (slow) loop and Diff2 bins are the inner (fast) loop:

```
# Diff1Type: EKIN  lo=0  hi=500  nbins=5
# Diff2Type: LET  lo=0.1  hi=100  nbins=5 LOG
# X Y Z EKIN LET FLUENCE
 0.0  0.0  10.0   50.0  0.200  0.000000e+00
 0.0  0.0  10.0   50.0  0.794  0.000000e+00
 0.0  0.0  10.0   50.0  3.162  3.654e-03
 ...  (next EKIN bin)
 0.0  0.0  10.0  150.0  0.200  0.000000e+00
 ...
```

## Scored quantities

| Keyword | Unit | Description |
|---------|------|-------------|
| `Dose` | MeV/g | Absorbed dose, SH12A-compatible |
| `DoseGy` | Gy | Absorbed dose in gray (`Dose * 1.602176634e-10`) |
| `Fluence` | 1/cm² | Particle fluence |
| `Energy` | MeV | Energy deposited in the voxel |
| `DLET` | MeV/cm | Dose-averaged LET |
| `TLET` | MeV/cm | Track-averaged LET |
| `DQEFF` | dim.less | Dose-averaged `(z_eff/beta)^2` |
| `TQEFF` | dim.less | Track-averaged `(z_eff/beta)^2` |

Quantities can be restricted to a material via `Settings`:

```text
Settings
    Name inWater
    Material Water

Output
    ...
    Quantity Dose
    Quantity Dose inWater
```

`Dose inWater` and `DoseGy inWater` score using the stopping power of water
regardless of the actual traversed material, equivalent to SH12A dose-to-water.
