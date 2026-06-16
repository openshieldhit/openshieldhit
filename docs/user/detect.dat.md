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

A `Quantity` line can be followed by `Diff1` and `Diff1Type` to produce a
differential (spectral) scorer.  The accumulator is expanded to
`geo_nbins × diff_nbins` bins.

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

`Diff1` applies to the most recently parsed `Quantity` line.

| Sub-keyword | Required | Description |
|-------------|----------|-------------|
| `Diff1 lo hi nbins [LOG]` | no | Activates differential mode.  `LOG` selects logarithmic binning (requires `lo > 0`). |
| `Diff1Type <kind>` | no | Physical quantity for the axis (default: `EKIN`). |

Supported `Diff1Type` values:

| Keyword | Synonyms | Axis quantity |
|---------|----------|---------------|
| `EKIN` | `E` | Kinetic energy [MeV] at step midpoint |
| `ENUC` | | Kinetic energy per nucleon [MeV/u] |
| `EAMU` | | Kinetic energy per atomic mass unit [MeV/u] |
| `LET` | `DEDX` | Electronic stopping power in transport medium [MeV/cm] |
| `QEFF` | `ZEFF2BETA2` | (z_eff/β)² |

Differential scoring is supported for `Fluence`, `Dose`, `DoseGy`, and `Energy`.
Averaged quantities (`DLET`, `TLET`, `DQEFF`, `TQEFF`) cannot carry a differential
axis.

In `TEXT` output the diff-axis bin centre is written as an extra column between the
spatial coordinates and the quantity values:

```
# Z R EKIN FLUENCE
 5.000000e-01  5.000000e-01  1.954e+00  3.456789e-03
 5.000000e-01  5.000000e-01  2.389e+00  8.765432e-03
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
