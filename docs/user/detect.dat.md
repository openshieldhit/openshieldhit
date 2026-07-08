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

### Zone

Score by `geo.dat` zone, selected **by zone name** — one output bin per listed
zone. Internal numeric zone indices are never used in input files.

```text
Geometry Zone
    Name MyZones
    Zone WaterBox
    Volume 13.37
    Zone Target
    Volume 2.0
```

Rules and semantics:

- `Zone <name>` selects one zone by its `geo.dat` name; repeat it for each zone to
  score. Unknown names are a hard error at setup.
- `Volume <cm3>` optionally follows a `Zone` card and sets the volume for that
  zone's volume-normalized quantities (`Dose`, `DoseGy`). Zone volumes cannot be
  computed from overlapping CSG primitives, so they must be given explicitly; a
  zone with no `Volume` warns and defaults to `1.0 cm3`.
- The output bin order is exactly the order of the `Zone` lines.
- Supported quantities: `Energy`, `Fluence`, `Dose`, `DoseGy`.
- `FileFormat BDO` (default) records which transport zone each bin is (`GEO_ZONES`
  tag) for labelling; `FileFormat TEXT` writes one row per zone with a numeric
  zone-index column. The per-zone volume is not stored in either — it is consumed by
  the ÷volume in postprocess, so the saved dose/fluence is already final.

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
| `Diff1Type <kind> [settings]` | no | Physical quantity for the first axis (default: `EKIN`).  Optional `settings` name overrides the stopping-power medium and/or density used for LET/QEFF axis binning (see below). |
| `Diff2 lo hi nbins [LOG]` | no | Activates double-differential mode (requires `Diff1`). |
| `Diff2Type <kind> [settings]` | no | Physical quantity for the second axis (default: `EKIN`).  Same optional `settings` override as `Diff1Type`. |

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
Additive differential pages are reported as differential quantities: after
spatial postprocessing, each bin is divided by the width of its `Diff1` bin, and
by the product of `Diff1` and `Diff2` widths for double-differential pages.  For
example, `Quantity Fluence` with `Diff1Type EKIN` is reported per MeV
(`1/cm²/MeV`), not as fluence accumulated per energy bin.  `LOG` binning uses
the actual logarithmic bin widths for this normalisation.

### Axis stopping-power override

For `LET` and `QEFF` axis types, the stopping power used for bin determination is normally
that of the transport medium at each step. An optional `Settings` reference on the
`Diff1Type`/`Diff2Type` line can override the evaluation medium, the density, or both:

```text
Settings
    Name in_Si
    Material Si

Output
    Filename spectra_Si.dat
    FileFormat TEXT
    Geo MyMesh
    Quantity Fluence
    Diff1     0  2000  1000
    Diff1Type DEDX in_Si       # bin by Si stopping power, not transport-medium SP
```

This is the correct way to compute dΦ/dLET_Si (fluence spectrum in silicon LET).
The override applies only to the axis binning; the fluence accumulation itself remains
material-independent.

A density-only override is also valid:

```text
Settings
    Name rho_half
    Density 0.5

Output
    Quantity Fluence
    Diff1     0  2000  1000
    Diff1Type DEDX rho_half    # bin by 0.5 g/cm3 * S_transport(E)
```

The same syntax works for `Dose` + `LET` axis, where the dose is converted to
dose-to-water while the axis is also binned in water LET:

```text
    Quantity Dose in_Water
    Diff1     0  2000  1000
    Diff1Type DEDX in_Water    # Quantity-level and Diff1Type overrides are independent
```

When a `Settings` override is active the ASCII output header reads:
`# Diff1Type: LET  lo=…  hi=…  nbins=…  (SP override active)`

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

## Statistical uncertainty (error bars)

Per-bin Monte-Carlo **standard error** is off by default.  Turn it on **per
estimator** by attaching a `Settings` block that carries `Variance On` to the
`Quantity` line — the same mechanism as `Quantity Dose inWater`:

```text
Settings
    Name withErr
    Variance On

Output
    Filename bragg.dat
    Geo Depth
    Quantity Dose withErr
```

Only the pages that reference such a block gain error columns; every other page
is scored exactly as before.  (There is deliberately no run-wide card or flag for
this yet — enabling it is a per-page decision.)

**How it works.** The run is split into contiguous, equal blocks of primary
histories (batches).  Each batch is scored independently and treated as one
observation of the per-primary mean; the spread *between* batches is the
uncertainty (the *batch-means* method).  A run therefore needs at least two
batches to report anything — a run with no batching has zero degrees of freedom
and its error columns are all zeros.  A plain single-threaded run with at least
one variance-tracking page is automatically split into a default number of
batches so it produces error bars out of the box.  Because each history's random
stream is a pure function of its global index, a fixed batch count is
bit-reproducible run to run.

**Output.** In `TEXT`/`ASCII` output each quantity gains a paired error column
immediately after its value column, e.g.

```text
# X Y Z DOSE DOSE_ERR DLET DLET_ERR
```

The error is the standard error of that cell's reported value, in the **same
units** (it already carries the per-primary or physical-mean scaling), so a plot
can use it directly as a `± ` bar.  For the averaged quantities (`DLET`, `TLET`,
`DQEFF`, `TQEFF`) the numerator and denominator errors are combined in
quadrature; ignoring their (strong, positive) correlation makes this a slightly
**conservative** over-estimate, never an under-estimate.

**Cost & scope.** Enabling `Variance On` roughly doubles the scoring memory of the
affected pages (a companion sum-of-squares array per accumulator) and adds one
merge per batch boundary; it does not change the scored values beyond
floating-point summation order.  It is currently written by the `TEXT`/`ASCII`
writer only — the `BDO2019` writer still emits the values (its standard-error field
is a planned addition).

The batch count also comes from any active checkpoint cadence: a run with
`--dump-every-primaries` / `NSTAT n step` (`nsave`), or the `--score-replicas N`
diagnostic, uses those batches instead of the internal default.
