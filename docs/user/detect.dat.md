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

### Multiple output formats — `FileFormat TEXT BDO`

One `Output` block can write its scored result in **several** formats at once by
listing more than one keyword on the `FileFormat` line:

```text
Output
    Filename NB_msh          # a stem when several formats are listed
    FileFormat TEXT BDO      # → NB_msh.dat  +  NB_msh.bdo
    Geo MyMesh
    Quantity Energy
    Quantity Fluence
```

This scores the geometry **once** and writes it out twice. It is strictly
cheaper than the old workaround of duplicating the whole `Output` block: the
scored arrays are shared, so **scoring memory does not grow with the number of
formats** (two formats or ten, the accumulator memory is the same), and each
history is deposited only once.

Filename rule:

- **A single format uses `Filename` verbatim** — every existing `detect.dat` is
  unaffected. `FileFormat TEXT` + `Filename NB_msh.dat` still writes exactly
  `NB_msh.dat`.
- **Several formats treat `Filename` as a stem.** A recognised trailing extension
  (`.dat .txt .bdo .bdz .bin .dcm .svg`) is stripped if present, then a canonical
  extension is appended per format: `.dat` (TEXT/ASCII/TXT/DAT), `.bdo`
  (BDO/BDO2019/BINARY/BIN), `.dcm` (RTDOSE), `.svg` (SVG). So both
  `Filename NB_msh` and `Filename NB_msh.dat` with `FileFormat TEXT BDO` yield
  `NB_msh.dat` + `NB_msh.bdo`.

#### Per-target filename overrides

By default the derived stem + canonical extension names every target. To give one
target a **verbatim** name instead — e.g. to keep a legacy filename that does not
match the canonical extension — write that name immediately after its format
keyword:

```text
Output
    Filename NB_msh
    FileFormat TEXT legacy.out BDO run.bdo   # → legacy.out  +  run.bdo
    Geo MyMesh
    Quantity Dose
```

Grammar and rules:

- A token that **is** a recognised format keyword (`TEXT/ASCII/TXT/DAT`,
  `BDO/BDO2019/BINARY/BIN`, `RTDOSE`, `SVG`) starts a new target.
- A token that is **not** a recognised keyword and immediately follows a format
  keyword is that target's override filename. It is used **verbatim** — stem
  stripping and canonical-extension derivation are skipped for that one target.
- Overrides can be mixed with derived names on the same line, and repeated
  `FileFormat` lines still accumulate.
- Collision detection runs after all names are resolved (derived or overridden):
  two targets resolving to the same path are rejected.
- **Ambiguity to be aware of:** because "not a known keyword ⇒ filename", a
  mistyped format (`FileFormat TEXT BDX`) is silently taken as an override
  filename (`BDX`) for the preceding format, not flagged as an unknown format.
  Likewise, a file whose name happens to be a bare format keyword (e.g. a file
  literally named `bdo`) cannot be expressed as an override. At most one override
  may follow each format keyword; a second trailing non-keyword token is an error.

Constraints:

- Two formats that resolve to the **same** file (e.g. `FileFormat TEXT DAT`, both
  `.dat`) are rejected with a clear error.
- `RTDOSE` **may** be combined with other formats in one block. The RTDOSE writer
  needs exactly one dose page, so the RTDOSE target consumes only the first
  `Dose`/`DoseGy` page of the shared page-set; the other pages exist for the
  TEXT/BDO/SVG targets and are **silently skipped** by the `.dcm` (do not expect
  `Fluence`/`LET` in it). A mixed block with **no** `Dose`/`DoseGy` page is
  rejected with a clear error.

Failure handling:

- Saving is **best-effort**: if one target cannot be written (e.g. its directory
  does not exist), the remaining targets are still written and the run's exit
  status reflects that a write failed. It is not all-or-nothing across targets.

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
- Supported quantities: `Energy`, `Fluence`, `Dose`, `DoseGy`, `DirtyDose`, `DirtyDoseGy`.
- `FileFormat BDO` (default) records which transport zone each bin is (`GEO_ZONES`
  tag) for labelling; `FileFormat TEXT` writes one row per zone with a numeric
  zone-index column. The per-zone volume is not stored in either — it is consumed by
  the ÷volume in postprocess, so the saved dose/fluence is already final.

### Native plot output — `FileFormat SVG`

> A **quick-look** plot for sanity-checking a result without leaving the
> terminal. The `openshieldhit` binary draws it itself — no Python, no plotting
> library — and the `.svg` opens in any web browser. For publication-quality
> figures, load the `.bdo`/`.dat` into
> [pymchelper](https://github.com/DataMedSci/pymchelper) instead.

Add `FileFormat SVG` to an `Output` and that output writes a 1-D line plot
instead of numeric data. To get the plot *and* the numbers from one block, list
both formats — `FileFormat SVG BDO` (or `SVG TEXT`) writes `bragg.svg` alongside
`bragg.bdo` from the same scored pages (see *Multiple output formats* above). Two
plot shapes are recognised automatically:

**Spatial profile** — a depth-dose / Bragg curve or any 1-D profile:

```text
Output
    Filename bragg          # ".svg" is appended automatically
    FileFormat SVG
    Geo MyMesh              # 1 × 1 × N mesh (one non-singleton axis)
    Quantity Dose
```

**Spectrum** — a differential page scored over a single spatial bin (a 0-D
voxel or a single `Zone`); x is the differential axis, log-scaled when the
binning is `LOG`:

```text
Output
    Filename spectrum
    FileFormat SVG
    Geo Spot               # 1 × 1 × 1 voxel, or a one-Zone geometry
    Quantity Fluence
    Diff1 10 600 40 LOG
    Diff1Type ENUC         # → x-axis "E/nucleon [MeV/u]", y "…/cm^2/(MeV/u)"
```

- Plotted values use the same normalisation as `TEXT` output, so a separate
  `FileFormat TEXT`/`BDO` `Output` for the same `Geo`/`Quantity` reports the
  same numbers.
- Each plotted `Quantity` is written to its **own** file. A single quantity keeps
  the `Filename` you gave (`bragg.svg`); with several, a `_p1`, `_p2`, … page
  suffix is added (`bragg_p1.svg`, `bragg_p2.svg`, …) so nothing is overwritten.
- The file is a plain SVG line plot — frame, grid, axis ticks, and one curve —
  that opens in any web browser.
- **Mixed outputs** are handled per page. When an `Output` holds several
  `Quantity` pages of which only some match the chosen shape — e.g. a 1-D depth
  mesh where one page adds a `Diff1` axis (making it 2-D spatial × energy) — only
  the matching pages are plotted and the rest are skipped.
- **Supported shapes:**
  - Spatial profile — `Geometry Mesh` with exactly one non-singleton axis, or
    `Geometry Cyl` with one non-singleton axis (R or Z).
  - Spectrum — a `Diff1` page over a single spatial bin (`1 × 1 × 1` mesh voxel
    or a one-`Zone` geometry).
- Anything that is not one of these 1-D shapes — 2-D/3-D meshes, 2-D
  `Diff1 × Diff2` spectra, spectra spanning several spatial bins, multi-zone
  profiles, rotated meshes — cannot be plotted, and that `Output` produces no
  file at all. Use `FileFormat TEXT`/`BDO` for that geometry instead, then
  visualise it with [pymchelper](https://github.com/DataMedSci/pymchelper),
  which handles 2-D maps, log axes, error bars, and multi-page output.

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
Averaged quantities (`DLET`, `TLET`, `DQEFF`, `TQEFF`, `DAVGE`, `TAVGE`, `DBETA`, `TBETA`)
cannot carry a differential axis.
Additive differential pages are reported as differential quantities: after
spatial postprocessing, each bin is divided by the width of its `Diff1` bin, and
by the product of `Diff1` and `Diff2` widths for double-differential pages.  For
example, `Quantity Fluence` with `Diff1Type EKIN` is reported per MeV
(`/cm²/MeV`, equivalent to `1/cm²/MeV`), not as fluence accumulated per energy
bin.  `LOG` binning uses the actual logarithmic bin widths for this
normalisation.

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
| `DirtyDose` | MeV/g | As `Dose`, but only charged particles above the mass-LET threshold (see below) |
| `DirtyDoseGy` | Gy | As `DoseGy`, but only charged particles above the mass-LET threshold |
| `Fluence` | 1/cm² | Particle fluence |
| `Energy` | MeV | Energy deposited in the voxel |
| `DLET` | MeV/cm | Dose-averaged LET |
| `TLET` | MeV/cm | Track-averaged LET |
| `DQEFF` | dim.less | Dose-averaged `(z_eff/beta)^2` |
| `TQEFF` | dim.less | Track-averaged `(z_eff/beta)^2` |
| `DAVGE` | MeV | Dose-averaged kinetic energy |
| `TAVGE` | MeV | Track-averaged kinetic energy |
| `DBETA` | dim.less | Dose-averaged relative speed `beta = v/c` |
| `TBETA` | dim.less | Track-averaged relative speed `beta = v/c` |

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

### Dirty dose

`DirtyDose` and `DirtyDoseGy` behave exactly like `Dose` and `DoseGy` but only
accumulate the contribution of **charged particles whose mass stopping power
(mass-LET) exceeds a fixed threshold** of `30 MeV·cm²/g` (equivalently `3 keV/µm`
in water at `ρ = 1 g/cm³`). The threshold is a compile-time constant
(`OSH_DIRTYDOSE_MASS_SP_THRESHOLD`). Neutral particles and any species without a
stopping-power table entry contribute nothing.

The mass-LET used for the gate follows the same medium as the dose: normally the
local (transport) material, or — when a `Settings` block overrides the scoring
medium — the override medium. So assuming a user-defined `inWater` settings override,
then `DirtyDose inWater` scores dose-to-water **and**
compares the mass-LET *in water* against the fixed threshold.
A density-only override does not change the mass-LET (Fano-invariant).

### Average kinetic energy and beta

`DAVGE`/`TAVGE` and `DBETA`/`TBETA` are dose- and track-averaged quantities built
the same way as `DLET`/`TLET` and `DQEFF`/`TQEFF`: each contributing step books a
per-step scalar — kinetic energy at the step midpoint, or the corresponding
relative speed `beta = v/c` — weighted by the energy deposited in the bin
(dose-averaged) or by the physical track length in the bin (track-averaged).

Unlike `DLET`/`TLET`/`DQEFF`/`TQEFF`, kinetic energy and beta are well-defined for
**any** particle, including neutrals and photons, so `DAVGE`/`TAVGE`/`DBETA`/`TBETA`
apply no charge gate and require no stopping-power table.

## Statistical uncertainty (error bars)

Per-bin Monte-Carlo **standard error** acquisition is turned off by default.
Turn it on **per estimator** by attaching a `Settings` block that carries
`Variance On` to the `Quantity` line — the same mechanism as `Quantity Dose inWater`:

```text
Settings
    Name withErr
    Variance On

Output
    Filename bragg.dat
    Geo MyDepthDoseCurve
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
`DQEFF`, `TQEFF`, `DAVGE`, `TAVGE`, `DBETA`, `TBETA`) the numerator and
denominator errors are combined in quadrature; ignoring their (strong, positive)
correlation makes this a slightly **conservative** over-estimate, never an
under-estimate.

**Cost & scope.** Enabling `Variance On` roughly doubles the scoring memory of the
affected pages (a companion sum-of-squares array per accumulator) and adds one
merge per batch boundary; it does not change the scored values beyond
floating-point summation order.  It is currently written by the `TEXT`/`ASCII`
writer only — the `BDO2019` writer still emits the values (its standard-error field
is a planned addition).

The batch count also comes from any active checkpoint cadence: a run with
`--dump-every-primaries` / `NSTAT n step` (`nsave`), or the `--score-replicas N`
diagnostic, uses those batches instead of the internal default.
