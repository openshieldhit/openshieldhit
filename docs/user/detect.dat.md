# detect.dat reference

!!! note "Work in progress"
    This page is a stub.  Full documentation coming soon.

`detect.dat` defines scoring detectors: geometry (mesh, cylindrical, DICOM CT grid,
DICOM RTDOSE), quantities (dose, LET, fluence), and output formats.

## Geometry types

### Mesh

Cartesian (X, Y, Z) scoring grid.

```
Geometry Mesh
    Name MyMesh
    X  -5.0   5.0   10    # lo  hi  nbins
    Y  -5.0   5.0   10
    Z   0.0  20.0  200
```

Each axis line gives the lower bound, upper bound, and number of equally-spaced
bins.  Axis order in the file does not matter.

### Cyl

Cylindrical (R, Z) scoring grid, rotationally symmetric around the Z axis.

```
Geometry Cyl
    Name MyCyl
    R   0.0   5.0    5    # lo  hi  nbins
    Z   0.0  20.0  200
```

`R` is the radial axis (inner radius to outer radius); `Z` is the axial axis
(along the beam direction).  Both axes must be present; order in the file does
not matter.  The axis origin is at (0, 0) in the universe frame unless a
rotation is specified.

Voxel volumes are computed exactly: V(ir) = π (r₁² − r₀²) Δz, so dose and
fluence are correctly normalised even for large radial bins.

## Scored quantities

| Keyword | Unit | Description |
|---------|------|-------------|
| `Dose` | MeV/g | Absorbed dose — no unit conversion, SH12A-compatible |
| `DoseGy` | Gy | Absorbed dose in gray (`Dose × 1.602176634 × 10⁻¹⁰`) |
| `Fluence` | 1/cm² | Particle fluence (ICRU definition) |
| `Energy` | MeV | Mean energy deposited per voxel |
| `DLET` | MeV/cm | Dose-averaged LET |
| `TLET` | MeV/cm | Track-averaged LET |
| `DQEFF` | dim.less | Dose-averaged (z_eff/β)² |
| `TQEFF` | dim.less | Track-averaged (z_eff/β)² |

Quantities can be restricted to a material via `Settings`:

```
Settings
    Name inWater
    Material Water

Output
    ...
    Quantity Dose
    Quantity Dose inWater
```

`Dose inWater` and `DoseGy inWater` score using the stopping power of water
regardless of the actual traversed material — equivalent to SH12A's dose-to-water.
