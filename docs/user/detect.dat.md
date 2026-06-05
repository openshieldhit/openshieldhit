# detect.dat reference

!!! note "Work in progress"
    This page is a stub.  Full documentation coming soon.

`detect.dat` defines scoring detectors: geometry (mesh, DICOM CT grid,
DICOM RTDOSE), quantities (dose, LET, fluence), and output formats.

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
