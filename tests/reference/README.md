# Reference benchmarks

Comparison of openshieldhit results against external Monte Carlo codes
(SH12A now, FLUKA planned).  Unlike the `cases::*` integration tests (exit
codes and exact-seed output) and their shape-only `.dat` comparison, the
benchmarks here compare **absolute per-primary** quantities, so global dose
deficits are caught.

## IDD benchmarks (`idd_*`, CTest label `reference`)

Integral depth-dose curves for pencil beams in a water phantom:

| case               | beam            | phantom              | mesh           |
|--------------------|-----------------|----------------------|----------------|
| `idd_water_70mev`  | 70 MeV protons  | R = 20 cm, L = 8 cm  | 160 × 0.5 mm   |
| `idd_water_150mev` | 150 MeV protons | R = 20 cm, L = 18 cm | 180 × 1 mm     |
| `idd_water_230mev` | 230 MeV protons | R = 20 cm, L = 36 cm | 360 × 1 mm     |

All use NUCRE 1, Gaussian straggling, Molière MCS, LOADDEDX water stopping
power, NSTAT = 200000, fixed seeds.

### MCS isolation set (issue #133)

A second family probes the multiple-scattering distal-edge problem (issue #133)
in isolation: **NUCRE off and STRAGG off**, so the distal edge is governed by MCS
alone. All three share a 200 MeV pencil beam on a water phantom (R = 20 cm,
L = 28 cm), NSTAT = 200000, fixed seeds, and differ only in the MSCAT switch:

| case                     | MSCAT          |
|--------------------------|----------------|
| `idd_water_200mev_scat0` | 0 (no scatter) |
| `idd_water_200mev_scat1` | 1 (Gaussian)   |
| `idd_water_200mev_scat2` | 2 (Molière)    |

Unlike the IDD cases above, the auto-compared `idd.dat` here is the **primary-proton
fluence** on a **1 cm² central column** (0.5 mm bins, *not* laterally integrated), so
it is sensitive to lateral escape; each case also writes a full-width laterally
integrated `ddc_wide.dat` and a radial `rad_cyl.dat` for manual comparison.

### Workflow

1. Reproduce the case setup in the reference code (see the per-case README;
   beam/geo/mat decks are SH12A-compatible as-is, only the scoring deck needs
   translating).
2. Export the laterally integrated depth-energy curve **per primary** as
   text: two columns `z [cm]  energy/bin` or openshieldhit mesh format
   (`X Y Z E ...`).  Binning may differ from the test mesh — the comparator
   converts to energy density and interpolates.
3. Drop it into `<case>/reference/idd_sh12a.dat` (any `reference/*.dat` is
   picked up; one test compares against every curve present).
4. Run: `ctest -L reference` — cases without reference data are SKIPPED.

### Pass criteria (`compare_idd.py`)

1. integral ratio within `--tol-integral` (default 5 %),
2. every bin within `--tol-bin` of the reference maximum (default 5 %,
   global-normalised so the distal falloff does not blow up relative errors),
3. peak position within `--tol-peak-mm` (default 1 mm; 95 %-of-max centroid).

Override per case in `<case>/compare_args.cmake`:

```cmake
set(COMPARE_ARGS --tol-integral 0.03 --tol-bin 0.03 --tol-peak-mm 0.5)
```

Start generous while physics gaps are known (untransported neutrons are
worth a few % of the integral at 150+ MeV), then tighten as they close to
pin against regressions.

## `shieldhit/`, `topas/`

`shieldhit/` contains both mirrored SH12A-compatible input decks and curated
SH12A gold-standard result fixtures for manual comparison. `topas/` currently
holds mirrored input decks only. Neither subtree is registered as a CTest
suite directly.
