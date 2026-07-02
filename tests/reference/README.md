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

The straggling isolation set (`idd_water_200mev_strag{0,1,2}`, issue #190) follows
the same recipe with **MSCAT off + NUCRE off**, toggling only the STRAGG switch.

### NUCRE isolation set (issue #212)

A family isolating the **nuclear-reaction channel** the same way: **MSCAT off +
STRAGG off**, 200 MeV proton pencil beam on the R = 20 cm, L = 28 cm water
phantom, NSTAT = 200000, fixed seeds, differing only in the NUCRE switch. With
scattering and straggling off, the nuclear channel is the only source of
secondaries, so the species decomposition and plateau spectrum isolate it.

| case                      | NUCRE                       | SH12A mirror |
|---------------------------|-----------------------------|--------------|
| `idd_water_200mev_nucre0` | 0 (off)                     | yes          |
| `idd_water_200mev_nucre1` | 1 (inelastic + pp-elastic)  | yes          |
| `idd_water_200mev_nucre2` | 2 (elastic only)            | no — OSH-only |
| `idd_water_200mev_nucre3` | 3 (inelastic only)          | no — OSH-only |

OpenShieldHIT's NUCRE takes 0–3; SHIELD-HIT12A only understands 0/1, so modes 2
and 3 (which decompose mode 1 into its elastic and inelastic parts) are
OpenShieldHIT-only diagnostics with no mirror deck.

`detect.dat` scores, on the 1 cm² × 0.5 mm central column, total `Dose` (the
auto-compared `idd.dat` col-4 headline) plus `Dose`/`Fluence` split by species —
primary protons, all protons, alphas, and heavy recoils (`Z > 2`; the C/O
recoils that p+A elastic will add) — and, in a thin mid-plateau slab, a
differential secondary spectrum dΦ/dE_kin vs E_kin (`spectrum.dat`, 0.1–300 MeV,
150 log bins).

Overlay OpenShieldHIT (run live, all four modes in parallel) against the
committed SH12A fixtures:

```bash
python3 tools/plot_nucre.py            # writes a multi-page nucre_report.pdf
```

The SH12A `idd.dat`/`spectrum.dat` fixtures are committed under
`shieldhit/idd_water_200mev_nucre{0,1}/` and are not re-run.  Note: OSH's
differential output is currently counts-per-bin, not a density, so the plot tool
divides by the log-bin width (SH12A already reports a density) — tracked as
issue #215.

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
