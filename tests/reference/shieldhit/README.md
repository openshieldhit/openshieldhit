# SH12A Reference Fixtures

This directory now serves two related purposes:

- small mirrored SH12A-compatible decks such as `00_minimal/` and
  `06_minimal_nucre/`
- curated SH12A gold-standard result fixtures for manual comparison against
  OpenShieldHIT

## MCS isolation decks (issue #133)

`idd_water_200mev_scat0/`, `scat1/`, `scat2/` are runnable SH12A-syntax mirrors
of the matching `tests/reference/idd_water_200mev_scat*` OpenShieldHIT cases:
200 MeV protons in water (R = 20 cm, L = 28 cm), NUCRE off + STRAGG off, differing
only in the MSCAT switch (0 = no scatter, 1 = Gaussian, 2 = Molière). They probe the
multiple-scattering distal-edge problem in isolation. `detect.dat` and `Water.txt`
are shared verbatim with the OSH cases; `beam.dat` (JPART0), `geo.dat` (zone→medium
map) and `mat.dat` (MEDIUM) carry the SH12A syntax.

Run with the installed `shieldhit`, then copy the produced `idd.dat` (its Page(0)
is the primary-proton fluence, which `compare_idd.py` reads directly) to
`../idd_water_200mev_scat<N>/reference/idd_sh12a.dat` to activate the CTest
`reference` comparison:

```bash
cd tests/reference/shieldhit/idd_water_200mev_scat2 && shieldhit .
cp idd.dat ../../idd_water_200mev_scat2/reference/idd_sh12a.dat
```

The curated fixture sets currently live in:

- `plan01_geoA_sobpcent/`
- `plan02_geoD_mono/`

Source repository:

- `https://github.com/APTG/2022_DCPT_LET/tree/main/data/sh12a`

Each fixture directory contains only the imported SH12A `*.dat` scorer outputs.
The original external SH12A inputs, plots, and scalar `NB_target*.txt` files are
not copied into `tests/reference/`.

Matching runnable OpenShieldHIT cases:

- `tests/cases/10_plan01_geoA_sobpcent`
- `tests/cases/11_plan02_geoD_mono`

## Page Map

The two curated fixture sets share the same scorer families and page numbering:

| File family | Page meaning |
|-------------|--------------|
| `NB_Z_narrow_dose_pN.dat` | `p1` Fluence; `p2` Dose; `p3` Dose `Protons`; `p4` Fluence `Primary`; `p5` Fluence `Protons` |
| `NB_Z_narrow_dose_water_pN.dat` | `p1` Fluence; `p2` Dose `in_Water`; `p3` Dose `Protons in_Water` |
| `NB_Z_narrow_LET_pN.dat` | `p1` DLET; `p2` DLET `Primary`; `p3` DLET `Protons`; `p4` TLET; `p5` TLET `Primary`; `p6` TLET `Protons` |
| `NB_Z_narrow_LET_water_pN.dat` | `p1` DLET `in_Water`; `p2` DLET `Primary in_Water`; `p3` DLET `Protons in_Water`; `p4` TLET `in_Water`; `p5` TLET `Primary in_Water`; `p6` TLET `Protons in_Water` |
| `NB_Z_narrow_QEFF_pN.dat` | `p1` DQEFF; `p2` DQEFF `Primary`; `p3` DQEFF `Protons`; `p4` TQEFF; `p5` TQEFF `Primary`; `p6` TQEFF `Protons` |
| `NB_target_diff_pN.dat` | `p1` Fluence vs `DEDX`; `p2` Fluence `Primary` vs `LET`; `p3` Fluence `in_Si` vs `DEDX`; `p4` Fluence `in_Si Primary` vs `DEDX` |
| `NB_target_water_diff_pN.dat` | `p1` Fluence `in_Water` vs `DEDX`; `p2` Fluence `Primary in_Water` vs `LET` |

## Completeness Of The Native Cases

The matching OpenShieldHIT cases intentionally emit more than the imported
reference fixture set at runtime:

- SH12A-style multi-page `*.bdo` files
- postprocessed page files such as `NB_Z_narrow_dose_p1.dat` after running
  `convertmc plotdata`

This keeps the native case definitions compact while still allowing the
generated page files to match the imported SH12A fixture naming.

## Known Comparison Caveat

The all-particle dose curves are currently expected to disagree visibly between
SH12A and OpenShieldHIT because secondary-particle production and transport are
not yet matching well enough. In practice this especially affects the manual
interpretation of pages such as:

- `NB_Z_narrow_dose_p2.dat`
- `NB_Z_narrow_dose_water_p2.dat`

For day-to-day development on this branch, it is usually more informative to
start with primary fluence, LET, QEFF, and the differential target scorers.

A residual configuration difference in some OpenShieldHIT example cases is
straggling: the original SH12A datasets used Vavilov straggling, and while
OpenShieldHIT now implements STRAGG 2 (Vavilov + Landau), several matching OSH
example inputs still default to STRAGG 1 (Gaussian) and can be switched to
STRAGG 2 for closer SH12A parity.

## Quick Start

1. Run OpenShieldHIT for one of the matching cases:

```bash
./build/bin/openshieldhit -v --outdir /tmp/plan01_geoA tests/cases/10_plan01_geoA_sobpcent
```

The curated scanned-field examples default to `NSTAT = 100000` primaries in
their `beam.dat` files.

2. Postprocess the generated BDO files into SH12A-style page files:

```bash
( cd /tmp/plan01_geoA && convertmc plotdata --many "*.bdo" )
```

3. Plot selected overlays against the imported SH12A fixtures:

```bash
python3 tools/plot_sh12a_reference.py \
    tests/reference/shieldhit/plan01_geoA_sobpcent \
    /tmp/plan01_geoA \
    NB_Z_narrow_dose_p4.dat NB_Z_narrow_LET_p1.dat NB_target_diff_p1.dat
```
