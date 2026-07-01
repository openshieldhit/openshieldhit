TEST 10_plan01_geoA_sobpcent
2026-06-16 / curated from external SH12A reference data

Scanned proton field in the Solid Water / PMMA slab phantom used in the
external SH12A dataset `plan01_field01_geoA_SOBPcent`.

Purpose:

- keep a runnable OpenShieldHIT version of the SH12A setup in-tree
- preserve the differential scorers used on this branch
- provide an XZ neutron-fluence map for manual inspection of transported
  secondary neutrons
- emit SH12A-style multi-page BDO outputs that can be postprocessed with
  `convertmc`

Reference data:

- `tests/reference/shieldhit/plan01_geoA_sobpcent/`

Notes:

- `USECBEAM` is intentional here: this is a monoenergetic scanned field, not a
  single pencil spot.
- The original SH12A setup used Vavilov straggling; OpenShieldHIT currently
  falls back to Gaussian straggling here because Vavilov is not implemented.
- The example defaults to `NSTAT = 100000` primaries in `beam.dat`.
- The directory passed via `--outdir` is created automatically if it does not
  already exist.
- The imported SH12A gold data kept in-tree are the 32 `*.dat` files only.
- This case is registered under `tests/cases` without `args.cmake`, so CTest
  only dry-runs it by default. Full transport is meant for manual inspection.

## Quick Start

1. Run OpenShieldHIT:

```bash
./build/bin/openshieldhit -v --outdir /tmp/plan01_geoA tests/cases/10_plan01_geoA_sobpcent
```

The `/tmp/plan01_geoA` directory does not need to exist beforehand.

2. Postprocess the BDO outputs with `convertmc`:

```bash
( cd /tmp/plan01_geoA && convertmc plotdata --many "*.bdo" )
```

This should generate files such as `NB_Z_narrow_dose_p1.dat`,
`NB_Z_narrow_LET_p1.dat`, `NB_XZ_neutron_fluence_p1.dat`, and the
corresponding `NB_target*.dat` files.

3. Generate a multi-page PDF overlay report against the imported SH12A reference:

```bash
python3 tools/plot_sh12a_reference.py tests/reference/shieldhit/plan01_geoA_sobpcent /tmp/plan01_geoA
```

This writes `/tmp/plan01_geoA/plan01_geoA_sobpcent_overlay.pdf` with one page
per reference scorer and descriptive titles.

Differential plots in the report are shown on log-log axes by default for
easier manual inspection.

The neutron-fluence map scores transported neutron steps only.  Neutron
kerma and point-like local deposits are intentionally not included yet; those
need the shared point-deposit scoring API.
