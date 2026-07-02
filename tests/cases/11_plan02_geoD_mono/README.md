TEST 11_plan02_geoD_mono
2026-06-16 / curated from external SH12A reference data

Scanned proton field in the Solid Water / PMMA slab phantom used in the
external SH12A dataset `plan02_field01_geoD_mono`.

Purpose:

- keep a runnable OpenShieldHIT version of the SH12A setup in-tree
- preserve the differential scorers used in this case
- emit SH12A-style multi-page BDO outputs that can be postprocessed with
  `convertmc`

Reference data:

- `tests/reference/shieldhit/plan02_geoD_mono/`

Notes:

- `mono` refers to the energy-layer setup, not to a single spot. `USECBEAM` is
  therefore part of the intended test coverage.
- The original SH12A setup used Vavilov straggling. This case is currently
  configured with `STRAGG 1` (Gaussian) and can now be switched to
  `STRAGG 2` for Vavilov/Landau straggling.
- The example defaults to `NSTAT = 100000` primaries in `beam.dat`.
- The directory passed via `--outdir` is created automatically if it does not
  already exist.
- The imported SH12A gold data kept in-tree are the 32 `*.dat` files only.
- This case is registered under `tests/cases` without `args.cmake`, so CTest
  only dry-runs it by default. Full transport is meant for manual inspection.

## Quick Start

1. Run OpenShieldHIT:

```bash
./build/bin/openshieldhit -v --outdir /tmp/plan02_geoD tests/cases/11_plan02_geoD_mono
```

The `/tmp/plan02_geoD` directory does not need to exist beforehand.

2. Postprocess the BDO outputs with `convertmc`:

```bash
( cd /tmp/plan02_geoD && convertmc plotdata --many "*.bdo" )
```

This should generate files such as `NB_Z_narrow_dose_p1.dat`,
`NB_Z_narrow_LET_p1.dat`, and the corresponding `NB_target*.dat` files.

3. Generate a multi-page PDF overlay report against the imported SH12A reference:

```bash
python3 tools/plot_sh12a_reference.py tests/reference/shieldhit/plan02_geoD_mono /tmp/plan02_geoD
```

This writes `/tmp/plan02_geoD/plan02_geoD_mono_overlay.pdf` with one page per
reference scorer and descriptive titles.

Differential plots in the report are shown on log-log axes by default for
easier manual inspection.
