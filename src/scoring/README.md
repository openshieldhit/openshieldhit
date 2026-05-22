# Scoring Module

This module owns scoring-domain data and logic.

## Structure

- `osh_scoring.c/h` — cold scoring workspace API (`struct osh_scoring_workspace`),
  public helpers for detector and filter definitions.
- `runtime/osh_scoring_compile.*` — compilation from cold scoring definitions into
  runtime accumulator buffers and geometry.  For geometries with
  `rtdose_template_path` set (set by the app when the output format is RTDOSE),
  the path is carried through to the geometry runtime so the save layer can
  write back to the template file.
- `runtime/osh_scoring_runtime.h` — runtime layout consumed by transport and simulation.
- `runtime/osh_scoring_step.*` — per-step scoring accumulation called from transport.
- `runtime/osh_scoring_postprocess.*` — post-run normalisation and derived quantity
  computation.
- `save/` — output writers; dispatched by `osh_scoring_save()` based on the
  `fileformat` field of each output definition:
  - `osh_scoring_save_ascii.*` — plain-text column format (`text`, `txt`, `ascii`, `dat`)
  - `osh_scoring_save_bdo2019.*` / `osh_scoring_save_bdo2019_raw.*` — binary BDO format
    (`bdo`, `bdo2019`, `binary`, `bin`); default when no format is specified
  - `osh_scoring_save_rtdose.*` — DICOM RTDOSE round-trip writer (`rtdose`): reads
    the RTDOSE template stored in `rtdose_template_path`, replaces the pixel data
    with scored values normalised by `nstat` and `dose_grid_scaling`, and writes
    the result as a new `.dcm` file.  All DICOM metadata is preserved unchanged.
    Only single-page (single quantity) outputs are supported; multi-page returns
    `OSH_ENOTSUP`.

## Output normalisation and multi-run merging

The two output formats take deliberately different positions on normalisation:

### ASCII (`text`, `txt`, `ascii`, `dat`)

Intended for quick inspection and single-run use.  Values are normalised
**per primary particle** (divided by `nstat`) at write time, except for
averaged quantities (DLET, TLET) which are already physical means after
`osh_scoring_postprocess()` and are written as-is.

ASCII output is **not** suitable for merging partial results from multiple
runs: once divided by `nstat` the absolute weight of each run is lost.

### BDO 2019 (`bdo`, `bdo2019`, `binary`, `bin`; default)

Intended for production use and multi-run accumulation.  Values are written
**exactly as accumulated** (raw sums), and the primary count is embedded as
the `OSHBDO_RT_NSTAT` tag.  The reader is responsible for normalisation.

This design supports two merging strategies:

1. **Native multi-threaded** (future): each thread accumulates into its own
   page buffers, then the simulation layer sums the per-thread pages and
   writes one BDO file with the total `nstat`.

2. **Embarrassingly-parallel / multi-node**: the user runs independent
   simulations and merges the resulting `.bdo` files.  Because each file
   carries its own `nstat`, the merge tool can apply the correct weighting
   per scorer kind:
   - `NORM` quantities (DOSE, FLUENCE, ENERGY, …): weighted average
     `X = (sum_j x_j) / (sum_j nstat_j)`
   - `AVER` quantities (DLET, TLET): weighted average
     `X = (sum_j x_j * nstat_j) / (sum_j nstat_j)`
   - `SUM` quantities (COUNT, …): plain sum `X = sum_j x_j`
   - `APPEND` quantities (MCPL): concatenation

   The `OSHBDO_PAG_NORMALIZE` tag in each page records the postproc mode so
   a merge tool does not need to re-derive it from the scorer type.

### Unit handling

`osh_scoring_postprocess()` applies physics transforms before saving:

- DOSE: converts MeV/g → Gy (`× OSH_MEVG2GY`) once per bin.
- DLET, TLET: computes the final ratio `data / data2` per bin; the
  intermediate `data2` denominator array is not written to any output file.

After postprocessing, page buffers hold values in their final physical units
(Gy, MeV, 1/cm², MeV/cm, …) but are still un-normalised raw sums.  Both
ASCII and BDO apply further normalisation as described above.

## Lifecycle

```
osh_scoring_workspace_create()      allocate cold workspace
osh_scoring_compile()               allocate accumulators, compile scoring geometry
osh_scoring_workspace_free()        release cold workspace
```

## Rules

- This directory owns cold scoring definitions, runtime compilation, postprocess,
  and output support.
- File-format loading and path resolution happen outside this module (app layer).
- Parsed scoring definitions stay separate from runtime scoring buffers.
- Save/output code consumes the scoring workspace plus scoring runtime, not
  transport internals.
- `osh_scoring_compile`, `osh_scoring_postprocess`, and `osh_scoring_save` are
  invoked by `src/simulation/` as part of `osh_simulation_create`,
  `osh_simulation_run`, and `osh_simulation_save`. App code never calls them
  directly.
