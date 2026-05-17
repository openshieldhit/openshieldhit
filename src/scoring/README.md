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
