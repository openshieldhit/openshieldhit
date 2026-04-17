# Material Module

This module owns material-domain data and logic.

## Structure

- `src/material/` — public material API and cold workspace type `struct osh_material_workspace`.
- `src/material/runtime/` — compiled runtime tables used during simulation (`osh_material_prepare`).
- `src/apps/osh/` — `mat.dat` parser (`osh_material_parse.*`) and
  path/file setup glue (`osh_material_setup_from_path`).

## Lifecycle

```
osh_material_workspace_create()     allocate cold workspace
osh_material_parse()                fill from mat.dat  [app layer]
osh_material_workspace_finalize()   validate and derive composition fields
osh_material_prepare()              compile runtime tables  [osh_simulation_create]
osh_material_workspace_free()       release cold workspace
```

`osh_material_setup_from_path()` (app layer) wraps parse + finalize into one call.

## Rules

- File-format parsing belongs to app code; it fills `struct osh_material_workspace`.
- Core owns workspace allocation/finalization and runtime-table preparation.
- Other modules may consume `struct osh_material_runtime`, but should not construct it.
