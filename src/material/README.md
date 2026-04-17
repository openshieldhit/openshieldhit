# Material Module

This module owns material-domain data and logic.

Structure:
- `src/material/` holds the public material API and setup-facing workspace types.
- `src/material/runtime/` holds the compiled runtime tables used during simulation.
- `src/apps/osh/` holds the current `mat.dat` parser (`osh_material_parse.*`) and
  path/file setup glue (`osh_material_setup_from_path`).

Rules:
- File-format parsing belongs to app code and fills the material workspace.
- Material core owns workspace allocation/finalization (`create` + `finalize`)
  and transport-table preparation.
- Runtime stopping-power and range tables are built by `runtime/`.
- Other modules may consume material runtime structs, but should not own their construction.
