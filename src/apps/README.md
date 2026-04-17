# Apps

This directory contains executable-facing application code built on top of
`osh_core`.

## Intent

- `include/openshieldhit/...` is the public library API.
- `src/...` contains the core implementation behind that API.
- `src/apps/...` contains format-specific parsers, CLI entry points, and other
  application layers that use the core library.

Application code here is private to this repository. It is not part of the
installed public API and should not be re-exported through `include/`.

## Current layout

- `src/apps/osh/`
  OpenShieldHIT-style parser/input handling plus the `openshieldhit` executable
  `main()`.

  The app layer provides one `setup_from_path` entry point per input file,
  each with a consistent signature:

  ```c
  enum osh_status osh_beam_setup_from_path    (path, lg, &beam_ws);
  enum osh_status osh_geometry_setup_from_path(path, lg, &geo_ws);
  enum osh_status osh_material_setup_from_path(path, lg, &mat_ws);
  enum osh_status osh_scoring_setup_from_path (path, lg, &scoring_ws);
  ```

  Each function opens a file, parses it, and prepares the cold workspace.
  Ownership stays with the caller; workspaces are freed independently via the
  corresponding `osh_X_workspace_free()` from the core library.

  `osh_run.c` orchestrates the top-level run: it resolves paths, calls the four
  `setup_from_path` functions, then delegates to the library via three calls:

  ```c
  osh_simulation_create(beam, geo, mat, scoring, &sim);
  osh_simulation_run(sim, out_dir);
  osh_simulation_free(sim);
  ```

  All runtime compilation, zone-to-material wiring, and transport execution are
  private to `src/simulation/`.  The app layer never touches runtime headers.

## Expected future layout

- `src/apps/osh_topas/`
- `src/apps/osh_fluka/`
- `src/apps/osh_json/`
- `src/apps/osh_webassembly/`

Each app directory may contain:

- parsers or importers for its source format
- glue code that fills the public core structs
- a `main()` when the app is a standalone executable

## Design rule

`osh_core` must not depend on any app.

Apps may depend on:

- the public library API in `include/openshieldhit/...`
- internal helpers from `src/...` when needed during migration

Over time, app code should rely on public helpers where that makes architectural
sense, but we do not need to force that separation all at once.
