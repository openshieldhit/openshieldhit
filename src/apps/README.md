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
  `setup_from_path` functions, rewrites output filenames to full paths, then
  delegates to the library via explicit run/save calls:

  ```c
  osh_simulation_create(beam, geo, mat, scoring, diag, &sim);
  osh_simulation_run(sim);
  osh_simulation_save(sim);
  osh_simulation_free(sim);
  ```

  On the current branch, the app layer may also provide a borrowed
  `struct osh_diag_sink` to the simulation/runtime path.  This lets the app
  decide whether runtime diagnostics go to stderr, a file, a GUI widget, a
  test buffer, or nowhere at all.

  All runtime compilation, zone-to-material wiring, and transport execution are
  private to `src/simulation/`.  The app layer never touches runtime headers.

  Save cadence and execution policy still live here. The app decides whether to
  run serially or with multiple simulations in parallel, when to save partial
  results, and what constitutes one completed chunk of work. The library
  itself selects the concrete writer per output block based on the parsed
  scoring format.

  Diagnostics ownership is explicit: app, simulation, and transport use a
  caller-owned sink, and deeper modules should either do the same or stay
  silent unless they sit at a real setup/compile boundary.

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
