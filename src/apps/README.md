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
  OpenShieldHIT-style parser/input handling plus the current `openshieldhit`
  executable `main()`.

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
