## Bootstrap Runtime Pipeline: CLI + Geometry + Beam + Parser Tests

We should establish a minimal end-to-end runtime pipeline so `openshieldhit` can load core inputs and run in a dry/skeleton mode.

### Goals
1. Add a CLI parser with basic functionality.
2. Add geometry reader integration.
3. Add geometry reader tests.
4. Add beam reader integration (already implemented), so `main` can handle two workspaces.

### Proposed TODO
- [x] Implement CLI skeleton in `main`:
  - [x] `--help`
  - [x] `--version`
  - [x] input args for geometry and beam files (e.g. `--geo`, `--beam`)
  - [x] optional `--dry-run` (parse/load only, no transport yet)
- [x] Wire geometry workspace loading into runtime flow.
- [ ] Wire beam workspace loading into runtime flow.
- [x] Add clear startup summary output (which files loaded, basic counts/status).
- [ ] Add tests for geometry parser using SHIELD-HIT reference-like inputs:
  - [ ] valid geometry cases
  - [ ] malformed geometry cases (expected failure + useful line/file error)
- [ ] Add at least one integration test:
  - [ ] load geometry + beam together in dry-run mode and exit `0`.

### Acceptance Criteria
- `openshieldhit --help` and `openshieldhit --version` behave consistently across platforms.
- Running with valid geometry + beam inputs initializes both workspaces successfully.
- Invalid geometry inputs fail with actionable diagnostics (file + line where possible).
- CI includes geometry parser tests and a dry-run integration test.

### Notes
This issue is intentionally scoped to bootstrap the executable contract and input loading path; transport/scoring/material physics logic will be handled in follow-up issues.

---

## Mini-RFC: Shared-Library-First API

### Motivation
- Keep `main` thin and stable.
- Make core logic reusable for tests, tools, and future bindings.
- Prepare packaging for both CLI + shared library (`.so`/`.dylib`/`.dll`).

### High-level Architecture
- `openshieldhit` executable: CLI parser + user interaction only.
- `libopenshieldhit` shared library: public runtime API.
- Internal modules remain modular, but are consumed through one public facade.

### Public API (first draft)
- `openshieldhit_context_create()`
- `openshieldhit_context_destroy()`
- `openshieldhit_geometry_load()`
- `openshieldhit_beam_load()`
- `openshieldhit_run_dry()`
- `openshieldhit_last_error()`

Notes:
- Public API uses opaque handles (no exposed internal structs).
- Error handling should return codes and provide a retrievable error message.
- CLI should only call into this API.

### CMake / Packaging Plan
- Build and install shared library target: `openshieldhit`.
- Install public headers under `include/openshieldhit/`.
- Keep existing internal static libs for now; compose them into the shared facade target.

### Naming Policy (new files from now on)
- Public API files use `openshieldhit_*` naming.
- Internal source/header files in `src/` use `osh_*` naming.
- Keep `src/main.c` as the thin CLI entry point.
- Install only public API headers; internal headers remain non-installed.

### First Implementation Slice
- [ ] Add `src/api/` module for facade implementation.
- [ ] Add public header: `include/openshieldhit/openshieldhit.h`.
- [ ] Add minimal context object and create/destroy API.
- [ ] Add geometry + beam load wrappers calling existing workspace loaders.
- [ ] Refactor `src/main.c` to call API (no direct module wiring in CLI).
- [ ] Add one integration test for `--dry-run` using geometry + beam inputs.
- [ ] Keep internal implementation files in `osh_*` style; expose only `openshieldhit_*` in installed/public API.

### Acceptance Criteria
- CLI builds and runs by linking only to `libopenshieldhit`.
- A program can load geometry + beam through the public API without depending on internal headers.
- New files added in this work follow the no-`osh` filename rule.
