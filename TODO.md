## Runtime / API Status

The repo has moved past the original bootstrap phase. We now have:

- [x] Internal CLI parser in `src/cli/`
- [x] Thin `src/main.c` that configures a public API context and calls the facade
- [x] Public header `include/openshieldhit/openshieldhit.h`
- [x] Opaque `openshieldhit_context_t`
- [x] Context lifecycle, deep-copying setters, version helpers, and `openshieldhit_last_error()`
- [x] `OPENSHIELDHIT_RUN_VALIDATE` mode for parse/validate-only execution
- [x] Geometry loading wired into `openshieldhit_run()`
- [x] Startup / validation summary output
- [x] Basic tests for the public API (`tests/test_osh_main.c`)
- [x] Basic tests for the internal CLI parser (`tests/test_osh_cli.c`)

## Current High-Priority TODO

- [ ] Wire beam loading into `openshieldhit_run()`
  - Current state: beam path resolution and file-existence reporting are present, but the actual loader is still marked "loader wiring pending".
  - Note: the internal beam parser/workspace code already exists under `src/beam/`; it is just not connected to the public facade yet.
- [ ] Wire material loading into `openshieldhit_run()`
- [ ] Wire detect/scoring loading into `openshieldhit_run()`
- [ ] Implement `OPENSHIELDHIT_RUN_NORMAL`
  - Current state: `OPENSHIELDHIT_RUN_NORMAL` returns `OPENSHIELDHIT_STATUS_NOT_SUPPORTED`.
- [ ] Add a validate-mode integration test using `tests/res/test01/`
  - Goal: run geometry + beam + mat + detect resolution/loading through the public API or CLI and assert expected exit status.
- [ ] Add tests for `openshieldhit_last_error()`
  - Suggested coverage: missing geometry file, unsupported normal mode, and error reset on a subsequent successful call.
- [ ] Add tests for context setter ownership semantics
  - Suggested coverage: strings are deep-copied, `NULL` clears a field, and failed replacement does not clobber the previous value.

## Packaging / Install TODO

- [ ] Decide whether the public facade should ship as static-only, shared-only, or both
  - Current state: `openshieldhit_api` is a static library target.
- [ ] Install the public library target
- [ ] Install public headers under `include/openshieldhit/`
- [ ] Keep internal headers and the internal CLI parser non-installed

## Cleanup / Follow-up

- [ ] Improve validation diagnostics where possible
  - Especially for malformed geometry / beam input: include useful file and line information when the underlying parser can provide it.
- [ ] Revisit the public API once beam/material/detect loading are wired
  - Decide whether dedicated load/validate helpers are needed in addition to the current context + `openshieldhit_run()` model.

## Notes

- The earlier `src/api/` idea is no longer necessary; the facade currently lives in `src/openshieldhit.c`.
- The public API is now context-based and opaque; the CLI parser remains internal in `src/cli/`.
- The old `--dry-run` CLI concept maps to `OPENSHIELDHIT_RUN_VALIDATE` in the public API.
