# Material Module

This directory owns material-domain data structures and setup logic used by the
core library.

## Structure

- `osh_material.c`, `osh_material.h` — cold material workspace API
  (`struct osh_material_workspace`), material validation/finalization, and
  public material helpers such as dE/dx overrides.
- `osh_material_atomic_data.*` — element/isotope lookup data and derived atomic
  helpers used during material setup.
- `osh_material_icru.*` — built-in ICRU material definitions and expansion
  helpers.
- `runtime/osh_material_compile.*` — compilation from cold material definitions
  into runtime stopping-power and range tables.
- `runtime/osh_material_runtime.h` — runtime table layout consumed by transport
  and simulation.

## Lifecycle

```
osh_material_workspace_create()     allocate cold workspace
osh_material_workspace_prepare()    validate and derive composition fields
osh_material_compile()              compile runtime tables
osh_material_workspace_free()       release cold workspace
```

## Rules

- This directory does not own file-format parsing or path-based setup.
- It defines and validates the cold material model used by the rest of the
  library.
- Setup-time validation/assembly diagnostics from
  `osh_material_workspace_prepare()` are emitted through a borrowed
  diagnostics sink supplied by the caller.
- Runtime transport tables are built here, then consumed elsewhere.
- Compile-time diagnostics for runtime table construction are emitted through
  a borrowed diagnostics sink supplied by the caller of
  `osh_material_compile()`.
- Other modules may consume `struct osh_material_runtime`, but should not
  construct it directly.
