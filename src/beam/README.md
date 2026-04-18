# Beam Module

This module owns beam-domain data structures and setup logic used by the core library.

## Structure

- `osh_beam.c/h` — cold beam workspace API (`struct osh_beam_workspace`) and public
  beam helpers, including `osh_beam_spots_set()`.
- `osh_beam_model.*` — beam model definitions (pencil, broad, spot-scanning).
- `osh_beam_spots.*` — spot-array allocation helpers and shared-value defaults.
- `osh_beamdef.h` — internal beam definition types shared within the module.
- `osh_beam_prepared.h` — opaque internal prepared state produced by
  `osh_beam_workspace_prepare()`.
- `runtime/osh_beam_runtime.*` — `osh_beam_compile()` and the hot runtime struct
  consumed by transport to sample primary particles.

## Lifecycle

```
osh_beam_workspace_create()         allocate cold workspace
osh_beam_workspace_prepare()        validate and derive beam model fields
osh_beam_compile()                  initialise primary source runtime
osh_beam_workspace_free()           release cold workspace
```

## Rules

- This directory owns cold beam data, beam validation, and beam runtime compilation.
- Setup-time validation and public beam summary output are emitted through
  borrowed diagnostics sinks supplied by the caller of
  `osh_beam_workspace_prepare()` and `osh_beam_print()`.
- Text-file loading and format-specific import happen outside this module.
- `osh_beam_compile()` is internal to `src/` and called only from `osh_simulation_create`.
- The runtime struct is owned by `struct osh_simulation`; nothing outside the
  simulation module holds a pointer to it.
