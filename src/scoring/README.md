# Scoring Module

This module owns scoring-domain data and logic.

## Structure

- `osh_scoring.c/h` — cold scoring workspace API (`struct osh_scoring_workspace`),
  public helpers for detector and filter definitions.
- `runtime/osh_scoring_compile.*` — compilation from cold scoring definitions into
  runtime accumulator buffers and geometry.
- `runtime/osh_scoring_runtime.h` — runtime layout consumed by transport and simulation.
- `runtime/osh_scoring_step.*` — per-step scoring accumulation called from transport.
- `runtime/osh_scoring_postprocess.*` — post-run normalisation and derived quantity
  computation.
- `save/` — output writers that consume workspace metadata plus runtime accumulators
  to produce scored result files.

## Lifecycle

```
osh_scoring_workspace_create()      allocate cold workspace
                                    (no workspace_prepare step — scoring has no derived fields)
osh_scoring_compile()               allocate accumulators, compile scoring geometry
osh_scoring_workspace_free()        release cold workspace
```

## Rules

- This directory owns cold scoring definitions, runtime compilation, postprocess,
  and output support.
- File-format loading happens outside this module.
- Parsed scoring definitions stay separate from runtime scoring buffers.
- Save/output code consumes the scoring workspace plus scoring runtime, not transport internals.
- Other modules may call into scoring runtime APIs, but should not own scorer compilation.
- `osh_scoring_compile`, `osh_scoring_postprocess`, and `osh_scoring_save` are
  invoked by `src/simulation/` as part of `osh_simulation_create`,
  `osh_simulation_run`, and `osh_simulation_save`. App code
  never calls them directly.
