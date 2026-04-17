# Scoring Module

This module owns scoring-domain data and logic.

Structure:
- `src/scoring/` holds the public scoring API and cold workspace types.
- `src/apps/osh/` holds the current `detect.dat` parser (`osh_scoring_parse.*`) and
  path/file setup glue (`osh_scoring_setup_from_path`).
- `src/scoring/runtime/` holds the compiled runtime representation used during simulation.
- `src/scoring/save/` holds cold-path file writers consuming workspace plus runtime.

Rules:
- File-format parsing belongs to app code and fills the scoring workspace.
- Parsed scoring definitions stay separate from runtime scoring buffers.
- Save/output code should consume the scoring workspace plus scoring runtime, not transport internals.
- Other modules may call into scoring runtime APIs, but should not own scorer compilation.
- Scoring runtime compilation (`osh_scoring_prepare`), postprocessing, and saving
  are invoked by `src/simulation/` as part of `osh_simulation_create` and
  `osh_simulation_run`.  App code never calls them directly.
