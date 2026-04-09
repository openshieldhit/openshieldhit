# Scoring Module

This module owns scoring-domain data and logic.

Structure:
- `src/scoring/` holds the public scoring API and cold workspace types.
- `src/scoring/parse/` holds `detect.dat` parsing.
- `src/scoring/runtime/` holds the compiled runtime representation used during simulation.
- `src/scoring/save/` holds cold-path file writers consuming workspace plus runtime.

Rules:
- Parsed scoring definitions stay separate from runtime scoring buffers.
- Save/output code should consume the scoring workspace plus scoring runtime, not transport internals.
- Other modules may call into scoring runtime APIs, but should not own scorer compilation.
