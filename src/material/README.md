# Material Module

This module owns material-domain data and logic.

Structure:
- `src/material/` holds the public material API and setup-facing workspace types.
- `src/material/runtime/` holds the compiled runtime tables used during simulation.

Rules:
- Parsing/setup data stays in the material workspace.
- Runtime stopping-power and range tables are built by `runtime/`.
- Other modules may consume material runtime structs, but should not own their construction.
