# Transport Module

This module owns transport-specific data flow and stepping orchestration.

Rules:
- `transport/` consumes runtime layers from other modules.
- `transport/` should not own preparation code for scoring, materials, or other domains.
- Shared simulation types that outgrow transport-specific ownership should move to a neutral shared header/module.
