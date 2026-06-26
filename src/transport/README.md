# Transport Module

This module owns transport-specific data flow and stepping orchestration.

Rules:
- `transport/` consumes runtime layers from other modules.
- `transport/` should not own preparation code for scoring, materials, or other domains.
- Shared simulation types that outgrow transport-specific ownership should move to a neutral shared header/module.

## Intended Design

The long-term transport driver is expected to own one queue or pool per
particle family:

- ions
- neutrons
- photons
- electrons

The outer loop should be scheduler-driven rather than hard-coded inside one
family kernel:

1. Seed the ion family with beam primaries.
2. Ask the scheduler which enabled family currently has work.
3. Drain that family's queue or pool with the matching transport kernel.
4. Let that kernel append secondaries into destination family queues.
5. Update queue state in the scheduler and repeat until every family is empty.

Current status:

- Ion transport is implemented as the main charged-particle CSDA wavefront.
- Neutron transport has a minimal pool-drain loop: it transports banked
  neutrons through geometry/material boundaries, samples free paths from the
  neutron macroscopic total cross section, and delegates reaction final states
  to `physics/neutron/`.
- Photon transport remains a stub.
- The scheduler seam exists now only to define family IDs, fixed-priority
  selection, and the dispatch ownership point in `osh_transport.c`.
- Cross-family feedback is still incomplete: neutron secondaries can be pushed
  back to the neutron pool, but charged secondaries from neutron reactions are
  not yet fed into the ion family and local neutron energy deposits are not yet
  scored.

Ownership notes:

- Family kernels own their step physics and local pool-drain logic.
- The scheduler owns cross-family orchestration policy.
- Material preparation, geometry compilation, and scoring setup stay outside
  `transport/`; transport should consume their runtime representations.

Diagnostics notes:

- The active runtime path receives a borrowed `struct osh_diag_sink` through
  `struct osh_transport_context`.
- Transport emits progress, unsupported-mode errors, and step-level runtime
  failures through that sink.
- `NULL` sink means silent transport diagnostics.
- Deeper modules outside `transport/` should either use their own explicit
  diagnostics sinks at real ownership boundaries or stay silent.

## Multithreading Direction

The intended parallelism model is not "one thread per particle family".
If ion transport dominates, most parallel work should happen inside the ion
family itself while the scheduler remains a lightweight top-level owner.

Expected shape:

1. The scheduler selects the next family with runnable work.
2. That family kernel submits many independent chunks to a shared worker pool.
3. Each worker drains its own local pool or scratch buffers to avoid hot
   contention on particle state.
4. Secondaries are buffered thread-locally by destination family.
5. Those thread-local buffers are merged back into the global family queues at
   explicit synchronization points.

Why this shape:

- It keeps worker threads busy even when one family dominates runtime.
- It avoids binding CPU resources to mostly idle families.
- It keeps cross-family orchestration separate from per-family inner-loop
  optimization.

Determinism notes for a future threaded implementation:

- RNG streams are keyed by global history index, not by position in one shared
  draw sequence.  Workers therefore own disjoint history ranges such as
  `[0, 250)` and `[250, 500)`, and each primary is seeded from
  `rndoffset + hist_lo + worker_local_index`.  This lets one history consume any
  number of random draws without shifting another worker's streams.
- Chunk ownership should be explicit and stable enough that runs can be debugged.
- Secondary merges should happen at clear boundaries rather than by fully
  shared push-on-every-step queues.
