# Simulation Module

This module owns the orchestration layer that compiles four cold workspaces
into a ready-to-run simulation and drives it to completion.

## Responsibility

`src/simulation/` sits between the public workspace API and the transport kernel.
It is the only place in the library that holds references to all four runtime
representations at once:

| Cold workspace              | Runtime compiled here             |
|-----------------------------|-----------------------------------|
| `osh_beam_workspace`        | `osh_beam_runtime`                |
| `osh_geometry_workspace`    | `osh_gemca_runtime`               |
| `osh_material_workspace`    | `osh_material_runtime`            |
| `osh_scoring_workspace`     | `osh_scoring_runtime`             |

All four runtime structs are private fields of the opaque `struct osh_simulation`.
Nothing outside this module ever sees them.

## Public API

```c
/* include/openshieldhit/simulation.h */

enum osh_status osh_simulation_create(
    struct osh_beam_workspace     *beam,
    struct osh_geometry_workspace *geo,
    struct osh_material_workspace *mat,
    struct osh_scoring_workspace  *scoring,
    struct osh_diag_sink const    *diag,
    struct osh_simulation        **sim_out);

enum osh_status osh_simulation_run(struct osh_simulation *sim);

enum osh_status osh_simulation_get_results(
    struct osh_simulation const *sim,
    struct osh_results const   **out);

unsigned long long osh_results_requested_nstat(struct osh_results const *results);
unsigned long long osh_results_completed_nstat(struct osh_results const *results);
int osh_results_has_completed_run(struct osh_results const *results);

enum osh_status osh_simulation_save(struct osh_simulation const *sim);

enum osh_status osh_simulation_free(struct osh_simulation *sim);
```

The cold workspaces are **borrowed**: the caller creates and owns them; the
simulation holds non-owning pointers.  `osh_simulation_free` releases only the
runtime resources compiled by `osh_simulation_create`.

The diagnostics sink is also borrowed.  On the current branch, this explicit
sink is the preferred path for runtime-facing diagnostics from simulation and
transport.  Passing `NULL` makes that path silent.

## What `osh_simulation_create` does

In order:

1. **Zone → material resolution** — walks `geo->prepared->zones[]` and resolves
   each zone's `material_name` string to a dense material index via
   `osh_material_by_name`.  This is the cross-domain wiring step that links
   geometry to materials.

2. **Geometry runtime** — calls `osh_gemca_compile` to compile the pointer-
   linked CSG trees from `osh_gemca_prepared` into flat, cache-friendly arrays
   (`osh_gemca_runtime`).

3. **Transport tables** — calls `osh_material_compile` to build stopping-power and
   CSDA-range tables for all material/projectile combinations up to the beam's
   maximum Z.

4. **Scoring runtime** — calls `osh_scoring_compile` to allocate accumulators and
   compile the scoring geometry definitions.

5. **Transport parameters** — copies scalar knobs from the beam workspace into
   `osh_transport_params` and translates beam-specific scatter/straggling enums
   into transport-owned enums.

6. **Beam runtime** — calls `osh_beam_compile` to initialise the primary
   source machinery.

The borrowed diagnostics sink is stored on the simulation and forwarded into
the transport context.  This keeps the main runtime path free of process-global
logger ownership.

## What `osh_simulation_run` does

1. **Transport** — drives `osh_transport_run_minimal` to completion.
2. **Postprocess** — calls `osh_scoring_postprocess` to finalise accumulators.

After a successful run, the simulation records both:

- the requested number of primaries
- the completed number of primaries actually represented by the current results

Today those values are typically equal, but keeping both counters distinct
prepares the API for future chunked, time-limited, or partially completed runs.

Saving is a separate explicit step. `osh_simulation_save` iterates the parsed
scoring outputs and selects the concrete writer for each output block based on
its configured format. The output paths themselves are expected to be fully
resolved by the application before save is called.

Calling `osh_simulation_save` before a completed run is invalid and returns
`OSH_ESTATE`.

## What Simulation Does Not Decide

This module deliberately does not own save cadence or chunk orchestration.
Those choices belong to the application layer.

Examples of app-owned policy:

1. Save every completed chunk of primaries.
2. Save every N seconds of wall-clock time.
3. Run one simulation serially or several simulations in parallel threads.
4. Merge partial results before saving.

The simulation layer exposes explicit `run` and `save_*` entry points, but it
does not currently prescribe when save calls should happen relative to run
beyond the obvious requirement that saving meaningful results implies a
completed run or run chunk. Format dispatch itself is owned by the library.

## Design rationale

### Why an opaque handle?

The four runtime structs are implementation details.  Exposing them would force
every caller to include internal headers from `beam/runtime/`, `gemca/runtime/`,
`material/runtime/`, and `scoring/runtime/`.  Changes to any of those structs
would break the calling code even when the observable behaviour does not change.

The opaque handle also makes it natural to change the internal layout — for
example, to add thread-local scratch buffers, a secondary particle queue, or an
event callback — without touching any code outside this module.

### Why does simulation start from cold workspaces?

This module assumes the four cold workspaces already exist and have been
prepared where needed. `src/simulation/` does not load files or interpret input
formats; it starts at the boundary where cold domain objects are available and
compiles them into runtime state.

### Lifetime rules

```
osh_X_workspace_create()   ← caller
osh_X_workspace_prepare()  ← caller
                           ← cold workspace alive
osh_simulation_create()    ← compiles runtimes, borrows cold workspaces
osh_simulation_run()       ← transport + scoring postprocess
osh_simulation_get_results() ← borrowed read-only handle to compiled results
osh_results_*_nstat()      ← inspect requested vs completed primaries
osh_simulation_save()      ← writes all configured output formats
osh_simulation_free()      ← frees runtime state only
osh_X_workspace_free()     ← caller
```

The cold workspaces must remain alive for the full lifetime of the simulation
object because result saving reads workspace metadata to annotate output files.

The same borrowed-lifetime rule applies to the diagnostics sink passed to
`osh_simulation_create()`.

## Transitional Note

The diagnostics sink currently covers the main runtime path:

- `osh_simulation_create`
- `osh_simulation_run`
- `osh_simulation_save`
- `transport/`

Older modules that have not yet been migrated may still emit through the
legacy global logger API.  That coexistence is intentional for now to keep the
migration branch narrow.
