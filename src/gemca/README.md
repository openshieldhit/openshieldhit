# GEMCA Module

This directory owns the geometry-compatibility layer that bridges the public
cold geometry model to the hot geometry runtime used by transport.

## Structure

- `osh_geometry_workspace.c` — geometry workspace preparation entry point.
- `osh_gemca2.*` — internal prepared geometry representation
  (`struct osh_gemca_prepared`) and related helpers.
- `prepare/` — conversion from cold geometry definitions into prepared GEMCA
  bodies, surfaces, transforms, and compiled CSG trees.
- `runtime/` — compilation from `osh_gemca_prepared` into flat runtime arrays
  used during transport.

## Lifecycle

```
osh_geometry_workspace_create()     allocate cold workspace
osh_geometry_workspace_prepare()    validate, build surfaces and CSG trees
osh_gemca_compile()                 compile runtime arrays
osh_geometry_workspace_free()       release cold workspace
```

## Rules

- This directory does not own file-format parsing or path-based setup.
- It owns the geometry preparation and compilation steps that turn cold
  geometry data into simulation-ready runtime arrays.
- `struct osh_gemca_prepared` is an internal compatibility type, not public API.

## Cold structs and compiled state coexist intentionally

After `osh_geometry_workspace_prepare()` succeeds, both representations live in
memory simultaneously:

- `osh_geometry_workspace` (cold) — the user-facing record: body names, types,
  raw argument lists, zone expressions, material name strings.
- `osh_gemca_prepared` (compiled) — derived state: surface lists, transformation
  matrices, compiled CSG ASTs, resolved material indices.

This is not redundancy to remove. The cold structs are the permanent,
stable description of the geometry as the user specified it. Post-runtime they
are the natural source for output headers, result annotations, and any
serialisation that needs to correlate scoring results back to zone or body
names. The compiled workspace is an ephemeral runtime artefact that could in
principle be rebuilt from the cold model at any time.
