# gemca — Layering Notes

This directory now has three distinct responsibilities:

1. `include/openshieldhit/geometry.h`
   The public cold geometry model exposed by `osh_core`.
   Callers and app adapters fill flat `osh_geometry_body[]` and
   `osh_geometry_zone[]` arrays without knowing anything about GEMCA internals.

2. `src/apps/osh/osh_geometry_parse.c`
   The current OpenShieldHIT `geo.dat` file parser.
   This is app-side code: it opens files, understands input-card syntax, and
   fills the public cold geometry structs.
   `osh_geometry_setup_from_path()` (in `osh_app_osh`) wraps parse + prepare
   into one call, matching the pattern used by beam and material.

3. `src/gemca/prepare/`
   Core-internal preparation and compilation helpers.
   This layer turns the cold geometry model into the internal compatibility
   workspace `osh_gemca_prepared` by:
   - allocating internal body/zone objects
   - building body surfaces and transforms
   - compiling raw zone expressions into pointer-linked CSG trees

4. `src/gemca/runtime/`
   The hot runtime representation used by transport.
   It compiles the internal compatibility workspace into flat arrays that are
   cache-friendly and SIMD-friendly.

The important boundary is:

`app parser -> cold geometry -> prepare/compile -> runtime -> simulation`

Two consequences follow from that:

- File parsing belongs in `src/apps/`, not in GEMCA core.
- Expression compilation still belongs in GEMCA core, because library users may
  construct cold geometry programmatically and still need the same prepare step
  before runtime compilation.

The current `osh_gemca_prepared` remains an internal compatibility type during the
migration. It is no longer the public geometry API.

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
