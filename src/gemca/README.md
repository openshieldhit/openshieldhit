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

3. `src/gemca/prepare/`
   Core-internal preparation and compilation helpers.
   This layer turns the cold geometry model into the internal compatibility
   workspace `gemca_workspace` by:
   - allocating internal body/zone objects
   - building body surfaces and transforms
   - compiling raw zone expressions into pointer-linked CSG trees

4. `src/gemca/runtime/`
   The hot runtime representation used by transport.
   It compiles the internal compatibility workspace into flat arrays that are
   cache-friendly and SIMD-friendly.

The important boundary is:

`app parser -> cold geometry -> prepare/compile -> runtime`

Two consequences follow from that:

- File parsing belongs in `src/apps/`, not in GEMCA core.
- Expression compilation still belongs in GEMCA core, because library users may
  construct cold geometry programmatically and still need the same prepare step
  before runtime compilation.

The current `gemca_workspace` remains an internal compatibility type during the
migration. It is no longer the public geometry API.
