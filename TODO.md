## Current State

The codebase now has a clear split between:

- `src/apps/` — file-format parsing, path resolution, and executable-facing glue
- `src/` core modules — cold domain models, prepare/compile steps, runtime data,
  transport, scoring, and simulation orchestration

Recent architectural cleanups already completed:

- [x] `LOADDEDX` import moved from core material code to app-side import
- [x] spotlist import moved from core beam code to app-side import
- [x] cold-to-runtime steps renamed from `prepare` to `compile` where they
      produce separate runtime objects
- [x] simulation orchestration split into `src/simulation/`
- [x] module READMEs updated to match the current layering
- [x] explicit diagnostics sink replaces the old library-owned global logger
      on the main runtime/setup path
- [x] legacy logger API removed; internal diagnostics helper renamed to
      `osh_diag`

## Near-Term API / Architecture

- [ ] Revisit remaining path/file metadata in public cold workspaces
  - `wdir`, `fname`, and similar fields should ideally disappear from public
    API as the remaining app-owned path handling is cleaned up
- [ ] Add result merge API for embarrassingly parallel runs
  - Likely shape: merge accumulated scoring/results state from multiple
    simulation instances before save/postprocess
- [ ] Add chunked / partial run control
  - Let the app run a simulation in batches, inspect partial results, merge,
    and save explicitly
- [ ] Finish naming consistency passes where worthwhile
  - remaining stale "prepare" wording in messages/comments
  - remaining internal compatibility aliases where they no longer help

## Beam

- [ ] Add ridge-modulator / ripple-filter support
- [ ] Add MCPL phase-space import
- [ ] Add EXTSPEC support
- [ ] Add parlev support

## Material

- [ ] Refresh material data with ICRU90 density / stopping-power /
      mean-excitation updates
- [ ] Derive number densities and electron densities
- [ ] Build optical-depth tables
- [ ] Build nuclear target-sampling tables per material
- [ ] Support lazy extension of material/projectile runtime tables at batch
      boundaries
- [ ] Add batch/SIMD lookup helpers for stopping power and CSDA range
- [ ] Add batch inverse-range-to-energy lookup for residual-range transport

## Geometry / GEMCA

- [ ] Add AVX2 batch path for `eval_distance` /
      `osh_gemca_runtime_get_distance_batch`
- [ ] Flatten `insns_flat[]` + `insn_begin[]` in GEMCA runtime to remove
      per-zone heap-pointer layout
- [ ] Add Jacobs voxel-traversal dispatch in `eval_distance`
- [ ] Add voxel geometry support in the GEMCA / geometry layer

## Transport

- [ ] Gaussian MCS mode
- [ ] Vavilov energy straggling
- [ ] Nuclear interactions and secondary particles
- [ ] Batch ion-step phases around runtime lookup hot spots
- [ ] Measure SIMD benefit for mixed-material pools vs same-material/species
      micro-batches

## Scoring

- [ ] Cylindrical `(R,Z)` mesh scoring
- [ ] Zone scoring
- [ ] Voxel scoring
- [ ] LET scoring (`DLET`, `TLET`)
- [ ] Alanine detector response
- [ ] MCPL phase-space output

## CT / Voxel Workflow

- [ ] `HU -> material_idx` segmentation in the voxel / CT layer
- [ ] `HU -> rho` density calibration as piecewise-linear table
- [ ] `HU -> WEPL` calibration table
- [ ] Prefer mass-normalized transport tables so local density scaling stays
      cheap

## Validation / Tests

- [ ] Add a validate-mode integration test path for the current executable flow
- [ ] Improve validation diagnostics with file / line propagation from parsers
- [ ] Add more direct API tests for the new data-only setters
  - `osh_material_dedx_set(...)`
  - `osh_beam_spots_set(...)`

## Packaging / Install

- [ ] Decide whether the public library should ship as static-only,
      shared-only, or both
- [ ] Install the public library targets and headers under
      `include/openshieldhit/`
- [ ] Keep internal headers and the internal CLI parser non-installed

## Cleanup / Follow-Up

- [ ] Revisit data-module structure (`material/`, `particle/`, embedded tables)
- [ ] Rename `struct ray` to `struct ray_v` throughout

## Notes

- Core modules should not parse files or own input-format-specific I/O.
- Diagnostics are caller-owned via borrowed `osh_diag_sink` instances.
- App code may parse any format it wants, then populate the public cold structs
  through data-only APIs.
- Cold workspaces are the stable user-facing model; runtime structs are derived
  compile products used by simulation and transport.
