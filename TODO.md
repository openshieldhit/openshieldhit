## Near-Term API / Architecture

- [ ] Revisit remaining path/file metadata in public cold workspaces
  (`wdir`, `fname`, and similar fields should disappear as app-owned path
  handling is cleaned up)
- [ ] Add result merge API for embarrassingly parallel runs
- [ ] Add chunked / partial run control (run in batches, inspect partial
  results, merge, save explicitly)
- [ ] Finish naming consistency: remaining stale "prepare" wording in
  messages/comments

## Beam

- [ ] Add ridge-modulator / ripple-filter support
- [ ] Add MCPL phase-space import
- [ ] Add EXTSPEC support
- [ ] Add parlev support

## Material

- [ ] Refresh material data with ICRU90 density / stopping-power /
      mean-excitation updates
- [ ] Derive number densities and electron densities
- [ ] Build optical-depth / nuclear target-sampling tables per material
- [ ] Support lazy extension of material/projectile runtime tables at batch
      boundaries
- [ ] Add batch/SIMD lookup helpers for stopping power and CSDA range

## Geometry / GEMCA

- [ ] Add AVX2 batch path for `eval_distance` /
      `osh_gemca_runtime_get_distance_batch`
- [ ] Flatten `insns_flat[]` + `insn_begin[]` in GEMCA runtime to remove
      per-zone heap-pointer layout
- [ ] Rename the remaining voxel distance helper to reflect its role as a
      geometry/runtime current-voxel boundary query

## Transport

- [ ] Gaussian MCS mode
- [ ] Vavilov energy straggling
- [ ] Nuclear fragmentation — secondary particle transport (SMM, Bondorf et al.)
- [ ] Batch ion-step phases around runtime lookup hot spots

## Scoring

- [ ] Cylindrical `(R,Z)` mesh scoring
- [ ] Zone scoring
- [ ] Alanine detector response
- [ ] MCPL phase-space output

## CT / Voxel / RTDOSE Workflow

Open:
- [ ] DOSE scoring with per-voxel density weighting for CT geometries (current
      implementation uses `st->rho` = transport zone density; CT voxel-level
      density correction requires per-voxel rho lookup in the scorer)
- [ ] Add axis-permuted row-major voxel layouts once index contract is settled
- [ ] Prefer mass-normalised transport tables so local density scaling is cheap

## Validation / Tests

- [ ] Add a validate-mode integration test path for the current executable flow
- [ ] Improve validation diagnostics with file/line propagation from parsers
- [ ] Add more direct API tests for data-only setters
  (`osh_material_dedx_set`, `osh_beam_spots_set`)

## Packaging / Install

- [ ] Decide whether the public library ships as static-only, shared-only, or both
- [ ] Install public library targets and headers under `include/openshieldhit/`
- [ ] Keep internal headers and the internal CLI parser non-installed

## Cleanup / Follow-Up

- [ ] Revisit data-module structure (`material/`, `particle/`, embedded tables)
- [ ] Rename `struct ray` to `struct ray_v` throughout

## DICOM Study Recalculator (future app)

A dedicated `apps/osh_dicom_study` (name TBD) will handle full-plan DICOM
recalculation without the file-parse overhead of `osh_sim`:

- Read CT + RTPLAN once; iterate over all beams/fields in RTPLAN.
- For each field: extract isocenter position and gantry/couch angle from
  RTPLAN, populate `beam_workspace` + geometry cold structs directly — no
  geo.dat / beam.dat round-trip.
- Spin up one `osh_sim` instance per field; merge dose results.

Positioning in `osh_sim` (manual / current design):
- The `DCM` card in geo.dat carries `gantry_deg couch_deg tx_cm ty_cm tz_cm`.
- `tx/ty/tz` should be set to `−isocenter_PCS [cm]` so that the DICOM
  isocenter maps to universe origin.
- The caller computes these from RTPLAN manually (or the future recalculator
  derives them automatically).
- This is intentional: `osh_sim` is a generic transport engine; RTPLAN
  semantics belong in the higher-level recalculator app.

## Notes

- Core modules must not parse files or own input-format-specific I/O.
- Diagnostics are caller-owned via borrowed `osh_diag_sink` instances.
- App code may parse any format it wants, then populate the public cold structs
  through data-only APIs.
- Cold workspaces are the stable user-facing model; runtime structs are derived
  compile products used by simulation and transport.
