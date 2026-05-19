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
- [ ] Nuclear interactions and secondary particles
- [ ] Batch ion-step phases around runtime lookup hot spots

## Scoring

- [ ] Cylindrical `(R,Z)` mesh scoring
- [ ] Zone scoring
- [ ] LET scoring (`DLET`, `TLET`)
- [ ] Alanine detector response
- [ ] MCPL phase-space output

## CT / Voxel / RTDOSE Workflow

Completed (branches 66-1 … 66-6):
- [x] Schneider 2000 + Permatassari 2020 HU calibration tables
- [x] `DCM` parser card: reads DICOM CT into a `VOX` body at parse time
- [x] Morton-8 and row-major voxel array layouts; `tile_order` field on grid
- [x] HU→bin and HU→rho runtime lookup tables (O(1), L1-resident)
- [x] Jacobs voxel traversal in `dist_voxel_body_rt()`
- [x] Current-medium transport: GEMCA returns zone/material/HU ref; transport
      queries material runtime for density and stopping power
- [x] `DicomRTDOSE` scoring geometry: app converts it to a plain Mesh using
      the RTDOSE grid dimensions and the patient→world coordinate offset
      stored in the CT body; library is agnostic
- [x] `DicomCT` scoring geometry: app converts the CT voxel body extent to a
      plain axis-aligned Mesh for unrotated cases; library is agnostic
- [x] RTDOSE write-back: `FileFormat RTDOSE` in detect.dat reads the DICOM
      template, overwrites pixel data, and saves a modified `.dcm`

Open:
- [ ] `DicomCT` scoring with gantry/couch rotation: the Mesh geometry needs an
      optional rotation matrix; `run_setup_voxel_scoring()` must copy the CT
      body's gantry/couch rotation (from `b->a[9..10]`) onto the scoring Mesh
      so the runtime maps universe → patient frame before binning.
      Plain `Geometry Mesh` cards are never rotated (universe frame only).
- [ ] `DicomRTDOSE` scoring with gantry/couch rotation: same rotation mechanism
      as DicomCT; RTDOSE grid lives in patient frame, so the CT body rotation
      applies here too.  Library stays DICOM-agnostic (rotated Mesh is general).
- [ ] DOSE scoring with density weighting (current RTDOSE output is ENERGY
      per primary; proper absorbed dose needs `rho × ds` weighting per voxel)
- [ ] Finish legacy `VOX` card parser path (`.hed`/`.ctx` workflow)
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
