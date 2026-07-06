## Near-Term API / Architecture

- [ ] Revisit remaining path/file metadata in public cold workspaces
  (`wdir`, `fname`, and similar fields should disappear as app-owned path
  handling is cleaned up)
- [ ] Add result merge API for embarrassingly parallel runs ([#230])
- [x] Add chunked / partial run control (run in batches, inspect partial
  results, save explicitly) — batch-aware checkpoint scheduler ([#195]/[#207]) plus
  periodic/on-demand family-exact dumps ([#193]). The per-worker *merge* step is
  still open (tracked in [#161]).
- [ ] Finish naming consistency: remaining stale "prepare" wording in
  messages/comments

## Beam

- [ ] Add ridge-modulator / ripple-filter support ([#42])
- [ ] Add MCPL phase-space import ([#41])
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
- [x] Flatten `insns_flat[]` + `insn_begin[]` in GEMCA runtime to remove
      per-zone heap-pointer layout — done together with hoisting the
      membership evaluators into `OSH_HD` headers (shared with the CUDA
      backend; exact-parity zone kernel in `src/gpu/`)
- [ ] Rename the remaining voxel distance helper to reflect its role as a
      geometry/runtime current-voxel boundary query

## Transport

- [x] Vavilov + Landau energy straggling (STRAGG 2, [#190]) — branch
      `190-physics-add-vavilov-and-landau-straggling`.  DONE: clean-room fit to
      the exact Vavilov (1957) distribution + universal Landau; pole-free
      polynomial+Chebyshev inverse-CDF evaluators; κ-dispatch (Gaussian ≥10 /
      Vavilov / Landau <0.01) wired into ion step.  Validated vs SH12A: distal
      80–20% fluence falloff width 0.4522 vs 0.4533 cm (0.2%); DLET at matched
      fluence within ~2%.  STRAGG 0/1 byte-identical.  No GEANT3/Thomsen numbers.
- [ ] Profile STRAGG 2 Vavilov overhead on the straggling benchmark ([#211]); current
      20k-history timing is still ~2x faster than SH12A, but OSH's relative
      STRAGG 2 cost over STRAGG 0 is larger.  Check
      `osh_physics_strag_vavilov_lambda()` for cheap wins before broader
      transport refactors.
- [ ] Urban energy-loss fluctuation (STRAGG 3, reserved) — needs δ-ray cut +
      restricted stopping power; tied to future delta-electron transport
- [ ] Nuclear fragmentation — secondary particle transport (SMM, Bondorf et al.) ([#176])
- [ ] Batch ion-step phases around runtime lookup hot spots

## Neutron Transport ([#154], merged via branch 154-neutron-transport-minimal-model)

Fast-neutron transport is implemented as a condensed model.  Neutrons produced
by abrasion and Fermi break-up are banked in the neutron pool and drained by
`osh_transport_neutron_run()` after each ion pass.

### What is implemented

- **Tier-1 cross sections**: condensed JEFF-4.0 PENDF0K tables (31-point irregular
  grid, 1 meV–20 MeV) in `osh_neutron_xsec_data.h` for 35 nuclides covering
  tissue, air, bone, detectors, and shielding materials.
  Condensing automated by `tools/condense_neutron_xsec.py`.
- **Tier-2 fallback**: Tripathi `σ_R` + geometric `σ_el` for any (Z,A) not in Tier-1.
- **Neutron pool**: full SoA (pos, dir, energy, weight, prim_idx, gen, RNG stream).
  Capacity = nstat; drained in wavefront batches sized to the geometry-scratch buffer.
- **Transport loop**: GEMCA boundary distances, free-path sampling, reaction dispatch.
- **Reaction channels**: elastic (isotropic CM; n-p recoil proton returned),
  `(n,γ)` capture, `(n,p)`/`(n,α)` two-body, compound nucleus → Fermi break-up
  or heavy-A sink.
- **Natural element isotope expansion**: A=0 material entries expanded to all
  naturally occurring isotopes with abundance-weighted number densities at
  nuclear-handler compile time.
- **Family scheduler integration**: ion pass fills the neutron pool; scheduler
  drains it; pool is sized to accumulate across all ion wavefront batches.

### Open items

- [ ] Charged secondaries from neutron reactions (n,p)/(n,α) fed back into the
      ion transport family (currently deposited locally)
- [ ] Local neutron energy deposits scored (point-deposit scoring path)
- [ ] Thermal-neutron physics below 1 eV (separate issue [#178])
- [ ] Validation: compare neutron fluence to FLUKA or TOPAS reference for
      130 MeV protons on water

### Separate issues (not in this branch)

- Neutron kerma scorer (fluence × kerma factor, like sh12A) — a scorer, not
  transport; independent of the above
- U-235/238 fission kinematics

## GPU backend ([#231])

Measured plan with milestones, revised decisions, and the post-mortem of
the first attempt: `docs/dev/gpu_port_plan.md`.  State: device-compilable
core slice (RNG, vector, material, atomic physics, GEMCA membership) with
nvcc compile guards; exact-parity zone kernel + RNG/atomics benches on
A100.  Next milestone: G1 tracer bullet (c1 depth-dose fully on device).

## Scoring

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

## DICOM Study Recalculator (future app, RTPLAN support [#92])

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

<!-- Issue / PR references (reference-style links; render as clickable #nnn on GitHub). -->
[#41]: https://github.com/openshieldhit/openshieldhit/issues/41
[#42]: https://github.com/openshieldhit/openshieldhit/issues/42
[#92]: https://github.com/openshieldhit/openshieldhit/issues/92
[#154]: https://github.com/openshieldhit/openshieldhit/issues/154
[#161]: https://github.com/openshieldhit/openshieldhit/issues/161
[#176]: https://github.com/openshieldhit/openshieldhit/issues/176
[#178]: https://github.com/openshieldhit/openshieldhit/issues/178
[#190]: https://github.com/openshieldhit/openshieldhit/issues/190
[#193]: https://github.com/openshieldhit/openshieldhit/issues/193
[#195]: https://github.com/openshieldhit/openshieldhit/issues/195
[#207]: https://github.com/openshieldhit/openshieldhit/pull/207
[#211]: https://github.com/openshieldhit/openshieldhit/issues/211
[#230]: https://github.com/openshieldhit/openshieldhit/issues/230
[#231]: https://github.com/openshieldhit/openshieldhit/issues/231
