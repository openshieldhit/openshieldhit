## Near-Term API / Architecture

- [ ] Revisit remaining path/file metadata in public cold workspaces
  (`wdir`, `fname`, and similar fields should disappear as app-owned path
  handling is cleaned up)
- [ ] Add result merge API for embarrassingly parallel runs
- [x] Add chunked / partial run control (run in batches, inspect partial
  results, save explicitly) — batch-aware checkpoint scheduler (#195/#207) plus
  periodic/on-demand family-exact dumps (#193). The per-worker *merge* step is
  still open (tracked in #161).
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

- [x] Vavilov + Landau energy straggling (STRAGG 2) — branch
      `190-physics-add-vavilov-and-landau-straggling`.  DONE: clean-room fit to
      the exact Vavilov (1957) distribution + universal Landau; pole-free
      polynomial+Chebyshev inverse-CDF evaluators; κ-dispatch (Gaussian ≥10 /
      Vavilov / Landau <0.01) wired into ion step.  Validated vs SH12A: distal
      80–20% fluence falloff width 0.4522 vs 0.4533 cm (0.2%); DLET at matched
      fluence within ~2%.  STRAGG 0/1 byte-identical.  No GEANT3/Thomsen numbers.
      Remaining nicety: render the `plot_straggling.py` strag2 overlay.
- [ ] Profile STRAGG 2 Vavilov overhead on the straggling benchmark; current
      20k-history timing is still ~2x faster than SH12A, but OSH's relative
      STRAGG 2 cost over STRAGG 0 is larger.  Check
      `osh_physics_strag_vavilov_lambda()` for cheap wins before broader
      transport refactors.
- [ ] Urban energy-loss fluctuation (STRAGG 3, reserved) — needs δ-ray cut +
      restricted stopping power; tied to future delta-electron transport
- [ ] Nuclear fragmentation — secondary particle transport (SMM, Bondorf et al.)
- [ ] Batch ion-step phases around runtime lookup hot spots

## Neutron Transport (issue #154, merged via branch 154-neutron-transport-minimal-model)

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
- [ ] Thermal-neutron physics below 1 eV (separate issue #178)
- [ ] Validation: compare neutron fluence to FLUKA or TOPAS reference for
      130 MeV protons on water

### Separate issues (not in this branch)

- Neutron kerma scorer (fluence × kerma factor, like sh12A) — a scorer, not
  transport; independent of the above
- U-235/238 fission kinematics

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
