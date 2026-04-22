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
- [x] legacy logger API removed; public diagnostics live in `diag.h`, internal
      diagnostics helpers in `osh_diag`, and fatal OOM aborts in `osh_abort`

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



# VOX/DCM Voxel CT Transport — Implementation Roadmap

## Context

Adding full voxel CT transport to OpenShieldHIT for proton/ion therapy. The simulation must:
- Load a DICOM CT series and build a geometry body from it
- Track particles through the CT volume using Jacobs voxel traversal with per-voxel HU→density mapping
- Split steps at material-bin boundaries (not density transitions) to satisfy the physics invariant that material index is immutable across a step
- Score dose onto an RTDOSE grid using the same Jacobs raytrace
- Apply gantry/couch rotation via a single precalculated 4×4 transform matrix; no IEC coordinate-system rewrite needed

**Preliminary work already in place:**
- `OSH_GEOMETRY_BODY_VOX` (type 5) reserved in `include/openshieldhit/geometry_defs.h`
- `src/gemca/voxel/`: `struct voxelct`, HU→density (`osh_gemca_voxel_hu2rho`), HU→bin (`osh_gemca_voxel_hu2idx`), Schneider-2000 24-bin table
- `_setup_vox()` skeleton in `src/gemca/osh_gemca2_calc_body.c:484` with TODOs at lines 501–532
- `gemca_rt_zone.voxel_body_idx` + `GEMCA_RT_PUSH_VOXEL_BODY` opcode, stub at `src/gemca/runtime/osh_gemca_runtime.c:1720`
- Jacobs raytrace fully implemented: `src/common/raytrace/osh_raytrace_jacobs_msh.c`
- DICOM CT + RTDOSE readers complete in `src/dicom/`

**DCM is a parser convenience card that fills VOX parameters from DICOM.** Internally there is only `OSH_GEOMETRY_BODY_VOX`. Cold and runtime structs have no knowledge of DICOM.
- `DCM <name> /path/to/ct/dir <gantry_deg> <couch_deg>` — parser reads DICOM CT, extracts bounding box + grid + HU array, populates the same fields as a `VOX` card.
- `VOX <name> <filename.ctx> <gantry_deg> <couch_deg>` — points to a legacy `.hed`/`.ctx` file; origin, grid dimensions, and spacing come from the file header (same philosophy as FLUKA's `VOXELS` card where `.vxl` carries all grid metadata). No inline corner position.
- Isocenter comes from CT DICOM origin; not user-provided for DCM.

**HU array:** one flat `const int16_t[]` per CT study, owned by the app layer for the run lifetime. `gemca_rt_body.hu` is a borrowed `const int16_t *` — no copies. `osh_dicom_ct` may be freed after setup once `hu[]` is pinned.

**HU→bin lookup:** precompute a `uint8_t lut[2601]` at startup (indexed by `hu + 1000`, clamped to [-1000, 1600]). Built once from breakpoints using the existing binary search as the build step. Replaces per-call binary search with a single array access (~2.6 KB, L1 resident). Both Schneider and Permatassari tables use the same LUT convention; only the breakpoints and bin count differ.

**HU→rho:** formula differs by table:
- Schneider 2000: global piecewise linear (Eqs. 20–24 from the paper), no bin index needed
- Permatassari 2020: `factor[bin] * (1000 + HU)` — requires the bin index; caller provides it from the LUT

**Two HU calibration tables supported:**

| | Schneider 2000 | Permatassari 2020 |
|---|---|---|
| Reference | Schneider et al., PMB 45 (2000) | Permatassari et al., PMB 65 (2020), Method C |
| Bins | 24 | 40 |
| Elements | 12 (H,C,N,O,Na,Mg,P,S,Cl,Ar,K,Ca) | 25 (adds He,Li,Be,B,F,Ne,Al,Si,Ti,Fe,Zn,I,Ba) |
| Density | Eqs. 20–24 piecewise linear | `factor[bin]*(1000+HU)` |
| Source | `osh_voxel_mat_schneider2000.h` | `osh_voxel_mat_permatassari2020.h` |
| Material prefix | `schneider_00`…`schneider_23` | `permatassari2020_00`…`permatassari2020_39` |

**Known issue — Permatassari air density:** the formula gives 0 at HU=-1000. This is intentional in TOPAS (scanner air = vacuum). Our `hu2rho_permatassari2020` must clamp to a physical minimum (~0.00121 g/cm³) or return the nominal bin density for the first bin. Verify against Permatassari 2020 paper once available.

**HUTABLE parser keyword:** a new top-level keyword in the material input file (NOT inside a `MATERIAL` block, since it registers N materials at once):
```
HUTABLE Schneider2000
```
or
```
HUTABLE Permatassari2020
```
This calls `osh_gemca_voxel_register_schneider_materials(wm)` or `osh_gemca_voxel_register_permatassari_materials(wm)`. A single `HUTABLE` card per material file; error on duplicate. The two-registration-function architecture avoids including both data headers in the same translation unit (variable name conflict: both define `_nmat`, `_nelm`, etc.). Each function lives in its own `.c` file.

---

## Milestones

### M1 — Schneider material registration + HU→bin LUT ✅ COMPLETE
**Goal:** Materialise all 24 Schneider bins into `osh_material_workspace` and build the O(1) HU→bin lookup table.

Status: implemented and tested on branch `66-1_schneider_lut`. All 45 tests pass.
- `_nelm` still set to 14 (should be 12 — fix in M1b below)

### M1b — Permatassari 2020 data + HU table selection
**Goal:** Add Permatassari 2020 as a second HU calibration table and a parser keyword to select it.

Files:
- `src/gemca/voxel/osh_voxel_mat_schneider2000.h`:
  - Fix `_nelm = 14` → `_nelm = 12` (only 12 non-zero Z values in `_ct_elmz`)
  - Add reference comment: Schneider W et al., PMB 45:459–478 (2000), Table 6 + Eqs. 20–24
- `src/gemca/voxel/osh_voxel_mat_permatassari2020.h` (new):
  - `_nmat = 40`, `_nelm = 25`, `unsigned int const _ct_elmz[25]` (Z=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,26,30,53,56)
  - `int16_t const _ct_hu[41]` — HU section breakpoints from TOPAS (–1024…3072)
  - `double const _ct_density_factor[40]` — density formula coefficients from TOPAS
  - `double const _ct_mean_excitation_energy[40]` — mean excitation energies [eV] from TOPAS
  - `float const _ct_relm[40][25]` — mass fractions (already 1.0-normalized in TOPAS source)
  - Reference comment: Permatassari et al., PMB 65:ab9702 (2020), Method C; data from TOPAS `SPRtoMaterial__Brain.txt`
- `src/gemca/voxel/osh_gemca2_voxel_hu_permatassari.c` (new, separate TU to avoid name conflicts):
  - `osh_gemca_voxel_register_permatassari_materials(wm)` — same realloc pattern as Schneider; uses `mat->mean_excitation_energy` from table; normalizes fractions by actual row sum as safety guard; bin 0 → GAS, rest → CONDENSED
  - `osh_gemca_voxel_build_hu_lut_permatassari2020(lut)` — same `build_hu_lut` helper as Schneider but with `_ct_hu[41]`, 41 breakpoints
  - `osh_gemca_voxel_hu2rho_permatassari2020(int16_t hu, int bin)` — returns `_ct_density_factor[bin] * (1000.0 + hu)`; clamps result to `max(rho, 0.0001)` for the near-zero air case
- `src/gemca/voxel/osh_gemca2_voxel_hu.h`: add declarations for all three new functions
- `src/gemca/CMakeLists.txt`: add `voxel/osh_gemca2_voxel_hu_permatassari.c`
- `src/apps/osh/osh_material_parse_keys.h`: add `OSH_MATERIAL_KEY_HUTABLE "hutable"`
- `src/apps/osh/osh_material_parse.c`: add `parse_hutable()` handler; must be called at top level (not inside a `MATERIAL` block); parse one token ("schneider2000" or "permatassari2020", case-insensitive); call the matching register function; error on unknown table name or duplicate `HUTABLE`
- `tests/unit/test_osh_voxel_hu.c`: add tests for Permatassari — registration count (42 = 2+40), LUT boundary values (`lut[-1000+1000]==0`, `lut[0+1000]==14`, `lut[1600+1000]==39`), `hu2rho_permatassari2020(0, 14) ≈ 0.998`

---

### M2 — DCM parser card + DICOM → VOX body parameters
**Goal:** Parse `DCM` card, read CT via DICOM, and fill the exact same VOX body fields that a hand-written `VOX` card would produce. No new body type needed — `OSH_GEOMETRY_BODY_VOX` is used throughout.

Files:
- `src/apps/osh/osh_geometry_parse.c` — in `_parse_bodies`, handle `"DCM"` keyword:
  - Read one string token (CT dir path) + 5 doubles (`gantry_deg couch_deg tx ty tz`)
  - Call `osh_dicom_ct_read(path, &ct, diag)` immediately (parse-time; DICOM is not deferred)
  - Derive and store VOX grid args `{x0,y0,z0, dx,dy,dz, nx,ny,nz}` and transform args `{gantry,couch,tx,ty,tz}`
  - For DCM, interpret tx/ty/tz as first-voxel-center and convert to corner-based placement by subtracting `0.5*spacing`
  - Keep legacy `VOX` card parsing disabled for now with explicit TODO parse error
- `src/gemca/osh_gemca2_calc_body.c` — complete the existing TODOs in `_setup_vox` (lines 501–532):
  - Accept pre-filled voxel-grid args (already set by parser for DCM and future VOX cards)
  - Build transform matrix from gantry + couch + isocenter using `osh_vect_rot_y` / `osh_vect_rot_z` (existing pattern); verify sign conventions match IEC 61217 once, document
  - Build enclosing surfaces from `x0 + d*n` extents
- `src/gemca/osh_gemca2.h` — add to `struct body`: `struct osh_raytrace_grid ct_grid; const int16_t *hu;`

HU ownership: the app layer (e.g. `osh_app_osh.c`) holds a `int16_t *` allocation for the CT volume. `b->hu` is a borrowed pointer into it. `struct body` does not own or free it.

Progress on branch `66-2`:
- [x] `DCM` body card wired in geometry parser
- [x] DICOM CT read at parse-time from `DCM` card
- [x] VOX arg model switched to `x0,y0,z0, dx,dy,dz, nx,ny,nz, gantry,couch, tx,ty,tz`
- [x] Explicit docs on corner-vs-center convention at parser/body setup level
- [x] `_setup_vox()` consumes the new arg layout and builds correct enclosure extents
- [x] Unit tests for DCM parsing + legacy VOX TODO behavior
- [x] `ct_grid` stored on cold `struct body` and copied to runtime `gemca_rt_body`
- [x] VOX transform matrix orthonormality covered by unit tests across representative gantry/couch angles
- [x] `hu` pointer ownership/borrowing wired for DCM (workspace-owned HU -> cold body -> runtime body)
- [x] Non-axial CT currently rejected with explicit parse error (future: apply `row_cosine`/`col_cosine` in placement)

Verify (current):
- parse DICOM series via `DCM`; assert spacing/counts and transformed placement args
- assert legacy `VOX` card returns TODO parse error

Verify (remaining):
- none (M2 acceptance checks covered by unit tests)

---

### Coordinate Conventions (DCM/VOX)

- Canonical reference moved to:
  - [docs/voxel_coordinates.md](docs/voxel_coordinates.md)
- Keep code comments concise and update the design note when conventions change.

---

### M3 — CT grid in `gemca_rt_body` (runtime compile)
**Goal:** Propagate `ct_grid` and `hu` pointer from `struct body` through `osh_gemca_compile()` into `gemca_rt_body`.

Status: completed for DCM-backed voxel bodies.

Files:
- `src/gemca/runtime/osh_gemca_runtime.h` — add to `struct gemca_rt_body`:
  ```c
  struct osh_raytrace_grid ct_grid; /* only for OSH_GEMCA_BODY_VOX/DCM */
  int16_t const *hu;               /* borrowed; lifetime: app run */
  ```
- `src/gemca/runtime/osh_gemca_runtime.c` — in `setup_bodies()`, copy `ct_grid` and `hu` for DCM body type

Verify: compile a geometry with one DCM body; assert `rt->bodies[i].ct_grid.n[2] == n_slices` and `hu != NULL`.

---

### M4 — Jacobs voxel traversal in transport hot path
**Goal:** Activate `GEMCA_RT_PUSH_VOXEL_BODY` to split steps at material-bin boundaries using Jacobs raytrace.

Status: not started (depends on M3 runtime `ct_grid` + `hu` propagation).

Files:
- `src/gemca/runtime/osh_gemca_runtime.c` — fill stub at line 1720:
  1. Call `osh_raytrace_traverse()` on `body->ct_grid` to get `crossings[]`
  2. Walk crossings; call `osh_gemca_voxel_hu2idx(hu[crossing.idx])` per voxel
  3. Return accumulated `path_len` up to (not including) the first bin change as the step boundary distance; if no bin change, return the full RPP boundary distance
- New: `src/gemca/runtime/osh_gemca_runtime_voxel.c` — factor out `dist_voxel_body_rt()` to keep the main switch clean

Physics: density is continuous (`hu2rho`); bin index is the material identity. The stopping-power correction is `rho_voxel / rho_bin_nominal` × nominal dE/dx for the bin's material (standard density-scaling approach). This means the transport kernel must receive `rho_voxel` per sub-step.

Verify: pencil beam through uniform-HU phantom; step lengths sum to expected traversal; bin-boundary crossings produce truncated steps.

---

### M5 — Per-voxel density to transport kernel
**Goal:** Pass per-voxel density from the Jacobs traversal to the transport kernel for correct WEPL/range calculation.

Status: not started (depends on M4 voxel traversal dispatch).

Files:
- Introduce `struct gemca_rt_voxel_step { double dist; double rho; }` returned from `dist_voxel_body_rt()` (avoids mutating `struct step` which is transport-owned)
- Wherever `osh_gemca_runtime_get_distance` result is consumed by the transport loop: use the returned `rho` to scale stopping power for WEPL accumulation instead of the zone material's nominal density

Verify: two-tissue phantom (bone + soft tissue); scored WEPL matches analytical expectation.

---

### M6 — RTDOSE scoring geometry (`OSH_SCORING_GEO_VOXEL`)
**Goal:** Wire up `OSH_SCORING_GEO_VOXEL` (already in enum as type 4) to score onto RTDOSE grid via Jacobs raytrace with per-voxel density weighting.

Status: not started (depends on M4/M5).

Files:
- `include/openshieldhit/scoring.h` — add `char *rtdose_path` to `struct osh_scoring_geometry_def` (for `kind == "Voxel"`)
- `src/scoring/runtime/osh_scoring_geometry_runtime.h` — add to `struct osh_scoring_geometry_runtime`:
  ```c
  struct osh_raytrace_grid rtdose_grid;
  int16_t const *ct_hu;           /* borrowed, for density lookup */
  double transform[16];           /* same rotation as DCM body */
  ```
- `src/scoring/runtime/osh_scoring_step.c` — new branch for `OSH_SCORING_GEO_VOXEL`:
  - Transform step coordinates using `geo->transform` (inverse of DCM body transform → RTDOSE frame)
  - Call `osh_raytrace_traverse()` on `rtdose_grid`
  - Accumulate dose per crossing bin, weighted by `hu2rho(ct_hu[idx])` × path_len
- `src/scoring/runtime/osh_scoring_compile.c` — populate `rtdose_grid` and `transform` from `osh_dicom_rtdose` struct; handle non-uniform z (frame_offsets array): treat as uniform if spacing is constant, otherwise use linear search for z-bin
- `src/apps/osh/osh_scoring_parse_geometry.c` — parse `Voxel /path/to/rtdose.dcm` keyword; store path in `geometry_def.rtdose_path`
- `src/apps/osh/osh_app_osh.c` — read RTDOSE via `osh_dicom_rtdose_read`, pass pointer into scoring compile

Coordinate dependency: RTDOSE lives in the same IEC beam frame as the CT body. The rotation matrix from the DCM body's `t[]` is shared — pass it to the scoring geometry runtime at compile time (store on cold geometry or retrieve from the geometry workspace).

Verify: known pencil-beam traversal; scored dose distribution matches analytical WEPL integral.

---

### M7 — Integration & CLI plumbing
**Goal:** End-to-end pipeline: user specifies `DCM` in geo.dat and `Geometry Voxel` in detect.dat.

Status: not started (depends on M2–M6 core pieces).

Files:
- `src/apps/osh/osh_geometry_parse.c` — ensure `DCM` keyword triggers M2 path:
  `DCM <name> <ct_dir> <gantry> <couch> <tx> <ty> <tz>`
- `src/apps/osh/osh_scoring_parse_geometry.c` — `Voxel` geometry kind accepts RTDOSE path, no axis definitions
- `src/apps/osh/osh_run.c` — orchestrate: detect DCM body → load CT → pass to geometry compile; load RTDOSE → pass to scoring compile; wire shared transform matrix

Verify: full system test with real DICOM CT + RTDOSE pair; dose grid dimensions match RTDOSE; no infinite loops; results reproducible.

---

## Dependency Graph

```
M1 (Schneider materials)     M2 (DCM parse + DICOM→body)
                                      │
                                      M3 (ct_grid in gemca_rt_body)
                                      │
                              M4 (Jacobs hot path)
                                      │
                              M5 (per-voxel density)
                                      │
M1 ──────────────────────── M6 (RTDOSE scoring)
                                      │
                    M2 + M6 ── M7 (CLI plumbing)
```

M1 and M2 are independent and can start in parallel on sub-branches.

---

## Open Design Questions (resolve before M2/M6)

1. **Shared transform matrix between geometry and scoring**: Store `double vox_matrix[16]` on the cold `struct body`; pass it to scoring compile step via the app orchestration layer, not via a separate cold scoring struct.

2. **LUT storage for hot path**: `uint8_t hu_lut[2601]` should live on `osh_gemca_prepared` (one per run, zero overhead for non-VOX geometries). Avoids a global/module-level singleton and keeps it with the geometry runtime.

---

## Critical Files

- [src/gemca/osh_gemca2_calc_body.c](src/gemca/osh_gemca2_calc_body.c)
- [src/gemca/runtime/osh_gemca_runtime.c](src/gemca/runtime/osh_gemca_runtime.c)
- [src/gemca/runtime/osh_gemca_runtime.h](src/gemca/runtime/osh_gemca_runtime.h)
- [src/gemca/voxel/osh_gemca2_voxel_hu.c](src/gemca/voxel/osh_gemca2_voxel_hu.c)
- [src/scoring/runtime/osh_scoring_step.c](src/scoring/runtime/osh_scoring_step.c)
- [src/scoring/runtime/osh_scoring_compile.c](src/scoring/runtime/osh_scoring_compile.c)
- [src/apps/osh/osh_geometry_parse.c](src/apps/osh/osh_geometry_parse.c)
- [src/apps/osh/osh_app_osh.c](src/apps/osh/osh_app_osh.c)
- [include/openshieldhit/geometry_defs.h](include/openshieldhit/geometry_defs.h)
