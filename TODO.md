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

### M3b — Array layout for cache-friendly voxel traversal
**Goal:** Make voxel array access cache-efficient for arbitrary ray directions, with infrastructure to compare layout strategies.

Status: not started (depends on M3; must land before M4 locks the index contract).

**Problem:** Row-major layout gives 524 KB stride per z-step for a 512×512 CT — guaranteed L3 miss. Same problem applies to RTDOSE dose accumulation (read-modify-write per particle crossing).

**Design: `tile_order` field on `struct osh_raytrace_grid`** controls which index formula Jacobs uses:

```c
uint8_t tile_order; /* 0=row-major (default/baseline), 8=Morton-8×8×8, 1..7=axis-perm (future) */
```

Jacobs line 146 switches on this — perfectly predicted branch, zero hot-path cost. Row-major stays available as benchmark baseline.

**Layout A — Morton-tiled (`tile_order=8`, primary candidate):**
`Tx×Ty×Tz` tiles of 512 voxels in row-major tile order, 3-way Morton order within each tile. Worst-case locality: ~1 KB per 8-voxel run in any axis.
```c
size_t Tx=(n[0]+7u)>>3u, Ty=(n[1]+7u)>>3u;   /* once before loop */
idx = ((ix>>3u) + Tx*((iy>>3u) + Ty*(iz>>3u))) * 512u
    + _m3x[ix&7u] | _m3y[iy&7u] | _m3z[iz&7u];
```
LUTs: `_m3x={0,1,8,9,64,65,72,73}`, `_m3y={0,2,16,18,128,130,144,146}`, `_m3z={0,4,32,36,256,260,288,292}`

**Layout B — Axis-permuted row-major (`tile_order=1..7`, "Bassler algorithm"):**
Permute storage axes so the beam-dominant direction is innermost (stride-1). Optimal for near-axial rays: a ray within ~30° of the inner axis stays stride-1; Morton's worst-case for the same ray is ~8 voxels within a tile. For oblique rays (near 45°) Morton wins.

**Runtime layout selection (future):** At geometry-compile time, when both the CT grid and the beam direction(s) are known, choose the layout automatically:
- Beam angle within ~30° of a principal axis → Layout B (Bassler): rotate axis assignment so dominant direction is innermost
- Beam angle near 45° or multi-angle arc → Layout A (Morton): balanced locality in all directions

Selection lives in the geometry-compile step (`osh_gemca_compile` or DCM parse), not in the hot path. Requires passing gantry/couch angles to the reorder function. The `tile_order` field on `osh_raytrace_grid` carries the result into Jacobs with no runtime branching cost.

**Shared infrastructure (CT + RTDOSE):**
- New `src/common/osh_voxel_order.h` — `tile_order` constants + inline index helpers
- New `src/common/osh_voxel_order.c` — `osh_voxel_reorder(src, nx, ny, nz, element_size, tile_order, out_len)`: generic reorder for any element size; used at CT load and RTDOSE array init
- `src/common/raytrace/osh_raytrace.h` — add `uint8_t tile_order` to `struct osh_raytrace_grid`
- `src/common/raytrace/osh_raytrace_jacobs_msh.c` line 146 — switch on `grid->tile_order`
- `src/apps/osh/osh_geometry_parse.c` — call `osh_voxel_reorder()` after pinning CT pixels; set `ct_grid.tile_order`; update `b->n_hu`
- New `tests/unit/test_osh_voxel_order.c` — LUT completeness; round-trip 16×16×16; non-power-of-8 dims
- `tests/unit/test_osh_geometry_dcm.c` — replace `hu[n_hu-1]` with interior voxel spot-check

Note: RTDOSE dose array (M6) uses same `osh_voxel_reorder` — defer to M6 but design the function generically now.

Verify: row-major baseline unchanged (`tile_order=0`); Morton round-trip correct; traversal spot-checks pass for both layouts.

---

### M3c — HU density lookup table on `osh_gemca_runtime`
**Goal:** Precompute a `float hu_rho_lut[2601]` at setup time so the Jacobs traversal hot path converts HU to density with a single array access — no formula, no branch, no function pointer.

Status: complete on this branch.

Note: this is now considered a transitional ownership choice. The longer-term
design target for M5+ is to keep GEMCA geometry-only and move HU→rho/property
resolution into `material/runtime`, with GEMCA returning geometric medium
identity (`zone_idx`, `material_idx`, optional `hu`) rather than physical
material properties directly.

**Result:** both runtime LUTs are now available in the hot path:

```c
int   bin_i = rt->hu_bin_lut[hu + 1000];  /* uint8_t — material bin  */
float rho_i = rt->hu_rho_lut[hu + 1000];  /* float   — density g/cm³ */
```

Same indexing for both tables: `hu + 1000`, with HU clamped to `[-1000, 1600]`. For non-VOX runs both pointers remain `NULL`.

Files:
- `src/gemca/voxel/osh_gemca2_voxel_hu.h` — rho LUT builders for Schneider and Permatassari
- `src/gemca/voxel/osh_gemca2_voxel_hu.c` / `osh_gemca2_voxel_hu_permatassari.c` — LUT construction
- `src/gemca/runtime/osh_gemca_runtime.h` / `.c` — owned `hu_bin_lut` and `hu_rho_lut` on the compiled runtime
- `src/simulation/osh_simulation.c` — build the active LUT pair from the selected `HUTABLE`
- `tests/unit/test_osh_voxel_hu.c` — LUT checks for both calibration tables

Verify: LUT values match the calibration functions; runtime copies/free paths are covered; both Schneider and Permatassari tables are exercised.

---

### M4 — Jacobs voxel traversal in `dist_voxel_body_rt()`
**Goal:** Implement the Jacobs traversal stub to produce a per-voxel segment list for one transport step.

Status: complete on this branch.

**Result:** `dist_voxel_body_rt()` now:

- transforms the ray to body-local coordinates via `osh_ray_transform()`
- calls `osh_raytrace_traverse()` on `body->ct_grid`
- walks the crossings in order until one of:
  - material-bin change (`bin_i != bin0`)
  - `segs_cap` reached
  - grid exit
- writes one segment per voxel crossing:
  - `ds` = path length through the voxel
  - `rho` = absolute voxel density [g/cm³] from `hu_rho_lut`
- returns the total traversed distance and the starting bin via `bin_out`

The segment payload now uses absolute density, not `rhocorr`. That matches the intended downstream use better: transport and scoring can work directly with `mass_stopping_power * rho` and `rho * ds`.

Files:
- `src/gemca/runtime/osh_gemca_runtime_voxel.c` — Jacobs traversal implementation
- `src/gemca/runtime/osh_gemca_runtime_voxel.h` — `gemca_rt_voxel_segment { ds, rho }`
- `src/gemca/runtime/CMakeLists.txt` — link against `osh_raytrace`
- `tests/unit/test_osh_voxel_body_rt.c` — direct unit coverage of the runtime voxel traversal

Verify:
- miss → `OSH_GEMCA_INFINITY`
- uniform grid → expected multi-voxel segment list
- bin change → early stop at first differing voxel
- `segs_cap` → truncation at caller capacity
- `segs == NULL` → distance-only path still works
- per-segment `rho` matches the LUT exactly

---

### M5 — Transport-facing medium query, with GEMCA/material separation
**Goal:** Finish voxel transport without leaking voxel-specific decisions into
`transport/`. GEMCA remains geometry-only; `material/runtime` owns all physical
property lookup (density, CSDA/range interpretation, radiation length, etc.).

Status: redesign in progress (replaces the earlier unified `segs[]` plan).

**Revised direction:**

- GEMCA should answer geometry questions only:
  - current `zone_idx`
  - current `material_idx`
  - optional current `hu`
  - distance to the next geometric discontinuity relevant for the current step
- transport remains a client of both GEMCA and `material/runtime`
- `material/runtime` resolves the physical properties from the GEMCA return
  value:
  - nominal density for analytic zones
  - HU→rho and any future HU-dependent property override for voxel zones
  - stopping/range/radiation-length data keyed by `material_idx`

**Key rule:** GEMCA must not depend on `material_rt`. The ownership boundary is:

- GEMCA: geometry, voxel traversal, zone/material identity, optional HU sample
- material/runtime: density/property lookup from `(material_idx, has_hu, hu)`
- transport: combines the two to perform stepping physics

**Working model for now:** transport treats each voxel like a zone. One step sees
exactly one current medium value (`rho` at the current point), not a segment
list spanning multiple voxels. This is the simplest pipeline that still leaves
room for future batch/SIMD work.

**Planned transport-facing data shape (name TBD):**

```c
struct osh_medium_ref {
    size_t zone_idx;
    size_t material_idx;
    char has_hu;
    int16_t hu;
};
```

Interpretation:
- analytic zone: `has_hu = 0`, `material_idx = zone material`
- voxel zone: `has_hu = 1`, `material_idx = voxel-bin material`,
  `hu = current voxel HU`

Then `material/runtime` provides scalar and later batched lookup helpers such as:
- density from `osh_medium_ref`
- CSDA/range lookup from `osh_medium_ref`
- radiation length / `Z_mean` / `Z/A` from `osh_medium_ref`

This keeps the future batch path natural: GEMCA can return `N` medium refs for
`N` particles, and material/runtime can evaluate `N` densities / ranges /
scattering properties without transport branching on voxel-ness.

**Immediate cleanup targets:**
- Remove the voxel/non-voxel branch from `src/transport/osh_transport_ion.c`
- Remove direct transport calls to `dist_voxel_body_rt()`
- Keep only one transport-visible geometry query at the GEMCA runtime layer for
  the current-step medium/boundary information
- Rename the voxel helper eventually: it is really a zone/medium query helper,
  not a transport-facing “body distance” API
- Replace the `n_step_segments == 0` fallback-to-vacuum path with a relookup /
  no-step outcome

**Physics policy for the basic M5 pipeline:**
- one current `rho` per step (the value at the entry point/current voxel)
- CSDA and straggling use that one `rho`
- if MCS remains enabled in voxel steps, the post-hinge tail must clip to the
  current voxel exit, not the whole same-bin zone
- more advanced multi-voxel / integrated-density algorithms are deferred until
  the ownership boundary is stable

Verify:
- voxel step uses the current voxel’s HU/rho only
- leaving a voxel volume triggers a clean relookup, not fake vacuum transport
- non-voxel transport remains unchanged
- transport has no voxel-specific knowledge; voxel handling stays entirely inside GEMCA/material

---

### M6 — RTDOSE scoring geometry (`OSH_SCORING_GEO_VOXEL`)
**Goal:** Wire up `OSH_SCORING_GEO_VOXEL` (already in enum as type 4) to score onto RTDOSE grid via Jacobs raytrace with per-voxel density weighting.

Status: not started (depends on the revised M5 ownership split).

Files:
- `include/openshieldhit/scoring.h` — add `char *rtdose_path` to `struct osh_scoring_geometry_def` (for `kind == "Voxel"`)
- `src/scoring/runtime/osh_scoring_geometry_runtime.h` — add to `struct osh_scoring_geometry_runtime`:
  ```c
  struct osh_raytrace_grid rtdose_grid;
  int16_t const *ct_hu;           /* borrowed, for density lookup */
  double transform[16];           /* same rotation as DCM body */
  ```
- `src/scoring/runtime/osh_scoring_step.c` — new branch for `OSH_SCORING_GEO_VOXEL`:
  - Basic M5-compatible version: score one current voxel step at a time
  - Later extension: accept a batched/current-step medium representation rather
    than calling voxel traversal from scoring
  - Intersect the current transport step against the RTDOSE grid (1D interval
    scan along the ray)
  - Accumulate dose per RTDOSE bin weighted by local `rho × sub_ds`
- `src/scoring/runtime/osh_scoring_compile.c` — populate `rtdose_grid` and `transform` from `osh_dicom_rtdose` struct; handle non-uniform z (frame_offsets array): treat as uniform if spacing is constant, otherwise use linear search for z-bin
- `src/apps/osh/osh_scoring_parse_geometry.c` — parse `Voxel /path/to/rtdose.dcm` keyword; store path in `geometry_def.rtdose_path`
- `src/apps/osh/osh_app_osh.c` — read RTDOSE via `osh_dicom_rtdose_read`, pass pointer into scoring compile

Coordinate dependency: RTDOSE lives in the same IEC beam frame as the CT body. The rotation matrix from the DCM body's `t[]` is shared — pass it to the scoring geometry runtime at compile time (store on cold geometry or retrieve from the geometry workspace).

Verify: known pencil-beam traversal; scored dose distribution matches
analytical WEPL integral.

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
                                   M3b (voxel array layout / Morton)
                                      │
                                   M3c (HU density LUT on osh_gemca_runtime)
                                      │
                              M4 (Jacobs + stop_reason in dist_voxel_body_rt)
                                      │
                              M5 (unified segs[] in transport kernel)
                                      │
M1 ──────────────────────── M6 (RTDOSE scoring)
                                      │
                    M2 + M6 ── M7 (CLI plumbing)
```

M1 and M2 are independent and can start in parallel on sub-branches.

---

## Open Design Questions (resolve before M2/M6)

1. **Shared transform matrix between geometry and scoring**: Store `double vox_matrix[16]` on the cold `struct body`; pass it to scoring compile step via the app orchestration layer, not via a separate cold scoring struct.

2. **Transitional LUT storage for hot path**: `uint8_t hu_bin_lut[2601]` and
   `float hu_rho_lut[2601]` currently live on `osh_gemca_runtime` (owned, freed
   by `osh_gemca_runtime_free()`). This is good enough for M3c/M4, but the
   target M5+ architecture is to move HU→rho/property resolution into
   `material/runtime` and keep GEMCA geometry-only.

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
