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
- [ ] Rename the remaining voxel distance helper to reflect its current role as
      a geometry/runtime current-voxel boundary query

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

- [ ] Finish legacy `VOX` card parser path if the old `.hed` / `.ctx` workflow
      should remain user-facing
- [ ] Add axis-permuted row-major voxel layouts once their index contract and
      selection heuristic are settled
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
- Track particles through the CT volume using Jacobs voxel traversal with current-voxel HU→density mapping
- Stop voxel steps at the current voxel boundary for the basic M5 transport path; multi-voxel density integration is deferred
- Score dose onto an RTDOSE grid using the same Jacobs raytrace
- Apply gantry/couch rotation via a single precalculated 4×4 transform matrix; no IEC coordinate-system rewrite needed

**Preliminary work already in place:**
- `OSH_GEOMETRY_BODY_VOX` (type 5) reserved in `include/openshieldhit/geometry_defs.h`
- `src/voxel/`: HU table constants, Schneider/Permatassari material registration, HU→bin LUT builders, HU→rho LUT builders, and WEPL utility
- DCM-backed voxel body setup in `src/apps/osh/` and `src/gemca/osh_gemca2_calc_body.c`
- `gemca_rt_zone.voxel_body_idx` + `GEMCA_RT_PUSH_VOXEL_BODY` opcode with current-voxel distance dispatch
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
This calls `osh_voxel_register_schneider_materials(wm)` or `osh_voxel_register_permatassari_materials(wm)`. A single `HUTABLE` card per material file; error on duplicate. The two-registration-function architecture avoids including both data headers in the same translation unit.

---

## Milestones

### M1 — Schneider material registration + HU→bin LUT ✅ COMPLETE
**Goal:** Materialise all 24 Schneider bins into `osh_material_workspace` and build the O(1) HU→bin lookup table.

Status: implemented. Schneider data and registration now live in `src/voxel/`; `_nelm` is 12.

### M1b — Permatassari 2020 data + HU table selection
**Goal:** Add Permatassari 2020 as a second HU calibration table and a parser keyword to select it.

Status: implemented in `src/voxel/`.

Files:
- `src/voxel/osh_voxel_mat_schneider2000.h` — Schneider material data
- `src/voxel/osh_voxel_mat_permatassari2020.h` — Permatassari material data
- `src/voxel/osh_voxel_hu_lut.c` — Schneider registration, HU→bin, HU→rho, WEPL utility
- `src/voxel/osh_voxel_hu_lut_permatassari.c` — Permatassari registration, HU→bin, HU→rho
- `src/apps/osh/osh_material_parse.c` — top-level `HUTABLE` handler
- `tests/unit/test_osh_voxel_hu.c` — registration and LUT checks for both calibration tables

---

### M2 — DCM parser card + DICOM → VOX body parameters
**Goal:** Parse `DCM` card, read CT via DICOM, and fill the exact same VOX body fields that a hand-written `VOX` card would produce. No new body type needed — `OSH_GEOMETRY_BODY_VOX` is used throughout.

Status: implemented for DCM-backed voxel bodies. Legacy `VOX` card parsing
still returns an explicit TODO parse error.

HU ownership: the app layer (e.g. `osh_app_osh.c`) holds a `int16_t *` allocation for the CT volume. `b->hu` is a borrowed pointer into it. `struct body` does not own or free it.

Files:
- `src/apps/osh/osh_geometry_parse.c` — parses `DCM name dir gantry couch tx ty tz`, reads CT metadata, reorders HU storage, and fills VOX-compatible body fields
- `src/gemca/osh_gemca2_calc_body.c` — consumes the VOX arg layout, builds the voxel enclosure surfaces, and builds the transform matrix
- `src/gemca/osh_gemca2.h` — cold body stores `ct_grid`, borrowed `hu`, and `n_hu`
- `tests/unit/test_osh_geometry_dcm.c` — parser, prepare/compile propagation, transform, and legacy VOX TODO coverage

Completed:
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

Status: implemented for row-major and Morton-8 layouts. Axis-permuted row-major
layouts remain future work.

**Problem:** Row-major layout gives 524 KB stride per z-step for a 512×512 CT — guaranteed L3 miss. Same problem applies to RTDOSE dose accumulation (read-modify-write per particle crossing).

**Design: `tile_order` field on `struct osh_raytrace_grid`** controls which index formula Jacobs uses:

```c
uint8_t tile_order; /* 0=row-major (default/baseline), 8=Morton-8×8×8, 1..7=axis-perm (future) */
```

Raytrace implementations switch on this layout field while emitting flat voxel
indices. Row-major stays available as benchmark baseline.

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
- `src/common/osh_voxel_order.h` — `tile_order` constants + inline index helpers
- `src/common/osh_voxel_order.c` — `osh_voxel_reorder(src, nx, ny, nz, element_size, tile_order, out_len)`: generic reorder for any element size
- `src/common/raytrace/osh_raytrace.h` — `uint8_t tile_order` on `struct osh_raytrace_grid`
- `src/common/raytrace/osh_raytrace_*_msh.c` — emit row-major or Morton indices according to `grid->tile_order`
- `src/apps/osh/osh_geometry_parse.c` — reorder CT pixels, set `ct_grid.tile_order`, and update `b->n_hu`
- `tests/unit/test_osh_voxel_order.c` — LUT completeness; round-trip 16×16×16; non-power-of-8 dims
- `tests/unit/test_osh_geometry_dcm.c` — interior voxel spot-check for reordered CT storage

Note: RTDOSE dose array (M6) uses same `osh_voxel_reorder` — defer to M6 but design the function generically now.

Verify: row-major baseline unchanged (`tile_order=0`); Morton round-trip correct; traversal spot-checks pass for both layouts.

---

### M3c — HU runtime lookup tables
**Goal:** Precompute HU lookup tables so voxel transport uses array lookups instead of formula calls in the hot path.

Status: complete.

Ownership after the M5 cleanup:
- GEMCA runtime owns `hu_bin_lut[2601]` because it resolves current voxel HU to
  material-bin identity.
- Material runtime owns `hu_rho_lut[2601]` because it resolves physical density
  from the current medium reference.

**Result:** both runtime LUTs are available in the hot path:

```c
int   bin_i = geom_rt->hu_bin_lut[hu + 1000];      /* uint8_t — material bin */
float rho_i = material_rt->hu_rho_lut[hu + 1000];  /* float — density g/cm³ */
```

Same indexing for both tables: `hu + 1000`, with HU clamped to `[-1000, 1600]`. For non-VOX runs both pointers remain `NULL`.

Files:
- `src/voxel/osh_voxel_hu_lut.h` / `.c` — Schneider builders and public voxel calibration API
- `src/voxel/osh_voxel_hu_lut_permatassari.c` — Permatassari builders
- `src/gemca/runtime/osh_gemca_runtime.h` / `.c` — owned `hu_bin_lut` on the compiled geometry runtime
- `src/material/runtime/osh_material_runtime.h` and `osh_material_compile.c` — owned `hu_rho_lut` on material runtime tables
- `src/simulation/osh_simulation.c` — passes the selected `HUTABLE` into both compile steps
- `tests/unit/test_osh_voxel_hu.c` — LUT checks for both calibration tables

Verify: LUT values match the calibration functions; runtime allocation/free paths are covered; both Schneider and Permatassari tables are exercised.

---

### M4 — Jacobs voxel traversal in `dist_voxel_body_rt()`
**Goal:** Implement voxel traversal and expose the current voxel boundary distance for one transport step.

Status: complete.

**Result:** `dist_voxel_body_rt()` now:

- transforms the ray to body-local coordinates via `osh_ray_transform()`
- calls `osh_raytrace_traverse()` on `body->ct_grid`
- uses only the first crossing for the current M5 transport policy
- fills at most one step segment:
  - `ds` = path length through the current voxel
  - `rho` = 0.0; density is resolved by material/runtime from `osh_zone_ref`
- returns the current voxel-exit distance and the current voxel bin via `bin_out`

Files:
- `src/gemca/runtime/osh_gemca_runtime_voxel.c` — Jacobs traversal implementation
- `src/gemca/runtime/osh_gemca_runtime_voxel.h` — current voxel distance helper
- `src/gemca/runtime/CMakeLists.txt` — link against `osh_raytrace`
- `tests/unit/test_osh_voxel_body_rt.c` — direct unit coverage of the runtime voxel traversal

Verify:
- miss → `OSH_GEMCA_INFINITY`
- uniform grid → one current-voxel segment
- bin change in the next voxel does not affect the current step distance
- `segs == NULL` → distance-only path still works
- Morton-backed HU storage is honored by zone-ref lookup

---

### M5 — Transport-facing medium query, with GEMCA/material separation
**Goal:** Finish voxel transport without leaking voxel-specific decisions into
`transport/`. GEMCA remains geometry-only; `material/runtime` owns all physical
property lookup (density, CSDA/range interpretation, radiation length, etc.).

Status: basic M5 behavior implemented; remaining work is cleanup/naming and
future integrated-density design.

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

**Completed cleanup:**
- Transport no longer calls `dist_voxel_body_rt()` directly.
- Transport asks GEMCA for current medium reference and current-medium boundary
  distance, then asks material/runtime for density.
- Voxel distance now stops at the current voxel boundary.
- `osh_gemca_runtime_get_zone_ref_batch()` resolves HU through the same
  layout-aware storage indexing used by raytrace.

**Remaining cleanup targets:**
- Rename the voxel helper eventually; it is really a geometry/runtime
  current-voxel boundary helper, not a transport-facing body API.
- Consider combining current medium reference and boundary distance into one
  GEMCA batch query once the interface stabilizes.

**Physics policy for the basic M5 pipeline:**
- one current `rho` per step (the value at the entry point/current voxel)
- CSDA and straggling use that one `rho`
- if MCS remains enabled in voxel steps, the post-hinge tail must clip to the
  current voxel exit
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

Status: not started (basic M5 ownership split is now in place).

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

Status: not started (mostly depends on M6 scoring).

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
                                   M3c (HU runtime LUTs)
                                      │
                              M4 (current voxel boundary distance)
                                      │
                              M5 (current-medium transport)
                                      │
M1 ──────────────────────── M6 (RTDOSE scoring)
                                      │
                    M2 + M6 ── M7 (CLI plumbing)
```

M1/M2/M3/M4 and the basic M5 path are now in place; M6/M7 remain the main
voxel-workflow milestones.

---

## Open Design Questions

1. **Shared transform matrix between geometry and scoring**: Store `double vox_matrix[16]` on the cold `struct body`; pass it to scoring compile step via the app orchestration layer, not via a separate cold scoring struct.

2. **Future integrated-density voxel transport**: The basic M5 path treats one
   voxel as the current medium for one step. A later design can restore
   multi-voxel density integration, but it should do so without putting
   material/runtime ownership back into GEMCA.

---

## Critical Files

- [src/gemca/osh_gemca2_calc_body.c](src/gemca/osh_gemca2_calc_body.c)
- [src/gemca/runtime/osh_gemca_runtime.c](src/gemca/runtime/osh_gemca_runtime.c)
- [src/gemca/runtime/osh_gemca_runtime.h](src/gemca/runtime/osh_gemca_runtime.h)
- [src/gemca/runtime/osh_gemca_runtime_voxel.c](src/gemca/runtime/osh_gemca_runtime_voxel.c)
- [src/voxel/osh_voxel_hu_lut.h](src/voxel/osh_voxel_hu_lut.h)
- [src/material/runtime/osh_material_runtime.h](src/material/runtime/osh_material_runtime.h)
- [src/scoring/runtime/osh_scoring_step.c](src/scoring/runtime/osh_scoring_step.c)
- [src/scoring/runtime/osh_scoring_compile.c](src/scoring/runtime/osh_scoring_compile.c)
- [src/apps/osh/osh_geometry_parse.c](src/apps/osh/osh_geometry_parse.c)
- [src/apps/osh/osh_app_osh.c](src/apps/osh/osh_app_osh.c)
- [include/openshieldhit/geometry_defs.h](include/openshieldhit/geometry_defs.h)
