# src/apps/osh

OpenShieldHIT-style input file parsers and the `openshieldhit` executable entry
point.  This directory is the policy boundary between the file system and the
pure simulation library.

## File overview

| File | Role |
|------|------|
| `osh_app_osh.c/h` | `main()` entry point; CLI argument parsing |
| `osh_run.c/h` | Top-level orchestration: resolve paths, load workspaces, run, save |
| `osh_beam_parse.c/h` | Parse `beam.dat` into `osh_beam_workspace` |
| `osh_geometry_parse.c/h` | Parse `geo.dat` into `osh_geometry_workspace`; handles `DCM` and `VOX` cards |
| `osh_material_parse.c/h` | Parse `mat.dat` into `osh_material_workspace`; handles `HUTABLE` keyword |
| `osh_scoring_parse*.c/h` | Parse `detect.dat` into `osh_scoring_workspace` |
| `osh_scoring_parse_geometry.c` | Per-geometry-kind keyword handlers (`Mesh`, `DicomCT`, `DicomRTDOSE`, …) |

## Orchestration (`osh_run.c`)

`osh_run()` calls four `setup_from_path` functions (one per input file), then
hands the cold workspaces to the library:

```c
osh_simulation_create(beam, geo, mat, scoring, diag, &sim);
osh_simulation_run(sim);
osh_simulation_save(sim);
osh_simulation_free(sim);
```

Between loading and `simulation_create`, `osh_run.c` performs two voxel-specific
setup steps in `run_setup_voxel_scoring()`:

### Phase 1 — DicomRTDOSE → Mesh conversion

For each scoring geometry with `kind == "dicomrtdose"`:

1. Read the RTDOSE DICOM file (from `vox_rtdose_path`) to obtain grid dimensions,
   pixel spacing, and frame offsets.
2. Look up the patient→world coordinate offset from the CT VOX body in the
   geometry workspace (`b->a[14..16]` = `tx_cm − ct_origin_cm` for each axis,
   stored by the `DCM` parser card).
3. Compute three Mesh axes (X, Y, Z) in simulation world coordinates:
   - RTDOSE DICOM origin is the first-voxel *centre* → convert to corner
   - mm → cm
   - Add patient→world offset
4. Mutate the cold geometry def: append the three axes, set `kind = "mesh"`.
   `vox_rtdose_path` is preserved on the def for the save step.

After this conversion the scoring library sees a plain Mesh geometry — no DICOM
awareness needed in the library hot path.

**Limitation:** the Mesh axes are axis-aligned in the patient frame.  Non-zero
gantry/couch angles require a coordinate transform in the scoring step before the
Mesh bin lookup, which is deferred.

### Phase 2 — DicomCT → Mesh conversion

For scoring geometries with `kind == "dicomct"`, read the world-corner
(`b->a[11..13]`), spacing (`b->a[3..5]`), and bin counts (`b->a[6..8]`) from
the CT VOX body and append X/Y/Z axes, then set `kind = "mesh"`.  The library
therefore only ever sees `"mesh"` geometries.

**Current limitation:** this conversion is axis-aligned.  It matches the
transport CT body only when gantry/couch rotation is zero.  Rotated `DicomCT`
scoring is deferred to a later follow-up that will apply the corresponding
coordinate transform during scoring.

## DCM body parameters (`b->a[]`)

The `DCM` parser card fills a VOX body with 17 numeric parameters:

| Index | Content |
|-------|---------|
| 0..2  | local voxel-grid corner (0,0,0 for DCM) |
| 3..5  | dx, dy, dz voxel spacing [cm] |
| 6..8  | nx, ny, nz voxel counts |
| 9..10 | gantry angle, couch angle [deg] |
| 11..13 | world corner of first voxel [cm] (tx/ty/tz − 0.5·spacing) |
| 14..16 | patient→world offset [cm] (tx/ty/tz − DICOM origin in cm) |

The gemca geometry engine requires at least indices 0..13 (`OSH_GEMCA_NARGS_VOX = 14`).
Indices 14..16 are extra metadata read only by `osh_run.c`.

## Output path rewriting

After loading the scoring workspace, `run_resolve_output_paths()` rewrites every
`output.filename` from the bare name in `detect.dat` to a full path under
`--outdir`.  The RTDOSE writer appends `.dcm` if the filename does not already
end in it.
