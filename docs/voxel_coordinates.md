# Voxel/CT Coordinate Conventions

This note is the canonical reference for `DCM`/`VOX` geometry coordinates.

## Scope

- `src/apps/osh/osh_geometry_parse.c` (`DCM` card mapping)
- `src/gemca/osh_gemca2_calc_body.c` (`_setup_vox` transform/surfaces)
- `src/gemca/runtime/osh_gemca_runtime.c` (M3: ct_grid/hu propagation, M4: Jacobs dispatch pending)

## Coordinate Frames

1. DICOM patient coordinates:
- Provided by CT tags (`ImagePositionPatient`, orientation cosines, spacing).
- `osh_dicom_ct.origin` is the **center** of the first voxel in the first slice.

2. VOX local coordinates:
- Internal voxel-grid frame used by VOX body parameters.
- `x0,y0,z0` is the **corner** of voxel `[0,0,0]`.
- `dx,dy,dz` is spacing in cm, `nx,ny,nz` are voxel counts.

3. Universe coordinates:
- Geometry frame used by GEMCA body transforms.
- `tx,ty,tz` places the VOX local corner in universe.

## Current VOX/DCM Argument Layout

The current internal VOX payload (`OSH_GEMCA_NARGS_VOX = 14`) is:

- `0..2`: `x0,y0,z0` (local voxel-corner) [cm]
- `3..5`: `dx,dy,dz` [cm]
- `6..8`: `nx,ny,nz`
- `9..10`: `gantry,couch` [deg]
- `11..13`: `tx,ty,tz` [cm], world position of local corner

For current `DCM` parsing:

- `x0=y0=z0=0`
- `dx,dy,dz` from CT spacing (`mm -> cm`)
- `nx,ny,nz` from CT dimensions (`cols,rows,n_slices`)
- input `tx,ty,tz` is interpreted as **first-voxel center** and converted to corner by subtracting `0.5*spacing` per axis

## Why Half-Voxel Shift Is Needed

DICOM `ImagePositionPatient` and typical RT-style placement values refer to a voxel center.
Internal VOX placement is corner-based.
Therefore:

- `tx_corner = tx_center - 0.5*dx`
- `ty_corner = ty_center - 0.5*dy`
- `tz_corner = tz_center - 0.5*dz`

## RTPLAN/Isocenter Chain (Future)

When RTPLAN parsing is added:

1. CT and RTPLAN must share FrameOfReferenceUID.
2. RTPLAN isocenter is read in patient coordinates.
3. Universe is typically chosen as:
   - `U = P_patient - I_patient`
   so isocenter is at universe origin.
4. VOX placement and gantry/couch transforms are then applied in that universe.

## Open Items

- Non-axial CT orientation is currently rejected in DCM parsing with a clear
  error (axial-only support for now). Future work will apply
  `row_cosine`/`col_cosine` to placement for full oblique support.
- `ct_grid` is propagated to cold/runtime body structs (`struct body`, `gemca_rt_body`).
- `hu` ownership: workspace-owned HU buffer is borrowed (no copy) by prepared/runtime GEMCA body structs.
- Jacobs voxel traversal dispatch (`GEMCA_RT_PUSH_VOXEL_BODY`) is pending (M4).
