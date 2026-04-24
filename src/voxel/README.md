# voxel — Shared CT calibration module

This module owns the HU (Hounsfield Unit) calibration tables that are needed by
both the geometry engine (`src/gemca/`) and the material runtime
(`src/material/runtime/`).

Placing this code here avoids a circular dependency: GEMCA and material must
not reach into each other, but both require the same calibration data.

## What lives here

| File | Purpose |
|------|---------|
| `osh_voxel_hu_lut.h` | Internal API: all `osh_voxel_*` functions |
| `osh_voxel_hu_lut.c` | Schneider 2000 (24 bins): register materials, build LUTs |
| `osh_voxel_hu_lut_permatassari.c` | Permatassari 2020 (40 bins): register materials, build LUTs |
| `osh_voxel_mat_schneider2000.h` | Static Schneider 2000 data tables |
| `osh_voxel_mat_permatassari2020.h` | Static Permatassari 2020 data tables |

Public constants (`OSH_HU_TABLE_*`, `OSH_VOXEL_HU_LUT_SIZE`) live in
`include/openshieldhit/voxel.h`.

## Responsibilities

- **Parse time** (`src/apps/osh/`): the application parser calls
  `osh_voxel_register_*_materials()` to add CT tissue bins to the material
  workspace, and sets `hu_table_type` on the material workspace.
- **Geometry compile** (`src/gemca/runtime/`): calls
  `osh_voxel_build_hu_bin_lut_*()` to build the HU→bin LUT inside the
  geometry runtime when `simulation_create()` passes a non-`NONE` HU table
  selector.
- **Material compile** (`src/material/runtime/`): calls
  `osh_voxel_build_hu_rho_lut_*()` to build the HU→density LUT inside the
  material runtime.

## Dependency rules

`osh_voxel` may include `openshieldhit/material.h` (public header) to
manipulate `struct osh_material_workspace`, but must not link against the
internal material or gemca libraries.
