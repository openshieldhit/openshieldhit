# gemca/voxel

This directory holds voxel/CT-specific helpers that sit beside GEMCA geometry
preparation but are not part of the analytic body runtime itself.

## Purpose

- HU-to-density and HU-to-material-bin calibration helpers
- material-table registration for CT transport workflows
- legacy voxel file parsing support retained from older code paths

## Files

- `osh_gemca2_voxel_hu.c` — Schneider 2000 calibration and registration
- `osh_gemca2_voxel_hu_permatassari.c` — Permatassari 2020 calibration and registration
- `osh_voxel_mat_schneider2000.h` — embedded Schneider material data
- `osh_voxel_mat_permatassari2020.h` — embedded Permatassari material data
- `osh_gemca2_voxel_parse.c` / `.h` — legacy `.hed` / `.ctx` voxel parser
- `osh_gemca2_voxel.h` / `osh_gemca2_voxel_defines.h` / `osh_gemca2_voxel_keys.h` — internal structs, constants, and parser keys

## Current role

The new DICOM-backed `DCM` path in `src/apps/osh/` fills voxel bodies directly
from CT studies. This directory provides the HU/material calibration pieces
needed by that workflow and keeps the legacy voxel parser available for future
`VOX` card work.
