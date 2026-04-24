# gemca/voxel

This directory holds legacy `.hed` / `.ctx` voxel-file parser support retained
for future VOX card work.

## Purpose

- legacy voxel file parsing support retained from older code paths

## Files

- `osh_gemca2_voxel_parse.c` / `.h` — legacy `.hed` / `.ctx` voxel parser
- `osh_gemca2_voxel.h` / `osh_gemca2_voxel_defines.h` / `osh_gemca2_voxel_keys.h` — internal structs, constants, and parser keys

## Current role

HU calibration and CT tissue material registration now live in `src/voxel/`.
The DICOM-backed `DCM` path in `src/apps/osh/` fills voxel bodies directly from
CT studies.
