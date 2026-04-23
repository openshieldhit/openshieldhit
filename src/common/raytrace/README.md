# common/raytrace

This directory provides the voxel-grid ray traversal layer used by scoring
today and by voxel CT transport work in progress.

## Purpose

- expose a small domain-neutral grid descriptor: `struct osh_raytrace_grid`
- return voxel crossings as flat `(idx, path_len)` pairs
- keep the traversal contract independent of GEMCA, scoring, or DICOM types

## Files

- `osh_raytrace.h` — public interface for the raytrace submodule
- `osh_raytrace_simple_msh.c` — pedagogical 3-D DDA reference
- `osh_raytrace_siddon_msh.c` — Siddon clip + DDA walk, current default
- `osh_raytrace_jacobs_msh.c` — alpha-parametric Jacobs implementation
- `CMakeLists.txt` — compile-time selection of the active implementation

## Notes

- The active algorithm is selected at configure time with
  `-DOSH_RAYTRACE_ALGORITHM=SIMPLE|SIDDON|JACOBS`.
- All implementations honor `grid->tile_order`, so row-major and Morton-tiled
  voxel layouts share the same caller contract.
- This layer only traces regular Cartesian grids. Geometry ownership,
  transforms, HU tables, and transport physics live elsewhere.
