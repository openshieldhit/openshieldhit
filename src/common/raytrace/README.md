# common/raytrace

This directory provides the voxel-grid ray traversal layer used by scoring and by voxel
CT transport.

## Concepts

Three distinct levels of abstraction share this directory:

**`osh_raytrace_grid`** — a low-level, domain-neutral descriptor for a regular 3-D grid.
It carries an origin, per-axis voxel spacing, and voxel counts.  Both the Cartesian mesh
and the cylindrical traversal algorithms accept this struct; their field conventions differ
(see below).  Point-location and first-crossing helpers also use it.

**Mesh geometry** (`Geometry Mesh` in detect.dat, `OSH_SCORING_GEO_MESH`) — a Cartesian
(X, Y, Z) scoring grid.  `osh_raytrace_grid` fields map to X/Y/Z origin, spacing, and
counts in the natural way.

**Cyl geometry** (`Geometry Cyl` in detect.dat, `OSH_SCORING_GEO_CYL`) — a cylindrical
(R, Z) scoring grid.  `osh_raytrace_grid` is reused with the following field convention:

| Field index | Mesh meaning | Cyl meaning |
|---|---|---|
| `[0]` | X origin / spacing / count | r_min / dr / nr |
| `[1]` | Y origin / spacing / count | unused (set 0 / 0 / 1) |
| `[2]` | Z origin / spacing / count | z_min / dz / nz |

Flat index for both geometries: `i0 + n[0]*(i1 + n[1]*i2)`.
For Cyl this reduces to `ir + n[0]*iz` since `n[1] = 1`.
Max crossings for Cyl: `2*n[0] + n[2]` (quadratic R-shell intersections give up to
two per R bin, plus one per Z plane).

## Purpose

- Expose a small, domain-neutral grid descriptor: `struct osh_raytrace_grid`.
- Return voxel crossings as `(idx, path_len, vol_inv)` triples via `struct osh_voxel_crossing`.
- Provide O(1) helpers for point location and first-voxel crossing queries.
- Keep the traversal contract independent of GEMCA, scoring, or DICOM types.

## `vol_inv` field

The `vol_inv` field in `osh_voxel_crossing` is **not filled by the traversal functions**.
It is filled by the scoring layer immediately after traversal, using geometry-type-specific
logic (uniform scalar for Mesh; per-R-bin lookup table for Cyl).  This keeps the
traversal functions geometry-agnostic and the scoring accumulation functions volinv-agnostic.

## Files

- `osh_raytrace.h` — public interface: `osh_raytrace_grid`, `osh_voxel_crossing`,
  `osh_raytrace_traverse`, `osh_raytrace_locate`, `osh_raytrace_first_crossing`
- `osh_raytrace_cyl.h` — CYL traversal interface: `osh_raytrace_cyl_traverse`
- `osh_raytrace_grid.c` — common point-location and first-crossing helpers
- `osh_raytrace_simple_msh.c` — pedagogical 3-D DDA reference (Mesh only)
- `osh_raytrace_siddon_msh.c` — Siddon clip + DDA walk, current default (Mesh only)
- `osh_raytrace_jacobs_msh.c` — alpha-parametric Jacobs implementation (Mesh only)
- `osh_raytrace_jacobs_cyl.c` — alpha-parametric Jacobs implementation for Cyl geometry
- `CMakeLists.txt` — compile-time selection of the active Mesh algorithm

## Notes

- The active Mesh algorithm is selected at configure time with
  `-DOSH_RAYTRACE_ALGORITHM=SIMPLE|SIDDON|JACOBS`.
- `osh_raytrace_jacobs_cyl.c` is always compiled in; it is not gated on
  `OSH_RAYTRACE_ALGORITHM`.
- All Mesh implementations honour `grid->tile_order`, so row-major and Morton-tiled
  voxel layouts share the same caller contract.  `tile_order` is ignored by the Cyl
  implementation (always row-major `ir + nr*iz`).
