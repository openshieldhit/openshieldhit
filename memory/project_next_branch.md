---
name: project_next_branch
description: tracks current branch goal and completed milestones
type: project
---

Branch `66-5`: M5 architecture redesign — complete.

**Completed in this session:**

`src/voxel/` refactor — the final piece of M5:

- `include/openshieldhit/voxel.h`: `OSH_HU_TABLE_*` constants, `OSH_VOXEL_HU_LUT_SIZE`
- `src/voxel/` module: `osh_voxel_hu_lut.c` (Schneider), `osh_voxel_hu_lut_permatassari.c`, shared data tables
- `osh_gemca_prepared.hu_table_type` replaces `hu_bin_lut` on cold storage
- `osh_geometry_workspace.hu_table_type` added; set in `osh_run.c` after both parsers run
- `osh_material_workspace` no longer holds LUTs — just `hu_table_type`
- `osh_gemca_compile` builds `hu_bin_lut` from `hu_table_type` (via `osh_voxel`)
- `osh_material_compile` builds `hu_rho_lut` from `hu_table_type` (via `osh_voxel`)
- Parser just sets `hu_table_type` + registers materials; no LUT allocation
- `simulation.c` is fully agnostic — no LUT copies

**Why:** Clean separation — gemca and material don't cross-depend; shared voxel module eliminates the wrong layering of material parser → gemca/voxel.

**All 43 tests pass.**
