TEST 14_gemca_bodies — 150 MeV protons through an elliptic cylinder + truncated cone

First integration case that exercises the `REC` (elliptic cylinder → `ELLZ`
surface) and `TRC` (truncated cone → `CONE` surface) bodies end-to-end. Before
issue #255 no case or example in the tree referenced these bodies, so the
elliptic-cylinder and cone ray-surface distance math was completely untested by
CTest.

Geometry (all along z):
  - body 1 `REC`  z∈[0,10], ellipse semi-axes 2×3 cm      → zone 001, Water
  - body 2 `TRC`  z∈[10,20], radius 4→2 cm                → zone 002, Water
  - body 3 `RCC`  z∈[-10,40], r=15 cm, minus 1,2          → zone 003, vacuum
  - body 4 `RCC`  z∈[-20,50], r=25 cm, minus 3            → zone 004, blackhole

A pencil beam enters at z=-5 along +z, crosses the elliptic cylinder and the
truncated cone, and stops in the water (150 MeV protons range ≈ 15.8 cm, so the
Bragg peak lands inside the `TRC` region). Electromagnetic only (NUCRE 0).

This is a smoke/coverage case: it validates that transport runs to completion
(exit 0) through these body types. Per-surface distance correctness is locked
by the unit tests test_osh_gemca_runtime_dist and test_osh_gemca2_dist.

Run:
  build_debug/bin/openshieldhit -v tests/cases/14_gemca_bodies/
