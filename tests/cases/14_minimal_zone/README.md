# 14_minimal_zone — end-to-end Zone scoring

150 MeV protons into a water target, scored with a **`Geometry Zone`** detector
(the Zone geometry's first end-to-end case; #245 acceptance criterion #5, #251).

Same beam and material as `00_minimal`, but the water target is split into two
CSG zones so the Zone output has two physically distinct non-empty bins:

| Zone (geo.dat) | Region | Role |
|---|---|---|
| `001` | water, z 0–10 cm | entrance plateau |
| `002` | water, z 10–20 cm | contains the ~15.8 cm Bragg peak |
| `003` | vacuum shell | beam entry (fluence, ~no energy) |
| `004` | blackhole shell | absorber (empty) |

`beam.dat` disables multiple scattering and straggling (`MSCAT 0`, `STRAGG 0`,
and `NUCRE 0`) so the zone-integrated result is deterministic and stable across
platforms — the reference in `expected/NB_zone.dat` is compared shape-wise by
`compare_dat.py` (score columns L1-normalised, 2 % tolerance). `coord_cols.txt`
tells the driver the single leading column is the zone id, not an X/Y/Z triple.

The run exercises the whole Zone pipeline: `detect.dat` parse → zone-name
resolution against `geo.dat` → compile → per-step Zone deposit → `postprocess`
÷volume → TEXT + BDO2019 save. Energy conservation is visible in the reference:
`ENERGY(001) + ENERGY(002) ≈ 150 MeV` per primary.
