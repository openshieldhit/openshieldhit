SH12A gold-standard fixture set: plan01 geoA SOBPcent
=====================================================

Curated from:
  https://github.com/APTG/2022_DCPT_LET/tree/main/data/sh12a/results/plan01_field01_geoA_SOBPcent

This directory intentionally keeps only the 32 ASCII `*.dat` scorer outputs
used as manual-reference fixtures. The SH12A input files, PNG plots, and scalar
`NB_target*.txt` files stay outside the repository.

Matching OpenShieldHIT case:

- `tests/cases/10_plan01_geoA_sobpcent`

## File inventory

Depth curves on the narrow Z mesh:

| File | Meaning |
|------|---------|
| `NB_Z_narrow_dose_p1.dat` | Fluence, all particles |
| `NB_Z_narrow_dose_p2.dat` | Dose in transport medium, all particles |
| `NB_Z_narrow_dose_p3.dat` | Dose in transport medium, all protons (primary + secondary) |
| `NB_Z_narrow_dose_p4.dat` | Fluence, primary protons only |
| `NB_Z_narrow_dose_p5.dat` | Fluence, all protons (primary + secondary) |
| `NB_Z_narrow_dose_water_p1.dat` | Fluence, all particles |
| `NB_Z_narrow_dose_water_p2.dat` | Dose to water, all particles |
| `NB_Z_narrow_dose_water_p3.dat` | Dose to water, all protons (primary + secondary) |
| `NB_Z_narrow_LET_p1.dat` | Dose-averaged LET, all particles |
| `NB_Z_narrow_LET_p2.dat` | Dose-averaged LET, primary protons only |
| `NB_Z_narrow_LET_p3.dat` | Dose-averaged LET, all protons (primary + secondary) |
| `NB_Z_narrow_LET_p4.dat` | Track-averaged LET, all particles |
| `NB_Z_narrow_LET_p5.dat` | Track-averaged LET, primary protons only |
| `NB_Z_narrow_LET_p6.dat` | Track-averaged LET, all protons (primary + secondary) |
| `NB_Z_narrow_LET_water_p1.dat` | Dose-averaged LET to water, all particles |
| `NB_Z_narrow_LET_water_p2.dat` | Dose-averaged LET to water, primary protons only |
| `NB_Z_narrow_LET_water_p3.dat` | Dose-averaged LET to water, all protons (primary + secondary) |
| `NB_Z_narrow_LET_water_p4.dat` | Track-averaged LET to water, all particles |
| `NB_Z_narrow_LET_water_p5.dat` | Track-averaged LET to water, primary protons only |
| `NB_Z_narrow_LET_water_p6.dat` | Track-averaged LET to water, all protons (primary + secondary) |
| `NB_Z_narrow_QEFF_p1.dat` | Dose-averaged QEFF, all particles |
| `NB_Z_narrow_QEFF_p2.dat` | Dose-averaged QEFF, primary protons only |
| `NB_Z_narrow_QEFF_p3.dat` | Dose-averaged QEFF, all protons (primary + secondary) |
| `NB_Z_narrow_QEFF_p4.dat` | Track-averaged QEFF, all particles |
| `NB_Z_narrow_QEFF_p5.dat` | Track-averaged QEFF, primary protons only |
| `NB_Z_narrow_QEFF_p6.dat` | Track-averaged QEFF, all protons (primary + secondary) |

Differential target scorers:

| File | Meaning |
|------|---------|
| `NB_target_diff_p1.dat` | Differential fluence in target, all particles, scored versus DEDX |
| `NB_target_diff_p2.dat` | Differential fluence in target, primary protons only, scored versus LET |
| `NB_target_diff_p3.dat` | Differential fluence in target using Si stopping power, all particles, scored versus DEDX |
| `NB_target_diff_p4.dat` | Differential fluence in target using Si stopping power, primary protons only, scored versus DEDX |
| `NB_target_water_diff_p1.dat` | Differential fluence in target using water stopping power, all particles, scored versus DEDX |
| `NB_target_water_diff_p2.dat` | Differential fluence in target using water stopping power, primary protons only, scored versus LET |

## Notes

- The same file inventory applies to `plan02_geoD_mono`.
- See `tests/reference/shieldhit/README.md` for the shared overview and known
  comparison caveats, especially for all-particle dose.
