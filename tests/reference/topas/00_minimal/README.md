TEST 00_minimal — 150 MeV protons in water, electromagnetic physics only
21.05.2026 / Niels Bassler

OpenTOPAS cross-validation reference case for tests/cases/00_minimal.
Physics: g4em-standard_opt3 only — no hadronic modules loaded, so nuclear
reactions are effectively OFF. Primary fluence should be flat with depth.

Run:
  cd tests/reference/topas/00_minimal
  topas run.txt

Output: NB_energy.csv, NB_fluence.csv (one value per Z bin, per primary)

Plot comparison with OSH and SH12A:
  python3 tools/plot_nucre.py
