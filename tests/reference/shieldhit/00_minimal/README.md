TEST 00_minimal — 150 MeV protons in water, no nuclear reactions
09.04.2026 / Niels Bassler

Baseline case: electromagnetic physics only (stopping, MCS, straggling).
Primary fluence is constant along depth — nothing removes primaries.

Run:
  build_rel/bin/openshieldhit -v tests/cases/00_minimal/

Compare with the tests/reference/idd_water_200mev_nucre* set to see the effect of nuclear reactions.
