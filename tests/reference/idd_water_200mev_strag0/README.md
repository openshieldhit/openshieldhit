Straggling distal-edge benchmark — 200 MeV protons in water, STRAGG 0 (off)
===========================================================================

Part of the issue #190 set (idd_water_200mev_strag0/strag1/strag2). Pencil beam
of 200 MeV protons on a water cylinder (R = 20 cm, L = 28 cm; CSDA range
~25.9 cm). MSCAT off and NUCRE off, so transport is pure EM stopping + energy
straggling: the distal-edge (Bragg-peak falloff) shape is governed by range
straggling alone.

  strag0: STRAGG 0  (off; ideal sharp-edge baseline, validates the scorer)
  strag1: STRAGG 1  (Gaussian / Bohr)
  strag2: STRAGG 2  (Vavilov / Landau)

Scoring (detect.dat), single 1 cm^2 central depth column, 0.5 mm bins:
  - idd.dat  columns: Dose (auto-compared), Fluence Primary, DLET, TLET.
             Dose is first, so the CTest "reference" harness compares the dose
             distal edge; DLET/TLET are for manual/overlay comparison of the
             LET tail.

Reference setup (SH12A / FLUKA):
  - identical phantom, beam, stopping-power table (Water.txt), and physics
    switches; beam/geo/mat decks are SH12A-compatible as-is, only the scoring
    deck needs translating;
  - score depth dose on the same narrow 1 cm^2 depth mesh;
  - export per-primary values as text: two columns (z [cm], dose per primary
    per bin) or openshieldhit mesh format (X Y Z D ...);
  - save as reference/idd_sh12a.dat (or idd_fluka.dat, ...).

The CTest entry (label "reference") is SKIPPED until a reference curve is present.
Tolerances can be tuned per case via compare_args.cmake.
