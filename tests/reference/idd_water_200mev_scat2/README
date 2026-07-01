MCS distal-edge benchmark — 200 MeV protons in water, MSCAT 2 (Moliere)
========================================================================

Part of the issue #133 set (idd_water_200mev_scat0/scat1/scat2). This is the case
the issue is about: the primary-proton fluence distal edge diverges from SH12A under
MSCAT 2 (Moliere) while matching with STRAGG 1. Pencil beam of 200 MeV protons on a
water cylinder (R = 20 cm, L = 28 cm; CSDA range ~25.9 cm). NUCRE off and STRAGG off,
so transport is pure EM stopping + multiple scattering: the distal-edge shape is
governed by MCS alone.

  scat0: MSCAT 0  (no scatter; ideal sharp-edge baseline, validates the scorer)
  scat1: MSCAT 1  (Gaussian)
  scat2: MSCAT 2  (Moliere; the case the issue is about)

Scoring (detect.dat), primary protons only (Filter GEN=0 Z=1 A=1):
  - idd.dat      1 cm^2 central depth column, 0.5 mm bins (fluence + energy).
                 This is the file the CTest "reference" harness auto-compares.
  - ddc_wide.dat full-width laterally-integrated depth curve (path-length/range
                 bias only, no lateral escape).
  - rad_cyl.dat  cylindrical radial profile, 25 radial bins x 14 coarse depth slabs.

Reference setup (SH12A / FLUKA):
  - identical phantom, beam, stopping-power table (Water.txt), and physics
    switches; beam/geo/mat decks are SH12A-compatible as-is, only the scoring
    deck needs translating;
  - score primary-proton fluence on the same narrow 1 cm^2 depth mesh;
  - export per-primary values as text: two columns (z [cm], fluence per primary
    per bin) or openshieldhit mesh format (X Y Z F ...);
  - save as reference/idd_sh12a.dat (or idd_fluka.dat, ...).

The CTest entry (label "reference") is SKIPPED until a reference curve is present.
Tolerances can be tuned per case via compare_args.cmake.
