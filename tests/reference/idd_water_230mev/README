IDD benchmark — 230 MeV protons in water
==========================================

Pencil beam of 230 MeV protons on a water cylinder (R = 20 cm, L = 36 cm),
scoring the laterally integrated energy deposition (IDD) on a 1D depth mesh
(360 bins over 0..36 cm).  Straggling: Gauss; MCS: Moliere; NUCRE on;
NSTAT = 200000; LOADDEDX water stopping power (Water.txt).

Reference setup (SH12A / FLUKA):
  - identical phantom, beam, stopping-power table, and physics switches;
  - score energy deposition on the same depth mesh, laterally integrated
    over the full phantom;
  - export per-primary values as text: two columns (z [cm], energy per
    primary per bin) or openshieldhit mesh format (X Y Z E ...);
  - save as reference/idd_sh12a.dat (or idd_fluka.dat, ...).

The CTest entry (label "reference") is SKIPPED until a reference curve is
present.  Tolerances can be tuned per case via compare_args.cmake, e.g.:
  set(COMPARE_ARGS --tol-integral 0.05 --tol-bin 0.05 --tol-peak-mm 1.0)
