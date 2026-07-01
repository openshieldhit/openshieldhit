SHIELD-HIT12A run decks — 200 MeV protons in water, STRAGG 0 (off)
==================================================================

These test decks are made for SHIELD-HIT12A, used as the external reference code
for openshieldhit.

SHIELD-HIT12A-syntax mirror of tests/reference/idd_water_200mev_strag0 (the
openshieldhit case), part of the issue #190 straggling-isolation set
(strag0/strag1/strag2; MSCAT off, NUCRE off, STRAGG 0/1/2).
Differences from the openshieldhit deck:
  - beam.dat uses JPART0 2 (proton) instead of PRIMARY 1 1;
  - geo.dat assigns materials via the trailing zone->medium map
        1    2    3        (zone numbers)
        1 1000    0        (medium: 1=Water, 1000=vacuum, 0=blackhole)
    instead of ASSIGNMAT;
  - mat.dat uses MEDIUM 1 / ICRU 276 / LOADDEDX / END.
detect.dat and Water.txt are identical to the openshieldhit case (depth Dose,
Fluence Primary, DLET, TLET on the 1 cm^2 central column).

Run, then export the depth-dose curve from the 1 cm^2 central column and drop it
into
  ../../idd_water_200mev_strag0/reference/idd_sh12a.dat
to activate the CTest "reference" comparison.
