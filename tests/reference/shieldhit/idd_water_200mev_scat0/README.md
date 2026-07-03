SHIELD-HIT12A run decks — 200 MeV protons in water, MSCAT 0 (no scatter)
========================================================================

These test decks are made for SHIELD-HIT12A, used as the external reference code
for openshieldhit.

SHIELD-HIT12A-syntax mirror of tests/reference/idd_water_200mev_scat0 (the
openshieldhit case), part of the issue #133 MCS-isolation set (scat0/scat1/scat2).
Differences from the openshieldhit deck:
  - beam.dat uses JPART0 2 (proton) instead of PRIMARY 1 1;
  - geo.dat assigns materials via the trailing zone->medium map
        1    2    3        (zone numbers)
        1 1000    0        (medium: 1=Water, 1000=vacuum, 0=blackhole)
    instead of ASSIGNMAT;
  - mat.dat uses MEDIUM 1 / ICRU 276 / LOADDEDX / END.
detect.dat and Water.txt are identical to the openshieldhit case.

Run, then export the primary-proton fluence depth curve from the 1 cm^2 central
column and drop it into
  ../../idd_water_200mev_scat0/reference/idd_sh12a.dat
to activate the CTest "reference" comparison.
