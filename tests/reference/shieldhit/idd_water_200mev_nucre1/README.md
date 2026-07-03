SH12A mirror — 200 MeV protons in water, NUCRE 1 (inel + pp-elastic)
===================================================================

SHIELD-HIT12A-syntax mirror of the OpenShieldHIT case
../../idd_water_200mev_nucre1 (issue #212 NUCRE isolation set).  MSCAT 0 +
STRAGG 0, NUCRE 1.  beam.dat uses JPART0; geo.dat carries the numeric
zone->medium map; mat.dat uses MEDIUM.  detect.dat and Water.txt are shared
verbatim with the OSH deck.

Run with the installed shieldhit, then copy the produced idd.dat to
../../idd_water_200mev_nucre1/reference/idd_sh12a.dat to activate the CTest
"reference" comparison:

  cd tests/reference/shieldhit/idd_water_200mev_nucre1 && shieldhit .
  cp idd.dat ../../idd_water_200mev_nucre1/reference/idd_sh12a.dat

Only NUCRE 0/1 have SH12A mirrors; modes 2 (elastic only) and 3 (inelastic only)
are OpenShieldHIT-only.
