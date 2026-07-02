SHIELD-HIT12A run decks — 200 MeV protons in water, STRAGG 2 (Vavilov)
======================================================================

SHIELD-HIT12A-syntax mirror of tests/reference/idd_water_200mev_strag2, part of
the issue #190 straggling-isolation set. See ./idd_water_200mev_strag0/README
(../strag0) for the shared syntax differences vs the openshieldhit deck.

This case: STRAGG 2 (Vavilov), MSCAT off, NUCRE off.

SHIELD-HIT12A implements Vavilov straggling, so this deck runs and produces the
reference curve. openshieldhit now implements STRAGG 2 (Vavilov + Landau) too,
so its strag2 case no longer aborts with OSH_ENOTSUP. Generate the reference
curve and drop it into
  ../../idd_water_200mev_strag2/reference/idd_sh12a.dat
so the CTest "reference" comparison (currently SKIPPED for lack of a curve)
becomes an active check.
