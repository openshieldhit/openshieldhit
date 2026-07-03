Straggling distal-edge benchmark — 200 MeV protons in water, STRAGG 1 (Gaussian)
================================================================================

Part of the issue #190 set (idd_water_200mev_strag0/strag1/strag2); see
../idd_water_200mev_strag0/README for the shared phantom, beam, scoring and
reference-setup description.

This case: STRAGG 1 (Gaussian / Bohr energy straggling), MSCAT off, NUCRE off.
The distal Bragg-peak edge is broadened by Gaussian range straggling relative to
the strag0 (off) baseline.

The CTest entry (label "reference") is SKIPPED until a reference curve is
present in reference/. Tolerances via compare_args.cmake.
