NUCRE isolation benchmark — 200 MeV protons in water, NUCRE 0 (nuclear off)
==========================================================================

Part of the issue #212 set (idd_water_200mev_nucre0/1/2/3).  These decks probe
the nuclear-reaction channel in isolation, the same way the #133 scat* decks
isolate multiple scattering and the #190 strag* decks isolate energy straggling.

Shared phantom / beam / scoring (all four cases)
------------------------------------------------
- 200 MeV proton pencil beam, BEAMPOS z = -5 cm, water RCC R = 20 cm, L = 28 cm
  (identical geometry to the scat/strag sets), NSTAT = 200000.
- MSCAT 0 and STRAGG 0: primaries follow a clean path with no lateral escape from
  the 1 cm^2 central column, so the nuclear channel is the ONLY source of
  secondaries.  The species decomposition and the plateau spectrum therefore
  isolate the nuclear contribution.
- detect.dat scores, on a 1 cm^2 x 0.5 mm central column (DDCnarrow):
    * total Dose (col 4 = the harness auto-compared observable, issue #212);
    * Dose and Fluence split by species — primary protons, all protons, alphas,
      and heavy recoils/fragments (Z >= 3, i.e. C/O recoils);
  and, in a thin mid-plateau slab (Plateau, z = 9.5..10.5 cm), a differential
  secondary spectrum dPhi/dEkin vs Ekin (0.1..300 MeV, 150 log bins) for all
  particles and for protons.

The four NUCRE modes
--------------------
- nucre0: NUCRE 0 — nuclear off (this case).  Baseline: no secondaries, so all
  species pages except the primary/all-proton curves are ~zero.
- nucre1: NUCRE 1 — inelastic (Tripathi) + pp-elastic; the physical "on" case.
- nucre2: NUCRE 2 — elastic only.  OSH-only diagnostic (SH12A has no such mode).
- nucre3: NUCRE 3 — inelastic only.  OSH-only diagnostic.
Modes 2 and 3 decompose mode 1 into its elastic and inelastic contributions.

Reference / comparison
----------------------
NUCRE 0/1 have SH12A mirror decks under
../shieldhit/idd_water_200mev_nucre{0,1}; modes 2/3 do not (SH12A only supports
NUCRE 0/1).  The CTest entry (label "reference") is SKIPPED until a reference
curve is present in reference/.  Tolerances via compare_args.cmake.

Overlay all four modes (OSH) against the SH12A mirrors (0/1):
  python3 tools/plot_nucre.py
