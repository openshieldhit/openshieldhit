# NUCRE production-spectrum scan

Stage-0 instrument of the fast nuclear reaction stage work
([#221](https://github.com/openshieldhit/openshieldhit/issues/221),
[#260](https://github.com/openshieldhit/openshieldhit/issues/260)).

`nucre_scan` drives the inelastic branch of the nuclear handler — abrasion
followed by Fermi break-up, exactly as wired in `osh_nuclear_handler_step()` —
for p + (Z, A) at a fixed incident energy, with no transport. It tabulates the
**production** (at-emission) observables that transport-level fluence scorers
cannot see directly:

- per-species (n, p, d, t, ³He, α) yields per inelastic event, mean energies,
  and kinetic-energy histograms (150 log bins, 0.1–300 MeV — the same
  *binning and axes* as the NUCRE reference-deck plateau spectra, for
  convenient shape comparison only: production-at-emission and transported
  plateau fluence are different observables, and the quantitative SH12A
  comparison is always deck-vs-deck at transport level);
- the prefragment excitation-energy distribution *before* de-excitation
  (the E\* supply feeding the break-up stage);
- leftover unprocessed fragments and their residual excitation after break-up.

Build with examples enabled and run, e.g. 200 MeV p + ¹⁶O:

```bash
cmake --preset release -DOSH_BUILD_EXAMPLES=ON
cmake --build --preset release --parallel --target nucre_scan
./build/bin/nucre_scan 8 16 200 > p_o16_200mev.dat
```

The E\*-per-hole sensitivity experiment of #260 rebuilds with, e.g.
`-DCMAKE_C_FLAGS="-DOSH_ABRASION_EXCITATION_PER_HOLE_MEV=30.0"` in a scratch
build directory and re-runs the scan; `plot_scan.py` overlays the outputs.
