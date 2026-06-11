# Fermi break-up validation against G4FermiBreakUp

Validates the openshieldhit Fermi break-up model
(`src/physics/nuclear/osh_nuclear_fermi_breakup.c`, a sequential-binary
development approximation) against the canonical **G4FermiBreakUp**
implementation.

The reference curves in `g4fbu_9.1_fixed/` were extracted from the standalone
Geant4 de-excitation test suite by **Igor Pshenichnov** (FermiTest, 2006;
Geant4 9.1 "fixed" build), kindly provided by the author and used with his
permission.  Each table gives, for a parent nucleus (C-12, C-13, N-12, N-13)
decayed **at rest** with E\* swept over 0–10 MeV/nucleon:

- `<nuc>_multiplicity.dat` — mean fragment multiplicity vs E\*/A
  (an un-broken parent counts as multiplicity 1),
- `<nuc>_zyield.dat` — mean per-event fragment yield for each charge Z.

## Running the comparison

```sh
cmake --build --preset debug --target fbu_scan
.venv/bin/python examples/05_fermi_breakup_validation/plot_comparison.py
```

This runs `fbu_scan <Z> <A>` (parent at rest, 100 E\* bins × 2000 events) for
all four nuclides and writes `fbu_multiplicity.png` and `fbu_zyield_C12.png`
next to this README.

## Results (2026-06)

- **E\*/A ≲ 2 MeV/nucleon** (the region populated by the abrasion stage,
  E\* ≈ 13.3 MeV per knocked-out nucleon): multiplicities agree with
  G4FermiBreakUp to ~1 %; thresholds and the C-12 → 3α channel are exact.
- **E\*/A ≳ 3 MeV/nucleon**: the sequential-binary scheme saturates near
  multiplicity ≈ 2.8 while the canonical simultaneous n-body model rises to
  ≈ 4.3 at 10 MeV/nucleon.  Cause: in a binary split all excitation converts
  to fragment kinetic energy — the products emerge *cold* (ground state), so
  the chain stops early.  The canonical model partitions the parent into n
  bodies (including excited fragment states) simultaneously, reaching
  high-multiplicity final states at high E\*.

Two anchor points from the validated region are pinned as unit tests in
`tests/unit/test_osh_nuclear_fermi_breakup.c` to guard against regressions.

## Regenerating the reference tables

The original ROOT files are committed in `tests/fixtures/g4fbu_9.1_fixed/`
(16 kB each).  The text tables can be re-extracted without a ROOT
installation:

```sh
.venv/bin/pip install uproot
.venv/bin/python extract_g4_reference.py ../../tests/fixtures/g4fbu_9.1_fixed g4fbu_9.1_fixed
```

## Known follow-ups

- Share excitation among binary-split products (or implement true n-body
  partitions with excited fragment states) to fix the high-E\* multiplicity
  deficit.
- The same test suite contains SMM material (`MultiFragTestHisto`, Fortran
  SMM outputs for Pb-208) — relevant once heavy-residue de-excitation (SMM)
  is implemented.
