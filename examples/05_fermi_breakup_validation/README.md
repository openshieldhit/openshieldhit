# Fermi break-up validation against G4FermiBreakUp

Validates the openshieldhit Fermi break-up model
(`src/physics/nuclear/osh_nuclear_fermi_breakup.c`, a microcanonical
statistical break-up: all N=2..6 ground-state partitions weighted by their
phase-space volume) against the canonical **G4FermiBreakUp** implementation.

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
cmake --build --preset debug --parallel --target fbu_scan
.venv/bin/python examples/05_fermi_breakup_validation/plot_comparison.py
```

This runs `fbu_scan <Z> <A>` (parent at rest, 100 E\* bins × 2000 events) for
all four nuclides and writes `fbu_multiplicity.png` and `fbu_zyield_C12.png`
next to this README.

The scan is **deterministic**: `fbu_scan` seeds PCG32 with a fixed value (4242)
and sweeps fixed bins, so the committed `*_osh.dat` tables are byte-stable and
the PNGs are visually stable across runs.  The `*_osh.dat` files and plots are
committed deliberately as a **visual regression reference** — regenerate them
after touching the model and a non-empty diff means the physics actually moved,
not RNG noise.

## Model

Each prefragment `(Z, A, E*)` is broken into `N = 2 .. 6` ground-state
fragments drawn from the full isotope database.  Partition `i` is sampled with
the microcanonical Fermi break-up weight

```
W_i ∝ [V/((2π)^{3/2} ħ³)]^{N−1} · (∏g_k / ∏n_j!) · (∏m_k)^{3/2}
        / Γ(3N/2−3/2) · (E*+Q_i)^{3N/2−5/2}
```

The leading free-volume / density-of-states factor (`V = (4/3)πr0³A`) sets the
multiplicity scale and makes the different-N weights dimensionally
commensurable.  `r0` is the **sole calibration knob**.  Two-body partitions use
exact two-body kinematics; `N ≥ 3` use Kopylov phase space; particle-unstable
products (He-5, Li-5, Be-8) are decayed further on a work stack.

## Calibration & results (2026-06)

`r0 = 0.50 fm` was fitted by minimising the RMS multiplicity deviation from
G4FermiBreakUp over all four nuclides (RMS ≈ 0.37, down from ≈ 1.0 for the old
sequential-binary model).  The effective `r0` is smaller than the geometric
nuclear radius (~1.2 fm) because the **ground-state-only** approximation routes
all of `E*+Q` into kinetic energy and would otherwise over-favour high
multiplicity; the smaller free volume absorbs that bias.

- **E\*/A ≳ 3 MeV/nucleon**: tracks G4 closely (within a few tenths of a unit)
  all the way to 10 MeV/nucleon, replacing the old ≈ 2.8 plateau.
- **Low E\* (~1–2 MeV/nucleon)**: a modest **overshoot** where a single
  high-N channel (e.g. C-12 → 3α) opens — its large phase space dominates
  because no energy is absorbed into fragment excitation.
- **Highest E\* (~10 MeV/nucleon)**: a small **undershoot** (osh ≈ 3.9 vs
  G4 ≈ 4.4 for C-12) — the tail that excited-state channels would fill.

Both residuals are the expected signature of the ground-state-only
approximation (see follow-up below), not a calibration deficiency.  Three
anchor points (1.05 / 5 / 10 MeV/nucleon) are pinned as unit tests in
`tests/unit/test_osh_nuclear_fermi_breakup.c` to guard the calibration.

## Regenerating the reference tables

The original ROOT files are committed in `tests/fixtures/g4fbu_9.1_fixed/`
(16 kB each).  The text tables can be re-extracted without a ROOT
installation:

```sh
.venv/bin/pip install uproot
.venv/bin/python extract_g4_reference.py ../../tests/fixtures/g4fbu_9.1_fixed g4fbu_9.1_fixed
```

## Known follow-ups

- **Excited fragment states (#196).**  Fragments are currently produced in their
  ground state only.  Adding systematic particle-unstable excited-state level
  tables would (a) fill the residual high-E\* multiplicity tail, (b) remove the
  low-E\* overshoot, and (c) let `r0` move back toward a physical ~1.2 fm.
- The same test suite contains SMM material (`MultiFragTestHisto`, Fortran
  SMM outputs for Pb-208) — relevant once heavy-residue de-excitation (SMM)
  is implemented.
