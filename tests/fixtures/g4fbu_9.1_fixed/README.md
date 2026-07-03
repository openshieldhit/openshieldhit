G4FermiBreakUp reference histograms (ROOT format)
==================================================

Original output of the standalone Geant4 de-excitation test FermiTest by
Igor Pshenichnov (2006; Geant4 9.1 "fixed" build), kindly provided by the
author and used with his permission.

One file per parent nucleus (C-12, C-13, N-12, N-13), decayed at rest with
E* sampled uniformly over 0..10 MeV/nucleon:
  multip  TH1D  mean fragment multiplicity vs E*/A
  Z, A    TH2D  mean per-event fragment charge/mass yields vs E*/A
  C10/C11 TH1D  specific isotope yields (carbon parents)

Plain-text extractions of these histograms are committed in
examples/05_fermi_breakup_validation/g4fbu_9.1_fixed/ and can be regenerated
without a ROOT installation:

  .venv/bin/pip install uproot
  .venv/bin/python examples/05_fermi_breakup_validation/extract_g4_reference.py \
      tests/fixtures/g4fbu_9.1_fixed examples/05_fermi_breakup_validation/g4fbu_9.1_fixed
