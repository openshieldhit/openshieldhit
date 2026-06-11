#!/usr/bin/env python3
"""Extract G4FermiBreakUp reference curves from Igor Pshenichnov's FermiTest
ROOT files into the plain-text tables committed in g4fbu_9.1_fixed/.

The ROOT files come from the standalone Geant4 de-excitation tests by
Igor Pshenichnov (FermiTest, 2006; Geant4 9.1 "fixed" build), used with his
kind permission.  Each file holds, for one parent nucleus decayed at rest
with E* sampled uniformly in 0..10 MeV/nucleon:
  - multip : TH1D, mean fragment multiplicity vs E*/A
  - Z, A   : TH2D, mean per-event fragment charge/mass yields vs E*/A
  - C10, C11 : TH1D, specific isotope yields (C12/C13 runs)

The original ROOT files are committed in tests/fixtures/g4fbu_9.1_fixed/.
Requires uproot (pip install uproot) — no ROOT installation needed:
    .venv/bin/python extract_g4_reference.py tests/fixtures/g4fbu_9.1_fixed g4fbu_9.1_fixed
"""

import os
import sys

import uproot

NUCLIDES = ["C12", "C13", "N12", "N13"]


def main():
    if len(sys.argv) != 3:
        print(__doc__, file=sys.stderr)
        return 2
    indir, outdir = sys.argv[1], sys.argv[2]
    os.makedirs(outdir, exist_ok=True)

    for nuc in NUCLIDES:
        path = os.path.join(indir, f"{nuc}.root")
        if not os.path.exists(path):
            print(f"skip {nuc}: {path} not found")
            continue
        f = uproot.open(path)

        m = f["multip"]
        edges = m.axis().edges()
        vals = m.values()
        with open(os.path.join(outdir, f"{nuc}_multiplicity.dat"), "w") as fh:
            fh.write(f"# G4FermiBreakUp (Geant4 9.1 fixed, I. Pshenichnov FermiTest) — {nuc} at rest\n")
            fh.write("# E*_per_A_MeV   mean_multiplicity\n")
            for i in range(len(vals)):
                fh.write(f"{0.5 * (edges[i] + edges[i + 1]):8.4f} {vals[i]:12.5f}\n")

        z2 = f["Z"]
        ze = z2.axes[1].edges()
        zv = z2.values()
        with open(os.path.join(outdir, f"{nuc}_zyield.dat"), "w") as fh:
            zlabels = [int(0.5 * (ze[j] + ze[j + 1])) for j in range(zv.shape[1])]
            fh.write(f"# G4FermiBreakUp (9.1 fixed) — {nuc}: mean fragments per event with charge Z\n")
            fh.write("# E*_per_A_MeV  " + "  ".join(f"Z={z}" for z in zlabels) + "\n")
            for i in range(zv.shape[0]):
                fh.write(
                    f"{0.5 * (edges[i] + edges[i + 1]):8.4f} "
                    + " ".join(f"{zv[i, j]:11.5f}" for j in range(zv.shape[1]))
                    + "\n"
                )
        print(f"{nuc}: extracted")
    return 0


if __name__ == "__main__":
    sys.exit(main())
