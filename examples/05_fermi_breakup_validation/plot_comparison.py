#!/usr/bin/env python3
"""Plot the openshieldhit Fermi break-up model against the G4FermiBreakUp
reference curves (Geant4 9.1 fixed; I. Pshenichnov's FermiTest, used with
his kind permission).

Runs the fbu_scan example binary for each nuclide, then produces:
  - fbu_multiplicity.png : mean fragment multiplicity vs E*/A (4 panels)
  - fbu_zyield_C12.png   : C-12 per-event charge yields vs E*/A

Usage:
    .venv/bin/python plot_comparison.py [--bin <path-to-fbu_scan>] [--outdir DIR]
"""

import argparse
import os
import subprocess
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
NUCLIDES = [("C12", 6, 12), ("C13", 6, 13), ("N12", 7, 12), ("N13", 7, 13)]


def load_table(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if s and not s.startswith("#"):
                rows.append([float(c) for c in s.split()])
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bin", default=os.path.join(HERE, "../../build_debug/bin/fbu_scan"), help="fbu_scan binary")
    p.add_argument("--outdir", default=HERE, help="output directory for .dat and .png files")
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # ---- run the model scan for each nuclide --------------------------------
    osh = {}
    for name, z, a in NUCLIDES:
        out = os.path.join(args.outdir, f"{name}_osh.dat")
        with open(out, "w") as fh:
            subprocess.run([args.bin, str(z), str(a)], stdout=fh, check=True)
        osh[name] = load_table(out)

    g4 = {name: load_table(os.path.join(HERE, "g4fbu_9.1_fixed", f"{name}_multiplicity.dat")) for name, _, _ in NUCLIDES}
    g4z = load_table(os.path.join(HERE, "g4fbu_9.1_fixed", "C12_zyield.dat"))

    # ---- figure 1: multiplicity, 2x2 panels ----------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True, sharey=True)
    for ax, (name, z, a) in zip(axes.flat, NUCLIDES):
        ax.plot([r[0] for r in g4[name]], [r[1] for r in g4[name]], "k-", lw=1.5, label="G4FermiBreakUp 9.1")
        ax.plot([r[0] for r in osh[name]], [r[1] for r in osh[name]], "r--", lw=1.5, label="osh statistical FBU")
        ax.set_title(f"{name} (Z={z}, A={a})")
        ax.grid(alpha=0.3)
    axes[0, 0].legend(fontsize=9)
    for ax in axes[1, :]:
        ax.set_xlabel("E* / A  [MeV/nucleon]")
    for ax in axes[:, 0]:
        ax.set_ylabel("mean fragment multiplicity")
    fig.suptitle("Fermi break-up: mean multiplicity vs excitation (parent at rest)")
    fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, "fbu_multiplicity.png"), dpi=130)
    print("wrote", os.path.join(args.outdir, "fbu_multiplicity.png"))

    # ---- figure 2: C-12 charge yields ----------------------------------------
    fig2, ax = plt.subplots(figsize=(9, 6))
    colors = {0: "tab:gray", 1: "tab:blue", 2: "tab:red", 3: "tab:green", 6: "tab:purple"}
    for zz, col in colors.items():
        ax.plot([r[0] for r in g4z], [r[1 + zz] for r in g4z], "-", color=col, lw=1.5, label=f"Z={zz} G4")
        ax.plot([r[0] for r in osh["C12"]], [r[2 + zz] for r in osh["C12"]], "--", color=col, lw=1.5, label=f"Z={zz} osh")
    ax.set_xlabel("E* / A  [MeV/nucleon]")
    ax.set_ylabel("mean fragments per event")
    ax.set_title("C-12 Fermi break-up: per-event charge yields")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, ncol=2)
    fig2.tight_layout()
    fig2.savefig(os.path.join(args.outdir, "fbu_zyield_C12.png"), dpi=130)
    print("wrote", os.path.join(args.outdir, "fbu_zyield_C12.png"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
