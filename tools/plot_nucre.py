#!/usr/bin/env python3
"""Compare dose and primary fluence for cases 00_minimal and 06_minimal_nucre.

Optionally overlays SHIELD-HIT12A reference results when
tests/reference/shieldhit/{00,06}_minimal*/NB_msh.dat are present.

Usage:
    python3 tools/plot_nucre.py [-o output.png]

Requires numpy and matplotlib.
"""

import argparse
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

CASE00 = "tests/cases/00_minimal/NB_msh.dat"
CASE06 = "tests/cases/06_minimal_nucre/NB_msh.dat"
CASE00_SH = "tests/reference/shieldhit/00_minimal/NB_msh.dat"
CASE06_SH = "tests/reference/shieldhit/06_minimal_nucre/NB_msh.dat"


def load_osh(path):
    """Load OpenShieldHIT NB_msh.dat — columns X Y Z energy fluence."""
    try:
        d = np.loadtxt(path, comments="#")
    except OSError:
        print(f"Cannot read {path} — run the case first.", file=sys.stderr)
        sys.exit(1)
    return d[:, 2], d[:, 3], d[:, 4]  # z, energy, fluence


def load_shieldhit(path):
    """Load SHIELD-HIT12A NB_msh.dat — columns Z energy fluence."""
    d = np.loadtxt(path, comments="#")
    return d[:, 0], d[:, 1], d[:, 2]  # z, energy, fluence


def version_string():
    try:
        return (
            "OpenShieldHIT "
            + subprocess
            .check_output(
                ["git", "describe", "--tags", "--dirty", "--always"],
                stderr=subprocess.DEVNULL,
            )
            .decode()
            .strip()
        )
    except Exception:
        return "OpenShieldHIT"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-o", metavar="FILE", help="save to file instead of displaying")
    args = parser.parse_args()

    z00, e00, f00 = load_osh(CASE00)
    z06, e06, f06 = load_osh(CASE06)

    has_sh00 = Path(CASE00_SH).exists()
    has_sh06 = Path(CASE06_SH).exists()
    if has_sh00:
        z00sh, e00sh, f00sh = load_shieldhit(CASE00_SH)
    if has_sh06:
        z06sh, e06sh, f06sh = load_shieldhit(CASE06_SH)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
    fig.suptitle(f"Nuclear reaction effect — 150 MeV protons in water\n{version_string()}", fontsize=11)

    ax1.plot(z00, e00, color="tab:blue", linestyle="-", label="OSH 00 — no nuclear (NUCRE=0)")
    if has_sh00:
        ax1.plot(z00sh, e00sh, color="tab:blue", linestyle="--", label="SH12A 00 — no nuclear (NUCRE=0)")

    ax1.plot(z06, e06, color="tab:orange", linestyle="-", label="OSH 06 — nuclear enabled (NUCRE=1)")
    if has_sh06:
        ax1.plot(z06sh, e06sh, color="tab:orange", linestyle="--", label="SH12A 06 — nuclear enabled (NUCRE=1)")
    ax1.set_ylabel("Energy deposit per primary (MeV)")
    ax1.legend()
    ax1.grid(True)

    ax2.plot(z00, f00, color="tab:blue", linestyle="-", label="OSH 00 — no nuclear (NUCRE=0)")
    if has_sh00:
        ax2.plot(z00sh, f00sh, color="tab:blue", linestyle="--", label="SH12A 00 — no nuclear (NUCRE=0)")

    ax2.plot(z06, f06, color="tab:orange", linestyle="-", label="OSH 06 — nuclear enabled (NUCRE=1)")
    if has_sh06:
        ax2.plot(z06sh, f06sh, color="tab:orange", linestyle="--", label="SH12A 06 — nuclear enabled (NUCRE=1)")
    ax2.set_xlabel("Depth (cm)")
    ax2.set_ylabel("Fluence per primary (cm⁻²)")
    ax2.legend()
    ax2.grid(True)

    plt.tight_layout()
    if args.o:
        plt.savefig(args.o, dpi=150)
        print(f"Saved to {args.o}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
