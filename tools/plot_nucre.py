#!/usr/bin/env python3
"""Compare dose and primary fluence for cases 00_minimal and 06_minimal_nucre.

Optionally overlays SHIELD-HIT12A and OpenTOPAS 4.0 reference results when
their output files are present.

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
CASE00_TP_ENERGY = "tests/reference/topas/00_minimal/NB_energy.csv"
CASE00_TP_FLUENCE = "tests/reference/topas/00_minimal/NB_fluence.csv"
CASE06_TP_ENERGY = "tests/reference/topas/06_minimal_nucre/NB_energy.csv"
CASE06_TP_FLUENCE = "tests/reference/topas/06_minimal_nucre/NB_fluence.csv"

_TOPAS_N_PRIMARIES = 100000
# phantom spans z=0..20 cm; beam enters at z=20 cm (bin N-1), exits at z=0
_TOPAS_PHANTOM_LENGTH_CM = 20.0


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


def load_topas(energy_path, fluence_path):
    """Load OpenTOPAS CSV output and return (depth_cm, energy_MeV, fluence_cm2).

    energy_path: NB_energy.csv  — EnergyDeposit Sum in MeV
    fluence_path: NB_fluence.csv — Fluence Sum in /mm²

    Returns values per primary.  Depth runs from 0 (entrance) to phantom
    length; TOPAS bin 0 is the low-z end (exit side), so we reverse here.
    """
    e_data = np.loadtxt(energy_path, comments="#", delimiter=",")
    f_data = np.loadtxt(fluence_path, comments="#", delimiter=",")
    z_bins = e_data[:, 2].astype(int)
    n_bins = len(z_bins)
    dz = _TOPAS_PHANTOM_LENGTH_CM / n_bins
    # depth from entrance: bin N-1 is entrance (depth≈0), bin 0 is exit
    depth = _TOPAS_PHANTOM_LENGTH_CM - (z_bins + 0.5) * dz
    energy = e_data[:, 3] / _TOPAS_N_PRIMARIES          # MeV per primary
    fluence = f_data[:, 3] / _TOPAS_N_PRIMARIES * 100.0  # /mm² → /cm²
    idx = np.argsort(depth)
    return depth[idx], energy[idx], fluence[idx]


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

    has_tp00 = (Path(CASE00_TP_ENERGY).exists() and Path(CASE00_TP_FLUENCE).exists()
                and Path(CASE00_TP_ENERGY).stat().st_size > 0)
    has_tp06 = (Path(CASE06_TP_ENERGY).exists() and Path(CASE06_TP_FLUENCE).exists()
                and Path(CASE06_TP_ENERGY).stat().st_size > 0)
    if has_tp00:
        z00tp, e00tp, f00tp = load_topas(CASE00_TP_ENERGY, CASE00_TP_FLUENCE)
    if has_tp06:
        z06tp, e06tp, f06tp = load_topas(CASE06_TP_ENERGY, CASE06_TP_FLUENCE)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
    fig.suptitle(f"Nuclear reaction effect — 150 MeV protons in water\n{version_string()}", fontsize=11)

    ax1.plot(z00, e00, color="tab:blue", linestyle="-", label="OSH 00 — no nuclear (NUCRE=0)")
    if has_sh00:
        ax1.plot(z00sh, e00sh, color="tab:blue", linestyle=":", label="SH12A 00 — no nuclear")
    if has_tp00:
        ax1.plot(z00tp, e00tp, color="tab:blue", linestyle="--", label="TOPAS 00 — no nuclear")

    ax1.plot(z06, e06, color="tab:orange", linestyle="-", label="OSH 06 — nuclear (NUCRE=1)")
    if has_sh06:
        ax1.plot(z06sh, e06sh, color="tab:orange", linestyle=":", label="SH12A 06 — nuclear")
    if has_tp06:
        ax1.plot(z06tp, e06tp, color="tab:orange", linestyle="--", label="TOPAS 06 — nuclear")
    ax1.set_ylabel("Energy deposit per primary (MeV)")
    ax1.legend()
    ax1.grid(True)

    ax2.plot(z00, f00, color="tab:blue", linestyle="-", label="OSH 00 — no nuclear (NUCRE=0)")
    if has_sh00:
        ax2.plot(z00sh, f00sh, color="tab:blue", linestyle=":", label="SH12A 00 — no nuclear")
    if has_tp00:
        ax2.plot(z00tp, f00tp, color="tab:blue", linestyle="--", label="TOPAS 00 — no nuclear")

    ax2.plot(z06, f06, color="tab:orange", linestyle="-", label="OSH 06 — nuclear (NUCRE=1)")
    if has_sh06:
        ax2.plot(z06sh, f06sh, color="tab:orange", linestyle=":", label="SH12A 06 — nuclear")
    if has_tp06:
        ax2.plot(z06tp, f06tp, color="tab:orange", linestyle="--", label="TOPAS 06 — nuclear")
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
