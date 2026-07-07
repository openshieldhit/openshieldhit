#!/usr/bin/env python3
"""Overlay nucre_scan production spectra for several configurations.

Stage-0 analysis of the fast nuclear reaction stage (issues #221 / #260):
takes one or more nucre_scan output files (typically the same reaction at
different compile-time knob values, e.g. OSH_ABRASION_EXCITATION_PER_HOLE_MEV,
or the same knob at different beam energies), prints a metrics table
(per-species yields, mean energies, prefragment E* statistics parsed from the
file headers), and writes a figure overlaying the per-species production
spectra and the prefragment E* distribution.

Usage:
  plot_scan.py --out scan.pdf LABEL1=file1.dat LABEL2=file2.dat ...
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/openshieldhit-matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SPECIES = ("n", "p", "d", "t", "he3", "alpha")
SPECIES_LABEL = {"n": "n", "p": "p", "d": "d", "t": "t", "he3": "He-3", "alpha": "alpha"}


def parse_scan(path: Path) -> dict:
    """Parse one nucre_scan output: header metrics + spectra table."""
    meta: dict = {"path": path, "yield": {}, "mean_e": {}}
    rows = []
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("#"):
                m = re.match(r"#\s+(n|p|d|t|he3|alpha)\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)$", s)
                if m:
                    meta["yield"][m.group(1)] = float(m.group(2))
                    meta["mean_e"][m.group(1)] = float(m.group(3))
                    continue
                m = re.search(r"T_lab=([\d.eE+-]+) MeV\s+sigma_inel=([\d.eE+-]+) mb\s+events=(\d+)", s)
                if m:
                    meta["t_lab"] = float(m.group(1))
                    meta["sigma_mb"] = float(m.group(2))
                    meta["events"] = int(m.group(3))
                m = re.search(r"Estar_mean=([\d.eE+-]+)\s+Estar_std=([\d.eE+-]+)\s+Estar_max=([\d.eE+-]+)", s)
                if m:
                    meta["estar_mean"] = float(m.group(1))
                    meta["estar_std"] = float(m.group(2))
                    meta["estar_max"] = float(m.group(3))
                m = re.search(r"leftover_fragments/event=([\d.eE+-]+)\s+leftover_Estar/event=([\d.eE+-]+)", s)
                if m:
                    meta["leftover_frag"] = float(m.group(1))
                    meta["leftover_estar"] = float(m.group(2))
                continue
            if s:
                rows.append([float(c) for c in s.split()])
    d = np.array(rows)
    meta["elo"], meta["ehi"] = d[:, 0], d[:, 1]
    meta["center"] = np.sqrt(d[:, 0] * d[:, 1])
    meta["width"] = d[:, 1] - d[:, 0]
    for i, sp in enumerate(SPECIES):
        meta[sp] = d[:, 2 + i]
    meta["estar"] = d[:, 8]
    return meta


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+", help="LABEL=path pairs of nucre_scan outputs")
    ap.add_argument("--out", type=Path, default=Path("nucre_scan.pdf"), help="output figure path")
    ap.add_argument("--title", default="nucre_scan production spectra", help="figure title")
    args = ap.parse_args()

    data = {}
    for spec in args.runs:
        label, _, path = spec.partition("=")
        if not path:
            ap.error(f"expected LABEL=path, got '{spec}'")
        data[label] = parse_scan(Path(path))

    # Metrics table on stdout.
    hdr = f"{'config':<16}" + "".join(f"  Y({s})/ev  <E({s})>" for s in ("alpha", "d", "t", "he3")) + "   E*mean  E*std"
    print(hdr)
    for label, m in data.items():
        cells = "".join(
            f"  {m['yield'].get(s, 0.0):8.4f}  {m['mean_e'].get(s, 0.0):7.3f}" for s in ("alpha", "d", "t", "he3")
        )
        print(f"{label:<16}{cells}   {m.get('estar_mean', 0.0):6.2f}  {m.get('estar_std', 0.0):5.2f}")

    fig, axes = plt.subplots(2, 3, figsize=(15, 8), sharex=True)
    panels = list(SPECIES[2:]) + ["p_n", "estar"]  # d, t, he3, alpha, nucleons, E*
    for ax, panel in zip(axes.flat, panels):
        for i, (label, m) in enumerate(data.items()):
            color = f"C{i}"
            if panel == "p_n":
                ax.plot(m["center"], m["p"] / m["width"], color=color, lw=1.2, label=f"{label} p")
                ax.plot(m["center"], m["n"] / m["width"], color=color, lw=1.2, ls="--", label=f"{label} n")
                ax.set_title("nucleons at emission")
            elif panel == "estar":
                ax.plot(m["center"], m["estar"] / m["width"], color=color, lw=1.2, label=label)
                ax.set_title("prefragment E* (before break-up)")
            else:
                ax.plot(m["center"], m[panel] / m["width"], color=color, lw=1.2, label=label)
                ax.set_title(f"{SPECIES_LABEL[panel]} at emission")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=6)
    for ax in axes[1]:
        ax.set_xlabel("Ekin (MeV)")
    for ax in axes[:, 0]:
        ax.set_ylabel("dN/dE per event (1/MeV)")
    fig.suptitle(args.title)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(args.out)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
