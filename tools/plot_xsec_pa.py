#!/usr/bin/env python3
"""Compare OpenShieldHIT's p+A cross sections against reference data (issue #277).

Inputs:
  * committed reference tables tests/reference/xsec/p_{O16,C12}.txt
    (EXFOR (P,NON)/(P,EL) points + ENDF/B-VIII.0 / TENDL-2023 nonelastic
    curves; produced by tools/fetch_pa_xsec.py);
  * the sigma_scan instrument (examples/06_nucre_validation), run live, which
    dumps exactly the sigma_R / sigma_el the transport uses;
  * the committed SH12A nucre1 fixture, from which the *implied* total removal
    rate is extracted as Sigma_eff(z) = -d/dz ln Phi_primary and mapped to
    energy via a Bragg-Kleeman range relation calibrated on the fixture's own
    range.

Figures (PDF pages, or PNGs via --png-dir):
  1. p+O-16: transport sigma_R vs EXFOR / LA150 / TENDL and ratio strip,
     plus the sigma_el curve.
  2. p+C-12: same, with the Garron (P,EL) integrated elastic point against
     OSH's sigma_el curve.
  3. Water removal: Sigma(E) of the 1 cm^2-column primary attenuation —
     SH12A-implied vs OSH model variants (sigma_R only / + sigma_el).

Usage:
    python tools/plot_xsec_pa.py [--out xsec_report.pdf]
                                 [--sigma-scan build/bin/sigma_scan]

To (re)generate the static figures embedded in the physics reference:
    python tools/plot_xsec_pa.py --png-dir docs/physics/img
The PNGs are deterministic (sigma_scan + committed reference data, no Monte
Carlo), so they are committed and regenerated deliberately, like the σ_R
data header.
"""

from __future__ import annotations

import argparse
import os
import subprocess
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/openshieldhit-matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

N_A = 6.02214076e23
RHO_WATER = 1.0  # g/cm^3
N_O_PER_CM3 = RHO_WATER * N_A / 18.015  # oxygen atoms
N_H_PER_CM3 = 2.0 * N_O_PER_CM3
MB_TO_CM2 = 1.0e-27

# Bragg-Kleeman exponent for protons in water (alpha is calibrated on the
# fixture's own range, so only the exponent is assumed).
BK_P = 1.77


def load_reference(path: Path) -> dict[str, np.ndarray]:
    """{source: array[[E_MeV, sigma_mb, err_mb], ...]} from a committed table."""
    blocks: dict[str, list[list[float]]] = {}
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            blocks.setdefault(parts[0], []).append([float(x) for x in parts[1:4]])
    return {k: np.array(v) for k, v in blocks.items()}


def run_sigma_scan(binary: Path, z: int, a: int) -> np.ndarray:
    out = subprocess.run([str(binary), str(z), str(a), "5", "250", "246"], capture_output=True, text=True, check=True).stdout
    rows = [[float(x) for x in ln.split()] for ln in out.splitlines() if ln.strip() and not ln.startswith("#")]
    return np.array(rows)  # E, sigma_R_mb, sigma_el_mb, theta_med


def implied_removal(idd_path: Path) -> tuple[np.ndarray, np.ndarray]:
    """(E_MeV, Sigma_eff_per_cm) implied by a nucre1 primary depth fluence.

    Sigma_eff(z) = -d/dz ln Phi_prim(z), smoothed; z -> E via Bragg-Kleeman
    R = alpha E^p with alpha calibrated so that E(0) = 200 MeV at the file's
    own dose-peak depth.  Works on both layouts: SH12A 1D (Z + 16 pages) and
    OpenShieldHIT mesh (X Y Z + 16 quantities).

    Note this is a *column observable*, not the microscopic event rate: a
    scattered primary keeps being counted until it drifts out of the 1 cm^2
    column, so Sigma_eff lags the event rate near the entrance.  Compare
    implied curves against each other (SH12A vs OSH), and event-rate model
    curves against each other — not across the two kinds.
    """
    rows = []
    with open(idd_path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith(("#", "!", "*")):
                continue
            try:
                rows.append([float(c) for c in s.split()])
            except ValueError:
                continue
    d = np.array(rows)
    if d.shape[1] >= 19:  # OSH mesh layout: X Y Z + 8 dose + 8 fluence
        z, phi, dose = d[:, 2], d[:, 12], d[:, 3]
    else:  # SH12A 1D layout: Z + 16 pages
        z, phi, dose = d[:, 0], d[:, 10], d[:, 1]

    r0 = z[np.argmax(dose)]  # range proxy: Bragg-peak depth
    alpha = r0 / 200.0**BK_P
    keep = z < (r0 - 1.5)  # avoid the distal falloff where Phi -> straggling-dominated
    z, phi = z[keep], phi[keep]

    # ln Phi slope over +-0.5 cm windows (fixture bins are 0.05 cm)
    lnphi = np.log(phi)
    win = 21
    kernel_x = np.arange(win) - win // 2
    slope = np.full(len(z), np.nan)
    half = win // 2
    for i in range(half, len(z) - half):
        y = lnphi[i - half : i + half + 1]
        slope[i] = np.polyfit(kernel_x * (z[1] - z[0]), y, 1)[0]
    ok = ~np.isnan(slope)
    e = ((r0 - z[ok]) / alpha) ** (1.0 / BK_P)
    return e, -slope[ok]


def sigma_water_per_cm(scan_o: np.ndarray, scan_h: np.ndarray, use_elastic: float) -> tuple[np.ndarray, np.ndarray]:
    """Sigma(E) [1/cm] for water from OSH curves; use_elastic scales sigma_el (0..1)."""
    e = scan_o[:, 0]
    sig_o = scan_o[:, 1] + use_elastic * scan_o[:, 2]
    sig_h = np.interp(e, scan_h[:, 0], use_elastic * scan_h[:, 2])
    return e, MB_TO_CM2 * (N_O_PER_CM3 * sig_o + N_H_PER_CM3 * sig_h)


def sigma_r_panel(ax, rax, ref: dict, scan: np.ndarray, label: str) -> None:
    e_osh, sr_osh = scan[:, 0], scan[:, 1]
    ax.plot(e_osh, sr_osh, "k-", lw=1.5, zorder=5, label="OSH $\\sigma_R$ (transport)")
    for src, arr in sorted(ref.items()):
        if src.startswith("exfor:"):
            ax.errorbar(arr[:, 0], arr[:, 1], yerr=arr[:, 2], fmt="o", ms=3, lw=0.8, capsize=2, label=f"EXFOR {src[6:]}")
    for src, style, lbl in (("endfb8-nonel", "C0-", "ENDF/B-VIII.0 (LA150/ICRU63)"), ("tendl23-nonel", "C1--", "TENDL-2023")):
        if src in ref:
            ax.plot(ref[src][:, 0], ref[src][:, 1], style, lw=1.2, label=lbl)
    ax.set(xlim=(0, 250), ylabel="$\\sigma$ (mb)", title=f"{label}: nonelastic / reaction cross section")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=6, ncol=2)

    if "endfb8-nonel" in ref:
        ev = ref["endfb8-nonel"]
        ratio = np.interp(ev[:, 0], e_osh, sr_osh) / ev[:, 1]
        rax.plot(ev[:, 0], ratio, "C0-", lw=1.2, label="OSH / LA150")
    for src, arr in sorted(ref.items()):
        if src.startswith("exfor:"):
            rax.plot(arr[:, 0], np.interp(arr[:, 0], e_osh, sr_osh) / arr[:, 1], "o", ms=3)
    rax.axhline(1.0, color="k", lw=0.5)
    rax.set(xlim=(0, 250), ylim=(0.7, 1.3), ylabel="OSH / ref")
    rax.tick_params(labelbottom=False)
    rax.grid(alpha=0.3)
    rax.legend(fontsize=6)


def build_element_figure(ref: dict, scan: np.ndarray, label: str):
    """σ_R (vs data) + ratio strip + σ_el (vs Garron) for one target element."""
    fig = plt.figure(figsize=(8, 7))
    gs = fig.add_gridspec(3, 1, height_ratios=[3, 1, 2], hspace=0.45)
    sigma_r_panel(fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), ref, scan, label)

    ax = fig.add_subplot(gs[2])
    ax.plot(scan[:, 0], scan[:, 2], "k-", lw=1.5, label="OSH $\\sigma_{el}$ (energy-dependent ratio)")
    for src, arr in sorted(ref.items()):
        if src.startswith("exfor-el:"):
            ax.errorbar(
                arr[:, 0],
                arr[:, 1],
                yerr=arr[:, 2],
                fmt="s",
                ms=5,
                color="C3",
                capsize=3,
                label=f"EXFOR (P,EL) {src[9:]}",
            )
    ax.set(xlim=(0, 250), xlabel="E (MeV)", ylabel="$\\sigma$ (mb)", title="integrated nuclear elastic")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=6)
    return fig


def build_removal_figure(scan_o: np.ndarray, scan_h: np.ndarray, sh12a_idd: Path, osh_idd: Path | None):
    """Primary removal rate Sigma(E) in water: SH12A-implied vs OSH model curves.

    The SH12A-implied curve and the OSH model curves are deterministic (from
    the committed fixture and sigma_scan); the optional OSH-implied overlay is
    from a Monte-Carlo run and is only drawn when @p osh_idd is given.
    """
    fig, ax = plt.subplots(figsize=(8, 5))
    e_sh, sig_sh = implied_removal(sh12a_idd)
    ax.plot(e_sh, sig_sh, "C7.", ms=2, label="SH12A implied (nucre1 fixture, $-d\\ln\\Phi_{prim}/dz$)")
    if osh_idd is not None:
        e_osh_i, sig_osh_i = implied_removal(osh_idd)
        ax.plot(e_osh_i, sig_osh_i, "C2.", ms=2, label="OSH implied (same extraction; compare to SH12A implied)")
    for use_el, style, lbl in (
        (0.0, "C0--", "OSH: $n_O\\sigma_R$ only"),
        (1.0, "C3-", "OSH: $n_O(\\sigma_R+\\sigma_{el}) + n_H\\sigma_{pp}$ (full escape)"),
    ):
        e, sig = sigma_water_per_cm(scan_o, scan_h, use_el)
        if use_el == 0.0:
            sig = MB_TO_CM2 * N_O_PER_CM3 * scan_o[:, 1]
        ax.plot(e, sig, style, lw=1.3, label=lbl)
    ax.set(
        xlim=(0, 210),
        ylim=(0, 0.05),
        xlabel="E (MeV)",
        ylabel="$\\Sigma$ (1/cm)",
        title="Primary removal rate in water, 1 cm$^2$ column (200 MeV deck)",
    )
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7)
    return fig


def main() -> None:
    root = Path(__file__).resolve().parent.parent
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=Path("xsec_report.pdf"), help="multi-page PDF report")
    ap.add_argument("--sigma-scan", type=Path, default=root / "build" / "bin" / "sigma_scan")
    ap.add_argument("--xsec-dir", type=Path, default=root / "tests" / "reference" / "xsec")
    ap.add_argument(
        "--sh12a-idd",
        type=Path,
        default=root / "tests" / "reference" / "shieldhit" / "idd_water_200mev_nucre1" / "idd.dat",
    )
    ap.add_argument(
        "--osh-idd",
        type=Path,
        default=None,
        help="optional OSH nucre1 idd.dat (e.g. from a plot_nucre.py workdir) to overlay its implied removal",
    )
    ap.add_argument(
        "--png-dir",
        type=Path,
        default=None,
        help="also export the individual figures as PNGs into this directory "
        "(e.g. docs/physics/img); the removal figure is the deterministic version",
    )
    args = ap.parse_args()

    ref_o = load_reference(args.xsec_dir / "p_O16.txt")
    ref_c = load_reference(args.xsec_dir / "p_C12.txt")
    scan_o = run_sigma_scan(args.sigma_scan, 8, 16)
    scan_c = run_sigma_scan(args.sigma_scan, 6, 12)
    scan_h = run_sigma_scan(args.sigma_scan, 1, 1)

    # (basename, figure) pairs, built once and reused for both the PDF and PNGs.
    figures = [
        ("xsec_p_O16", build_element_figure(ref_o, scan_o, "p + O-16")),
        ("xsec_p_C12", build_element_figure(ref_c, scan_c, "p + C-12")),
        # PNG export omits the Monte-Carlo overlay so the committed docs figure
        # is reproducible from committed data + sigma_scan alone; the PDF keeps
        # whatever --osh-idd was passed.
        ("xsec_removal_water", build_removal_figure(scan_o, scan_h, args.sh12a_idd, args.osh_idd)),
    ]

    with PdfPages(args.out) as pdf:
        for _name, fig in figures:
            pdf.savefig(fig)
    print(f"wrote {args.out}")

    if args.png_dir is not None:
        args.png_dir.mkdir(parents=True, exist_ok=True)
        for name, fig in figures:
            path = args.png_dir / f"{name}.png"
            fig.savefig(path, dpi=110, bbox_inches="tight", metadata={"Software": None})
            print(f"wrote {path}")

    for _name, fig in figures:
        plt.close(fig)


if __name__ == "__main__":
    main()
