#!/usr/bin/env python3
"""Overlay OpenShieldHIT vs SHIELD-HIT12A for the issue #212 NUCRE cases.

The cases tests/reference/idd_water_200mev_nucre{0,1,2,3} (200 MeV protons in
water, MSCAT off, STRAGG off) isolate the nuclear-reaction channel: the only
source of secondaries is NUCRE, so a species decomposition and the plateau
secondary spectrum directly expose the nuclear contribution to dose and fluence.

The four NUCRE modes:
  * nucre0  NUCRE 0 - nuclear off (baseline, no secondaries)
  * nucre1  NUCRE 1 - inelastic (Tripathi) + pp-elastic (the physical "on" case)
  * nucre2  NUCRE 2 - elastic only            (OSH-only; SH12A has no such mode)
  * nucre3  NUCRE 3 - inelastic only          (OSH-only)
Modes 2 and 3 decompose mode 1.

OpenShieldHIT is run live (all cases in parallel, at low NSTAT) so the overlay
reflects the current physics.  SHIELD-HIT12A is NOT re-run: its full-statistics
idd.dat / spectrum.dat are committed under tests/reference/shieldhit/ as fixed
reference fixtures (they only change if the geometry/scoring changes).  Only
NUCRE 0/1 have SH12A fixtures; modes 2/3 are plotted OpenShieldHIT-only.

Per case, on the 1 cm^2 central column (idd.dat) and the mid-plateau slab
(spectrum.dat) this writes a multi-page PDF of: depth dose by species, depth
fluence by species, and the plateau secondary spectrum dPhi/dEkin (log-log).
A final summary page overlays all-particle dose and the plateau proton spectrum
across the four modes.

Column layouts:
  idd.dat quantity order (detect.dat): Dose, Dose Primary, Dose Protons,
    Dose Alphas, Dose HeavyRec, Fluence, Fluence Primary, Fluence Protons,
    Fluence Alphas, Fluence HeavyRec
    OpenShieldHIT mesh : X Y Z <10 quantities>   (Z=col2, quantities col3..12)
    SHIELD-HIT12A 1D   : Z <10 quantities>        (Z=col0, quantities col1..10)
  spectrum.dat (Ekin 0.1..300 MeV, 150 log bins, pages: Fluence, Fluence Protons)
    OpenShieldHIT : X Y Z EKIN Phi_all Phi_prot   -- Phi PER BIN (see issue #215)
    SHIELD-HIT12A : Page Diff1bin Value           -- Value already /cm^2/MeV
  Issue #215: OSH differential output is counts-per-bin, not a density, so this
  tool divides the OSH spectrum by the log-bin width; SH12A is already a density.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import math
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/openshieldhit-matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

CASES = ["nucre0", "nucre1", "nucre2", "nucre3"]
MODEL = {
    "nucre0": "NUCRE 0 (off)",
    "nucre1": "NUCRE 1 (inel + pp-el)",
    "nucre2": "NUCRE 2 (elastic only)",
    "nucre3": "NUCRE 3 (inelastic only)",
}
# SHIELD-HIT12A only implements NUCRE 0/1; modes 2/3 are OpenShieldHIT-only.
SH12A_CASES = {"nucre0", "nucre1"}

# idd.dat quantity order (must match detect.dat).
QUANTITIES = (
    "dose", "dose_prim", "dose_prot", "dose_alpha", "dose_heavy",
    "flu", "flu_prim", "flu_prot", "flu_alpha", "flu_heavy",
)

# Differential spectrum axis (must match the Diff1 line in detect.dat).
SPEC_LO, SPEC_HI, SPEC_NBINS = 0.1, 300.0, 150


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parent.parent
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--nstat", type=int, default=20000, help="OpenShieldHIT primaries per run (default 20000)")
    p.add_argument("--repo-root", type=Path, default=root, help="openshieldhit repo root")
    p.add_argument("--osh", type=Path, default=root / "build" / "bin" / "openshieldhit", help="OpenShieldHIT binary")
    p.add_argument("--workdir", type=Path, default=None, help="working dir for runs (default: a temp dir)")
    p.add_argument("--out", type=Path, default=None, help="output PDF (default: <workdir>/nucre_report.pdf)")
    p.add_argument("--cases", nargs="+", default=CASES, choices=CASES, help="cases to include")
    p.add_argument("--no-run", action="store_true", help="reuse existing OpenShieldHIT outputs in --workdir")
    return p.parse_args()


def load_numeric(path: Path) -> np.ndarray:
    """Load whitespace-separated numeric rows, skipping comments / label lines."""
    rows = []
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith(("#", "!", "*")):
                continue
            try:
                rows.append([float(c) for c in s.split()])
            except ValueError:
                continue  # stray axis-label line
    if not rows:
        raise ValueError(f"{path}: no numeric data")
    return np.array(rows)


def log_bin_centers_widths(lo: float, hi: float, nbins: int):
    """Geometric bin centres and widths for a logarithmic axis (matches the OSH scorer)."""
    edges = np.logspace(math.log10(lo), math.log10(hi), nbins + 1)
    centers = np.sqrt(edges[:-1] * edges[1:])
    widths = np.diff(edges)
    return centers, widths


def set_nstat(case_dir: Path, dst: Path, nstat: int) -> None:
    dst.mkdir(parents=True, exist_ok=True)
    for f in case_dir.iterdir():
        if f.is_file():
            shutil.copy(f, dst / f.name)
    beam = dst / "beam.dat"
    text = beam.read_text().splitlines()
    out = []
    for ln in text:
        if ln.lstrip().startswith("NSTAT"):
            out.append(f"NSTAT      {nstat}    -1")
        else:
            out.append(ln)
    beam.write_text("\n".join(out) + "\n")


def format_elapsed(seconds):
    if seconds is None:
        return "n/a"
    if seconds < 60.0:
        return f"{seconds:.2f} s"
    minutes, seconds = divmod(seconds, 60.0)
    return f"{int(minutes)}m {seconds:.1f}s"


def _run(cmd, cwd=None):
    """Run a command; return (ok, last_stderr_line, elapsed_s)."""
    start = time.perf_counter()
    r = subprocess.run(cmd, cwd=cwd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    elapsed_s = time.perf_counter() - start
    if r.returncode == 0:
        return True, "", elapsed_s
    errs = [ln.strip() for ln in r.stderr.splitlines() if ln.strip()]
    for ln in errs:
        if "not implemented" in ln.lower():
            return False, ln.split("]", 1)[-1].strip(), elapsed_s
    return False, (errs[-1] if errs else f"exit {r.returncode}"), elapsed_s


def run_osh(args, case: str):
    src = args.repo_root / "tests" / "reference" / f"idd_water_200mev_{case}"
    deck = args.workdir / "osh_deck" / case
    out = args.workdir / "osh_out" / case
    set_nstat(src, deck, args.nstat)
    out.mkdir(parents=True, exist_ok=True)
    ok, err, elapsed_s = _run([str(args.osh), "-v", "--outdir", str(out), str(deck)])
    return (out, "", elapsed_s) if ok else (None, err, elapsed_s)


def sh12a_ref_dir(args, case: str) -> Path:
    return args.repo_root / "tests" / "reference" / "shieldhit" / f"idd_water_200mev_{case}"


def depth_quantities(out_dir: Path, code: str):
    """Return dict {z, <QUANTITIES>} for the 1D idd.dat depth file, or None."""
    path = out_dir / "idd.dat"
    if not path.exists():
        return None
    d = load_numeric(path)
    if code == "osh":  # X Y Z <10 quantities>
        z = d[:, 2]
        base = 3
    else:  # SH12A: Z <10 quantities>
        z = d[:, 0]
        base = 1
    if d.shape[1] < base + len(QUANTITIES):
        return None
    cols = {q: d[:, base + i] for i, q in enumerate(QUANTITIES)}
    cols["z"] = z
    return cols


def osh_spectrum(out_dir: Path):
    """OpenShieldHIT plateau spectrum -> dPhi/dEkin (dividing Phi(bin) by bin width, issue #215)."""
    path = out_dir / "spectrum.dat"
    if not path.exists():
        return None
    d = load_numeric(path)  # X Y Z EKIN Phi_all Phi_prot
    if d.shape[1] < 6:
        return None
    ekin = d[:, -3]
    _, widths = log_bin_centers_widths(SPEC_LO, SPEC_HI, SPEC_NBINS)
    n = min(len(ekin), len(widths))
    return {
        "ekin": ekin[:n],
        "dphi_all": d[:n, -2] / widths[:n],
        "dphi_prot": d[:n, -1] / widths[:n],
    }


def sh12a_spectrum(path: Path):
    """SHIELD-HIT12A committed plateau spectrum (Page Diff1bin Value; already /cm^2/MeV)."""
    if not path.exists():
        return None
    d = load_numeric(path)  # page bin value
    if d.shape[1] < 3:
        return None
    centers, _ = log_bin_centers_widths(SPEC_LO, SPEC_HI, SPEC_NBINS)
    out = {"ekin": centers}
    for page, key in ((0, "dphi_all"), (1, "dphi_prot")):
        rows = d[d[:, 0].astype(int) == page]
        rows = rows[np.argsort(rows[:, 1])]  # order by bin index
        vals = rows[:, 2]
        n = min(len(centers), len(vals))
        out[key] = vals[:n]
    return out


def main() -> int:
    args = parse_args()
    if args.workdir is None:
        args.workdir = Path(tempfile.mkdtemp(prefix="nucre_"))
    args.workdir.mkdir(parents=True, exist_ok=True)
    if args.out is None:
        args.out = args.workdir / "nucre_report.pdf"

    try:
        version = subprocess.check_output(
            ["git", "-C", str(args.repo_root), "describe", "--tags", "--dirty", "--always"],
            stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        version = "unknown"

    # --- Run OpenShieldHIT for every case concurrently (SH12A is a fixed fixture). ---
    osh_res = {}
    if args.no_run:
        for case in args.cases:
            osh_res[case] = (args.workdir / "osh_out" / case, "", None)
    else:
        tasks = {}
        with concurrent.futures.ThreadPoolExecutor(max_workers=len(args.cases)) as ex:
            for case in args.cases:
                print(f"[{case}] launching OpenShieldHIT (nstat={args.nstat}) ...", flush=True)
                tasks[ex.submit(run_osh, args, case)] = case
            for fut in concurrent.futures.as_completed(tasks):
                case = tasks[fut]
                out, err, elapsed = fut.result()
                osh_res[case] = (out, err, elapsed)
                if err:
                    print(f"[{case}] OpenShieldHIT unavailable after {format_elapsed(elapsed)}: {err}", flush=True)
                else:
                    print(f"[{case}] OpenShieldHIT done in {format_elapsed(elapsed)}", flush=True)

    data = {}
    for case in args.cases:
        osh_out, osh_err, osh_s = osh_res[case]
        osh = depth_quantities(osh_out, "osh") if osh_out is not None else None
        osh_spec = osh_spectrum(osh_out) if osh_out is not None else None
        sh = sh_spec = None
        if case in SH12A_CASES:
            ref = sh12a_ref_dir(args, case)
            sh = depth_quantities(ref, "sh12a")
            sh_spec = sh12a_spectrum(ref / "spectrum.dat")
        data[case] = {"osh": osh, "osh_err": osh_err, "osh_s": osh_s, "osh_spec": osh_spec,
                      "sh": sh, "sh_spec": sh_spec}

    DOSE_SPECIES = [
        ("dose", "all", "C0"), ("dose_prim", "primary p", "C1"),
        ("dose_prot", "all p", "C2"), ("dose_alpha", "alphas", "C3"),
        ("dose_heavy", "heavy rec (Z>=3)", "C4"),
    ]
    FLU_SPECIES = [
        ("flu", "all", "C0"), ("flu_prim", "primary p", "C1"),
        ("flu_prot", "all p", "C2"), ("flu_alpha", "alphas", "C3"),
        ("flu_heavy", "heavy rec (Z>=3)", "C4"),
    ]

    with PdfPages(args.out) as pdf:
        for case in args.cases:
            d = data[case]
            osh, sh = d["osh"], d["sh"]
            if osh is None and sh is None:
                continue
            ref = osh or sh
            zend = max((c["z"].max() for c in (osh, sh) if c), default=ref["z"].max())

            fig, ax = plt.subplots(2, 2, figsize=(11, 8))
            note = ""
            if osh is None:
                note = f"\n[OpenShieldHIT unavailable: {d['osh_err']}]"
            elif case not in SH12A_CASES:
                note = "\n[OpenShieldHIT-only NUCRE mode]"
            elif sh is None:
                note = "\n[SHIELD-HIT12A fixture missing]"
            fig.suptitle(f"200 MeV p -> water, {MODEL[case]}   (MSCAT off, STRAGG off; "
                         f"OSH nstat={args.nstat})\nOpenShieldHIT {version} (solid) vs "
                         f"SHIELD-HIT12A committed fixture (dashed){note}", fontsize=11)

            def overlay_species(a, families):
                for c, style, tag in ((osh, "-", "OSH"), (sh, "--", "SH12A")):
                    if c is None:
                        continue
                    for key, label, color in families:
                        a.plot(c["z"], c[key], style, color=color,
                               label=f"{tag} {label}", drawstyle="steps-mid")
                a.set_xlim(0.0, zend)

            a = ax[0, 0]
            overlay_species(a, DOSE_SPECIES)
            a.set(title="Depth dose by species, 1 cm^2 column", xlabel="Depth (cm)", ylabel="Dose per primary (MeV/g)")
            a.grid(True)
            a.legend(fontsize=7, ncol=2)

            a = ax[0, 1]
            for c, style, tag in ((osh, "-", "OSH"), (sh, "--", "SH12A")):
                if c is None:
                    continue
                a.plot(c["z"], c["dose"] - c["dose_prim"], style, color="C0",
                       label=f"{tag} all-primary", drawstyle="steps-mid")
                a.plot(c["z"], c["dose_alpha"], style, color="C3", label=f"{tag} alphas", drawstyle="steps-mid")
                a.plot(c["z"], c["dose_heavy"], style, color="C4", label=f"{tag} heavy rec", drawstyle="steps-mid")
            a.set(title="Secondary dose (nuclear-added)", xlabel="Depth (cm)", ylabel="Dose per primary (MeV/g)")
            a.set_xlim(0.0, zend)
            a.grid(True)
            a.legend(fontsize=7)

            a = ax[1, 0]
            overlay_species(a, FLU_SPECIES)
            a.set(title="Depth fluence by species, 1 cm^2 column", xlabel="Depth (cm)", ylabel="Fluence per primary (1/cm^2)")
            a.grid(True)
            a.legend(fontsize=7, ncol=2)

            a = ax[1, 1]
            for spec, style, tag in ((d["osh_spec"], "-", "OSH"), (d["sh_spec"], "--", "SH12A")):
                if spec is None:
                    continue
                a.plot(spec["ekin"], spec["dphi_all"], style, color="C0", label=f"{tag} all")
                a.plot(spec["ekin"], spec["dphi_prot"], style, color="C2", label=f"{tag} protons")
            a.set(title="Plateau secondary spectrum (z=9.5-10.5 cm)",
                  xlabel="Ekin (MeV)", ylabel="dPhi/dEkin per primary (1/cm^2/MeV)")
            a.set_xscale("log")
            a.set_yscale("log")
            a.grid(True, which="both", alpha=0.3)
            handles, _ = a.get_legend_handles_labels()
            if handles:
                a.legend(fontsize=7)

            fig.tight_layout(rect=(0, 0, 1, 0.92))
            pdf.savefig(fig)
            plt.close(fig)

        # Summary: all-particle dose and plateau proton spectrum across NUCRE modes.
        fig, ax = plt.subplots(1, 2, figsize=(13, 6))
        for i, case in enumerate(args.cases):
            for c, style, alpha, tag in ((data[case]["osh"], "-", 1.0, "OSH"),
                                         (data[case]["sh"], "--", 0.7, "SH12A")):
                if c is None:
                    continue
                ax[0].plot(c["z"], c["dose"], style, color=f"C{i}", alpha=alpha,
                           label=f"{tag} {MODEL[case]}", drawstyle="steps-mid")
            for spec, style, alpha, tag in ((data[case]["osh_spec"], "-", 1.0, "OSH"),
                                            (data[case]["sh_spec"], "--", 0.7, "SH12A")):
                if spec is None:
                    continue
                ax[1].plot(spec["ekin"], spec["dphi_prot"], style, color=f"C{i}", alpha=alpha,
                           label=f"{tag} {MODEL[case]}")
        ax[0].set(title="All-particle depth dose vs NUCRE mode", xlabel="Depth (cm)", ylabel="Dose per primary (MeV/g)")
        ax[0].grid(True)
        ax[1].set(title="Plateau proton spectrum vs NUCRE mode", xlabel="Ekin (MeV)",
                  ylabel="dPhi/dEkin per primary (1/cm^2/MeV)")
        ax[1].set_xscale("log")
        ax[1].set_yscale("log")
        ax[1].grid(True, which="both", alpha=0.3)
        for a in ax:
            handles, _ = a.get_legend_handles_labels()
            if handles:
                a.legend(fontsize=8)
        fig.suptitle(f"OpenShieldHIT {version} (solid) vs SHIELD-HIT12A committed fixtures (dashed)", fontsize=11)
        fig.tight_layout(rect=(0, 0, 1, 0.95))
        pdf.savefig(fig)
        plt.close(fig)

    print(f"\nWrote {args.out}")
    print(f"\nOpenShieldHIT wall-clock (nstat={args.nstat}); SHIELD-HIT12A = committed fixture:")
    for case in args.cases:
        print(f"  {case:<8}{format_elapsed(data[case]['osh_s'])}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
