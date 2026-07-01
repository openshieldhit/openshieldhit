#!/usr/bin/env python3
"""Overlay OpenSHIELDHIT vs SHIELD-HIT12A for the issue #190 straggling cases.

The cases tests/reference/idd_water_200mev_strag{0,1,2} (200 MeV protons in water,
MSCAT off, NUCRE off; STRAGG 0/1/2) isolate the energy-straggling distal-edge
behaviour.  This tool runs both codes at low NSTAT and writes a multi-page PDF
comparing, per case, on the single 1 cm^2 central column (idd.dat):

  * depth dose per primary (Bragg curve) -- the auto-compared observable;
  * primary-proton fluence vs depth;
  * dose-averaged LET (DLET) vs depth;
  * track-averaged LET (TLET) vs depth.

A final summary page overlays the dose distal edge and the DLET tail of all three
STRAGG modes for both codes.

NOTE: STRAGG 2 (Vavilov) is not yet implemented in OpenSHIELDHIT (issue #190);
that run is reported as unavailable until it lands.

Column layouts handled (idd.dat = Dose, Fluence, DLET, TLET in detect.dat order):
  OpenSHIELDHIT mesh : X Y Z DOSE FLUENCE DLET TLET   (Z=col2, quantities col3..6)
  SHIELD-HIT12A 1D   : Z DOSE FLUENCE DLET TLET        (Z=col0, quantities col1..4)
"""

from __future__ import annotations

import argparse
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

CASES = ["strag0", "strag1", "strag2"]
MODEL = {
    "strag0": "STRAGG 0 (off)",
    "strag1": "STRAGG 1 (Gaussian)",
    "strag2": "STRAGG 2 (Vavilov)",
}


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parent.parent
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--nstat", type=int, default=20000, help="primaries per run (default 20000)")
    p.add_argument("--repo-root", type=Path, default=root, help="openshieldhit repo root")
    p.add_argument("--osh", type=Path, default=root / "build" / "bin" / "openshieldhit", help="OpenSHIELDHIT binary")
    p.add_argument("--sh12a", default="shieldhit", help="SHIELD-HIT12A binary (in PATH or absolute)")
    p.add_argument("--workdir", type=Path, default=None, help="working dir for runs (default: a temp dir)")
    p.add_argument("--out", type=Path, default=None, help="output PDF (default: <workdir>/straggling_report.pdf)")
    p.add_argument("--cases", nargs="+", default=CASES, choices=CASES, help="cases to include")
    p.add_argument("--no-run", action="store_true", help="reuse existing outputs in --workdir, do not run the codes")
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
                continue  # stray axis-label line (e.g. SH12A "Z [cm]")
    if not rows:
        raise ValueError(f"{path}: no numeric data")
    return np.array(rows)


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


def format_speedup(osh_s, sh_s):
    if osh_s is None or sh_s is None or osh_s <= 0.0:
        return "n/a"
    return f"{sh_s / osh_s:.2f}x"


def _run(cmd, cwd=None):
    """Run a command; return (ok, last_stderr_line, elapsed_s)."""
    start = time.perf_counter()
    r = subprocess.run(cmd, cwd=cwd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    elapsed_s = time.perf_counter() - start
    if r.returncode == 0:
        return True, "", elapsed_s
    errs = [ln.strip() for ln in r.stderr.splitlines() if ln.strip()]
    # Prefer the most informative line (e.g. "Vavilov ... not implemented").
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


def run_sh12a(args, case: str):
    src = args.repo_root / "tests" / "reference" / "shieldhit" / f"idd_water_200mev_{case}"
    if not src.is_dir():
        return None, "no SH12A deck in tests/reference/shieldhit/", None
    deck = args.workdir / "sh12a_deck" / case
    set_nstat(src, deck, args.nstat)
    ok, err, elapsed_s = _run([args.sh12a, "."], cwd=deck)
    return (deck, "", elapsed_s) if ok else (None, err, elapsed_s)


# idd.dat quantity order (must match detect.dat): Dose, Fluence, DLET, TLET.
QUANTITIES = ("dose", "fluence", "dlet", "tlet")


def depth_quantities(out_dir: Path, code: str):
    """Return dict {z, dose, fluence, dlet, tlet} for the 1D idd.dat depth file."""
    d = load_numeric(out_dir / "idd.dat")
    if code == "osh":  # X Y Z DOSE FLUENCE DLET TLET
        z = d[:, 2]
        cols = {q: d[:, 3 + i] for i, q in enumerate(QUANTITIES)}
    else:  # SH12A: Z DOSE FLUENCE DLET TLET
        z = d[:, 0]
        cols = {q: d[:, 1 + i] for i, q in enumerate(QUANTITIES)}
    cols["z"] = z
    return cols


def main() -> int:
    args = parse_args()
    if args.workdir is None:
        args.workdir = Path(tempfile.mkdtemp(prefix="straggling_"))
    args.workdir.mkdir(parents=True, exist_ok=True)
    if args.out is None:
        args.out = args.workdir / "straggling_report.pdf"

    try:
        version = subprocess.check_output(
            ["git", "-C", str(args.repo_root), "describe", "--tags", "--dirty", "--always"],
            stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        version = "unknown"

    def load_all(out_dir, code):
        if out_dir is None or not (out_dir / "idd.dat").exists():
            return None
        return depth_quantities(out_dir, code)

    data = {}
    for case in args.cases:
        if args.no_run:
            osh_out, osh_err, osh_s = args.workdir / "osh_out" / case, "", None
            sh_out, sh_err, sh_s = args.workdir / "sh12a_deck" / case, "", None
        else:
            print(f"[{case}] running OpenSHIELDHIT (nstat={args.nstat}) ...", flush=True)
            osh_out, osh_err, osh_s = run_osh(args, case)
            if osh_err:
                print(f"[{case}] OpenSHIELDHIT unavailable after {format_elapsed(osh_s)}: {osh_err}", flush=True)
            else:
                print(f"[{case}] OpenSHIELDHIT done in {format_elapsed(osh_s)}", flush=True)
            print(f"[{case}] running SHIELD-HIT12A (nstat={args.nstat}) ...", flush=True)
            sh_out, sh_err, sh_s = run_sh12a(args, case)
            if sh_err:
                print(f"[{case}] SHIELD-HIT12A unavailable after {format_elapsed(sh_s)}: {sh_err}", flush=True)
            else:
                print(f"[{case}] SHIELD-HIT12A done in {format_elapsed(sh_s)}", flush=True)
        data[case] = {"osh": load_all(osh_out, "osh"), "osh_err": osh_err,
                      "osh_s": osh_s, "sh": load_all(sh_out, "sh12a"),
                      "sh_err": sh_err, "sh_s": sh_s}

    with PdfPages(args.out) as pdf:
        for case in args.cases:
            d = data[case]
            osh, sh = d["osh"], d["sh"]
            if osh is None and sh is None:
                continue
            ref = osh or sh  # at least one present; use it to pick the peak depth
            zpeak = ref["z"][np.argmax(ref["dose"])]
            zend = max((c["z"].max() for c in (osh, sh) if c), default=ref["z"].max())

            fig, ax = plt.subplots(2, 2, figsize=(11, 8))
            note = ""
            if osh is None:
                note = f"\n[OpenSHIELDHIT unavailable: {d['osh_err']}]"
            elif sh is None:
                note = f"\n[SHIELD-HIT12A unavailable: {d['sh_err']}]"
            timing = (f"\nTiming: OpenSHIELDHIT {format_elapsed(d['osh_s'])}, "
                      f"SHIELD-HIT12A {format_elapsed(d['sh_s'])}")
            fig.suptitle(f"200 MeV p -> water, {MODEL[case]}   (MSCAT off, NUCRE off, "
                         f"nstat={args.nstat})\nOpenSHIELDHIT {version} vs SHIELD-HIT12A{note}{timing}",
                         fontsize=11)

            def overlay(a, quantity):
                for c, style, name in ((osh, "-", "OpenSHIELDHIT"), (sh, "--", "SHIELD-HIT12A")):
                    if c is None:
                        continue
                    color = "C0" if name.startswith("Open") else "C1"
                    a.plot(c["z"], c[quantity], style, color=color, label=name, drawstyle="steps-mid")

            a = ax[0, 0]
            overlay(a, "dose")
            a.set(title="Depth dose, 1 cm^2 column", xlabel="Depth (cm)", ylabel="Dose per primary")
            a.set_xlim(zpeak - 4, zend)
            a.grid(True)
            a.legend(loc="upper left")

            a = ax[0, 1]
            overlay(a, "fluence")
            a.set(title="Primary fluence, 1 cm^2 column", xlabel="Depth (cm)", ylabel="Fluence (1/cm^2)")
            a.set_xlim(zpeak - 4, zend)
            a.grid(True)
            a.legend(loc="lower left")

            a = ax[1, 0]
            overlay(a, "dlet")
            a.set(title="Dose-averaged LET (DLET)", xlabel="Depth (cm)", ylabel="DLET (MeV/cm)")
            a.set_xlim(zpeak - 4, zend)
            a.grid(True)
            a.legend(loc="upper left")

            a = ax[1, 1]
            overlay(a, "tlet")
            a.set(title="Track-averaged LET (TLET)", xlabel="Depth (cm)", ylabel="TLET (MeV/cm)")
            a.set_xlim(zpeak - 4, zend)
            a.grid(True)
            a.legend(loc="upper left")

            fig.tight_layout(rect=(0, 0, 1, 0.93))
            pdf.savefig(fig)
            plt.close(fig)

        # Summary: dose distal edge and DLET tail across STRAGG modes.
        fig, ax = plt.subplots(1, 2, figsize=(13, 6))
        zends = []
        for i, case in enumerate(args.cases):
            for c, style, alpha, tag in ((data[case]["osh"], "-", 1.0, "OSH"),
                                         (data[case]["sh"], "--", 0.7, "SH12A")):
                if c is None:
                    continue
                zends.append(c["z"].max())
                ax[0].plot(c["z"], c["dose"], style, color=f"C{i}", alpha=alpha,
                           label=f"{tag} {MODEL[case]}", drawstyle="steps-mid")
                ax[1].plot(c["z"], c["dlet"], style, color=f"C{i}", alpha=alpha,
                           label=f"{tag} {MODEL[case]}", drawstyle="steps-mid")
        if zends:
            for a in ax:
                a.set_xlim(max(zends) - 5, max(zends))
        ax[0].set(title="Dose distal edge vs STRAGG mode", xlabel="Depth (cm)", ylabel="Dose per primary")
        ax[1].set(title="DLET tail vs STRAGG mode", xlabel="Depth (cm)", ylabel="DLET (MeV/cm)")
        for a in ax:
            a.grid(True)
            handles, _ = a.get_legend_handles_labels()
            if handles:
                a.legend(fontsize=8)
        fig.suptitle(f"OpenSHIELDHIT {version} (solid) vs SHIELD-HIT12A (dashed)", fontsize=11)
        timing_lines = ["Wall-clock timings:"]
        for case in args.cases:
            d = data[case]
            timing_lines.append(f"{MODEL[case]}: OSH {format_elapsed(d['osh_s'])}, "
                                f"SH12A {format_elapsed(d['sh_s'])}")
        fig.text(0.01, 0.01, "  |  ".join(timing_lines), fontsize=8)
        fig.tight_layout(rect=(0, 0.05, 1, 0.95))
        pdf.savefig(fig)
        plt.close(fig)

    print(f"\nWrote {args.out}")
    print(f"\nTiming summary (wall clock, nstat={args.nstat}; speedup = SH12A / OSH):")
    print(f"{'case':<8}{'OpenSHIELDHIT':>16}{'SHIELD-HIT12A':>16}{'speedup':>12}")
    for case in args.cases:
        d = data[case]
        print(f"{case:<8}{format_elapsed(d['osh_s']):>16}{format_elapsed(d['sh_s']):>16}"
              f"{format_speedup(d['osh_s'], d['sh_s']):>12}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
