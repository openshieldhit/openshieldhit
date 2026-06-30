#!/usr/bin/env python3
"""Overlay OpenSHIELDHIT vs SHIELD-HIT12A for the issue #133 MCS-isolation cases.

The cases tests/reference/idd_water_200mev_scat{0,1,2} (200 MeV protons in water,
NUCRE off, STRAGG off; MSCAT 0/1/2) isolate the multiple-scattering distal-edge
behaviour.  This tool runs both codes at low NSTAT (no nuclear reactions and no
straggling -> the curves converge quickly) and writes a multi-page PDF comparing,
per case:

  * primary-proton fluence vs depth on the 1 cm^2 central column (idd.dat) and
    full-width laterally integrated (ddc_wide.dat) -- the issue #133 observable;
  * energy deposition per primary vs depth (Bragg curve);
  * radial primary-fluence profile at the Bragg-peak depth (rad_cyl.dat).

A final summary page overlays the narrow-column fluence distal edges of all three
MSCAT modes for both codes.

Column layouts handled:
  OpenSHIELDHIT mesh : X Y Z FLUENCE ENERGY          (Z=col2, fluence=col3)
  OpenSHIELDHIT cyl  : Z R FLUENCE ENERGY            (Z=col0, R=col1)
  SHIELD-HIT12A 1D   : Z FLUENCE ENERGY              (Z=col0, fluence=col1)
  SHIELD-HIT12A cyl  : R Z FLUENCE ENERGY            (R=col0, Z=col1)
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

CASES = ["scat0", "scat1", "scat2"]
MODEL = {"scat0": "MSCAT 0 (no scatter)", "scat1": "MSCAT 1 (Gaussian)", "scat2": "MSCAT 2 (Moliere)"}


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parent.parent
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--nstat", type=int, default=8000, help="primaries per run (default 8000)")
    p.add_argument("--repo-root", type=Path, default=root, help="openshieldhit repo root")
    p.add_argument("--osh", type=Path, default=root / "build" / "bin" / "openshieldhit", help="OpenSHIELDHIT binary")
    p.add_argument("--sh12a", default="shieldhit", help="SHIELD-HIT12A binary (in PATH or absolute)")
    p.add_argument("--workdir", type=Path, default=None, help="working dir for runs (default: a temp dir)")
    p.add_argument("--out", type=Path, default=None, help="output PDF (default: <workdir>/mcs_scat_report.pdf)")
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
    # Prefer the most informative line over the generic failure tail.
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
    deck = args.workdir / "sh12a_deck" / case
    set_nstat(src, deck, args.nstat)
    ok, err, elapsed_s = _run([args.sh12a, "."], cwd=deck)
    return (deck, "", elapsed_s) if ok else (None, err, elapsed_s)


def depth_fluence_energy(out_dir: Path, fname: str, code: str):
    """Return (z, fluence, energy) for a 1D depth file from either code."""
    d = load_numeric(out_dir / fname)
    if code == "osh":  # X Y Z FLUENCE ENERGY
        return d[:, 2], d[:, 3], d[:, 4]
    return d[:, 0], d[:, 1], d[:, 2]  # SH12A: Z FLUENCE ENERGY


def radial_at_depth(out_dir: Path, code: str, z_target: float):
    """Return (r, fluence, z_slab) for the last depth slab at/before z_target.

    The radial bins are coarse (2 cm); the slab whose centre is just *past* the
    Bragg peak is essentially empty (beam stopped), so pick the deepest slab
    centred at or before z_target, where there is still primary fluence to show.
    """
    d = load_numeric(out_dir / "rad_cyl.dat")
    if code == "osh":  # Z R FLUENCE ENERGY
        z, r, fl = d[:, 0], d[:, 1], d[:, 2]
    else:  # SH12A: R Z FLUENCE ENERGY
        r, z, fl = d[:, 0], d[:, 1], d[:, 2]
    centers = np.unique(z)
    at_or_before = centers[centers <= z_target + 1e-9]
    zsel = at_or_before.max() if at_or_before.size else centers[np.argmin(np.abs(centers - z_target))]
    m = np.isclose(z, zsel)
    order = np.argsort(r[m])
    return r[m][order], fl[m][order], zsel


def main() -> int:
    args = parse_args()
    if args.workdir is None:
        args.workdir = Path(tempfile.mkdtemp(prefix="mcs_scat_"))
    args.workdir.mkdir(parents=True, exist_ok=True)
    if args.out is None:
        args.out = args.workdir / "mcs_scat_report.pdf"

    try:
        version = subprocess.check_output(
            ["git", "-C", str(args.repo_root), "describe", "--tags", "--dirty", "--always"],
            stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        version = "unknown"

    def load_all(out_dir, code):
        if out_dir is None or not (out_dir / "idd.dat").exists():
            return None
        return {
            "narrow": depth_fluence_energy(out_dir, "idd.dat", code),
            "wide": depth_fluence_energy(out_dir, "ddc_wide.dat", code),
            "dir": out_dir,
        }

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
            zr_, _, er_ = ref["narrow"]
            zpeak = zr_[np.argmax(er_)]
            zend = max((c["narrow"][0].max() for c in (osh, sh) if c), default=zr_.max())

            fig, ax = plt.subplots(2, 2, figsize=(11, 8))
            note = ""
            if osh is None:
                note = f"\n[OpenSHIELDHIT unavailable: {d['osh_err']}]"
            elif sh is None:
                note = f"\n[SHIELD-HIT12A unavailable: {d['sh_err']}]"
            timing = (f"\nTiming: OpenSHIELDHIT {format_elapsed(d['osh_s'])}, "
                      f"SHIELD-HIT12A {format_elapsed(d['sh_s'])}")
            fig.suptitle(f"200 MeV p -> water, {MODEL[case]}   (NUCRE off, STRAGG off, "
                         f"nstat={args.nstat})\nOpenSHIELDHIT {version} vs SHIELD-HIT12A{note}{timing}",
                         fontsize=11)

            def overlay(a, key, idx):
                for c, style, name in ((osh, "-", "OpenSHIELDHIT"), (sh, "--", "SHIELD-HIT12A")):
                    if c is None:
                        continue
                    z = c[key][0]
                    y = c[key][idx]
                    color = "C0" if name.startswith("Open") else "C1"
                    a.plot(z, y, style, color=color, label=name, drawstyle="steps-mid")

            a = ax[0, 0]
            overlay(a, "narrow", 1)
            a.set(title="Primary fluence, 1 cm^2 column", xlabel="Depth (cm)", ylabel="Fluence (1/cm^2)")
            a.set_xlim(zpeak - 4, zend)
            a.grid(True)
            a.legend(loc="lower left")

            a = ax[0, 1]
            overlay(a, "wide", 1)
            a.set(title="Primary fluence, full-width (lat. integrated)", xlabel="Depth (cm)", ylabel="Fluence (1/cm^2)")
            a.set_xlim(zpeak - 4, zend)
            a.grid(True)
            a.legend(loc="lower left")

            a = ax[1, 0]
            overlay(a, "narrow", 2)
            a.set(title="Energy deposition, 1 cm^2 column", xlabel="Depth (cm)", ylabel="Energy deposition per primary")
            a.grid(True)
            a.legend(loc="upper left")

            a = ax[1, 1]
            zsel_used = zpeak
            for c, mk, name in ((osh, "o-", "OpenSHIELDHIT"), (sh, "s--", "SHIELD-HIT12A")):
                if c is None:
                    continue
                code = "osh" if name.startswith("Open") else "sh12a"
                r, fr, zsel_used = radial_at_depth(c["dir"], code, zpeak)
                color = "C0" if code == "osh" else "C1"
                a.plot(r, fr, mk, color=color, ms=3, label=name, drawstyle="steps-mid")
            a.set(title=f"Radial primary fluence at z = {zsel_used:.0f} cm (near Bragg peak)",
                  xlabel="Radius (cm)", ylabel="Fluence (1/cm^2)")
            a.grid(True)
            a.legend(loc="upper right")

            fig.tight_layout(rect=(0, 0, 1, 0.93))
            pdf.savefig(fig)
            plt.close(fig)

        # Summary: distal edges of the narrow fluence across MSCAT modes
        fig, ax = plt.subplots(figsize=(9, 6))
        zends = []
        for i, case in enumerate(args.cases):
            for c, style, alpha, tag in ((data[case]["osh"], "-", 1.0, "OSH"),
                                         (data[case]["sh"], "--", 0.7, "SH12A")):
                if c is None:
                    continue
                z, fl, _ = c["narrow"]
                zends.append(z.max())
                ax.plot(z, fl, style, color=f"C{i}", alpha=alpha, label=f"{tag} {MODEL[case]}",
                        drawstyle="steps-mid")
        if zends:
            ax.set_xlim(max(zends) - 5, max(zends))
        ax.set(title="Primary-fluence distal edge (1 cm^2 column) vs MSCAT mode\n"
               f"OpenSHIELDHIT {version} (solid) vs SHIELD-HIT12A (dashed)",
               xlabel="Depth (cm)", ylabel="Fluence (1/cm^2)")
        ax.grid(True)
        handles, _ = ax.get_legend_handles_labels()
        if handles:
            ax.legend(fontsize=8)
        timing_lines = ["Wall-clock timings:"]
        for case in args.cases:
            d = data[case]
            timing_lines.append(f"{MODEL[case]}: OSH {format_elapsed(d['osh_s'])}, "
                                f"SH12A {format_elapsed(d['sh_s'])}")
        fig.text(0.01, 0.01, "\n".join(timing_lines), fontsize=8)
        fig.tight_layout(rect=(0, 0.11, 1, 1))
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
