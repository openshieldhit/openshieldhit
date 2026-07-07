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
    Dose Deuterons, Dose Tritons, Dose He3, Dose Alphas, Dose HeavyRec,
    Fluence, Fluence Primary, Fluence Protons, Fluence Deuterons,
    Fluence Tritons, Fluence He3, Fluence Alphas, Fluence HeavyRec
    OpenShieldHIT mesh : X Y Z <16 quantities>   (Z=col2, quantities col3..18)
    SHIELD-HIT12A 1D   : Z <16 quantities>        (Z=col0, quantities col1..16)
  spectrum.dat (Ekin 0.1..300 MeV, 150 log bins, pages: Fluence, Fluence Protons,
    Fluence Deuterons, Fluence Tritons, Fluence He3, Fluence Alphas)
    OpenShieldHIT : X Y Z EKIN Phi_all Phi_p Phi_d Phi_t Phi_He3 Phi_alpha
      -- Phi PER BIN (see issue #215)
    SHIELD-HIT12A : Page Diff1bin Value
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
from typing import Any

os.environ.setdefault("MPLCONFIGDIR", "/tmp/openshieldhit-matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

# Reading aid: draw the SHIELD-HIT12A reference as a THICK, LIGHTENED underlay and
# OpenShieldHIT as a THIN, saturated line on top of it, both as steps (bins).  So
# the eye reads species by colour and code by weight, with OSH sitting "inside"
# the fat reference band.
REF_LW, OSH_LW = 3.0, 1.3
REF_LIGHTEN = 0.6  # 0 = original colour, 1 = white


def lighten(color: Any, amount: float = REF_LIGHTEN) -> tuple[float, float, float]:
    r, g, b = mcolors.to_rgb(color)
    return (r + (1.0 - r) * amount, g + (1.0 - g) * amount, b + (1.0 - b) * amount)


def code_style(code: str, color: Any) -> dict[str, Any]:
    """Line kwargs for a series: 'sh' = thick lightened underlay, 'osh' = thin saturated overlay."""
    if code == "sh":
        return dict(color=lighten(color), lw=REF_LW, zorder=1, drawstyle="steps-mid")
    return dict(color=color, lw=OSH_LW, zorder=3, drawstyle="steps-mid")

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
    "dose", "dose_prim", "dose_prot", "dose_deut", "dose_trit", "dose_he3", "dose_alpha", "dose_heavy",
    "flu", "flu_prim", "flu_prot", "flu_deut", "flu_trit", "flu_he3", "flu_alpha", "flu_heavy",
)

# let.dat quantity order (must match detect.dat).
LET_QUANTITIES = (
    "dlet", "dlet_prim", "dlet_prot", "dlet_deut", "dlet_trit", "dlet_he3", "dlet_alpha", "dlet_heavy",
    "tlet", "tlet_prim", "tlet_prot", "tlet_deut", "tlet_trit", "tlet_he3", "tlet_alpha", "tlet_heavy",
)

# Differential spectrum axis (must match the Diff1 line in detect.dat).
SPEC_LO, SPEC_HI, SPEC_NBINS = 0.1, 300.0, 150

# Differential spectrum page/column order (must match detect.dat).
SPECTRUM_SPECIES = (
    ("dphi_all", "all", "C0"),
    ("dphi_prot", "p", "C2"),
    ("dphi_deut", "d", "C3"),
    ("dphi_trit", "t", "C4"),
    ("dphi_he3", "He3", "C5"),
    ("dphi_alpha", "alpha", "C6"),
)


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parent.parent
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-n", "--nstat", type=int, default=50000, help="OpenShieldHIT primaries per run (default 50000)")
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


def format_nstat(nstat):
    """Human-readable primary count label for plot titles."""
    if nstat is None:
        return "n/a"
    if nstat >= 1.0e6:
        return f"{nstat / 1.0e6:g}M"
    if nstat >= 1.0e3:
        return f"{nstat / 1.0e3:g}k"
    return f"{nstat:g}"


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


def read_nstat(beam_path: Path):
    """Return the NSTAT (number of primaries) from a beam.dat, or None."""
    if not beam_path.exists():
        return None
    for line in beam_path.read_text().splitlines():
        toks = line.split()
        if toks and toks[0].upper() == "NSTAT":
            try:
                return float(toks[1])
            except (IndexError, ValueError):
                return None
    return None


def load_depth(out_dir: Path, code: str, filename: str, quantities):
    """Return dict {z, <quantities>} for a 1D depth file, or None if absent.

    OpenShieldHIT mesh output is `X Y Z <quantities>` (Z=col2); SHIELD-HIT12A 1D
    is `Z <quantities>` (Z=col0).
    """
    path = out_dir / filename
    if not path.exists():
        return None
    d = load_numeric(path)
    base = 3 if code == "osh" else 1
    z = d[:, 2] if code == "osh" else d[:, 0]
    if d.shape[1] < base + len(quantities):
        return None
    cols = {q: d[:, base + i] for i, q in enumerate(quantities)}
    cols["z"] = z
    return cols


def depth_quantities(out_dir: Path, code: str):
    """Species-resolved dose/fluence depth curves from idd.dat."""
    return load_depth(out_dir, code, "idd.dat", QUANTITIES)


def depth_let(out_dir: Path, code: str):
    """Species-resolved DLET/TLET depth curves from let.dat (may be absent)."""
    return load_depth(out_dir, code, "let.dat", LET_QUANTITIES)


def rel_diff(osh, sh, key):
    """(OSH - SH12A)/SH12A vs depth for one quantity, masked in the low-reference
    tail where the ratio explodes.  Returns (z, reldiff) or (None, None)."""
    if osh is None or sh is None or key not in osh or key not in sh:
        return None, None
    n = min(len(osh[key]), len(sh[key]), len(osh["z"]))
    o = np.asarray(osh[key][:n], dtype=float)
    s = np.asarray(sh[key][:n], dtype=float)
    z = np.asarray(osh["z"][:n], dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        rd = (o - s) / s
    smax = np.nanmax(np.abs(s)) if n else 0.0
    if smax > 0.0:
        rd = np.where(np.abs(s) > 0.01 * smax, rd, np.nan)  # hide the noisy tail
    return z, rd


def draw_residual_strip(ax, osh, sh, families, zend, ylim_pct=30.0):
    """Per-species (OSH-SH12A)/SH12A in % under a main panel, sharing its x-axis.

    One line per species, coloured to match the main panel above (whose legend
    labels them), so the strip is 'coupled' to the plot it sits under.
    """
    drew = False
    for key, _label, color in families:
        z, rd = rel_diff(osh, sh, key)
        if rd is None:
            continue
        ax.plot(z, 100.0 * rd, color=color, lw=1.0, drawstyle="steps-mid")
        drew = True
    if drew:
        ax.axhline(0.0, color="k", lw=0.5)
        ax.set_ylim(-ylim_pct, ylim_pct)
    else:
        ax.text(0.5, 0.5, "no SH12A reference", ha="center", va="center",
                transform=ax.transAxes, fontsize=6, color="gray")
    ax.set_xlim(0.0, zend)
    ax.set_ylabel("(O-S)/S [%]", fontsize=6)
    ax.tick_params(labelsize=6)
    ax.grid(True, alpha=0.3)


def osh_spectrum(out_dir: Path):
    """OpenShieldHIT plateau spectrum -> dPhi/dEkin (dividing Phi(bin) by bin width, issue #215).

    New decks write columns X Y Z EKIN Phi_all Phi_p Phi_d Phi_t Phi_He3 Phi_alpha.
    Older three-page outputs (all, p, alpha) are still read correctly while the
    reference fixtures are being regenerated.
    """
    path = out_dir / "spectrum.dat"
    keys = ("dphi_all", "dphi_prot", "dphi_deut", "dphi_trit", "dphi_he3", "dphi_alpha")
    if not path.exists():
        return None
    d = load_numeric(path)
    if d.shape[1] < 6:
        return None
    _, widths = log_bin_centers_widths(SPEC_LO, SPEC_HI, SPEC_NBINS)
    n = min(d.shape[0], len(widths))
    out = {"ekin": d[:n, 3]}
    if d.shape[1] == 7:
        keys = ("dphi_all", "dphi_prot", "dphi_alpha")
    for i, key in enumerate(keys):
        col = 4 + i
        if col < d.shape[1]:
            out[key] = d[:n, col] / widths[:n]
    return out


def sh12a_spectrum(ref_dir: Path):
    """SHIELD-HIT12A committed plateau spectrum -> per-primary dPhi/dEkin.

    Despite the `[/cm^2/MeV]` header label, SH12A's spectrum.dat "Value" is the
    fluence summed into each bin over ALL primaries — absolute, not per-primary
    and not per-MeV (verified: sum over bins == NSTAT * per-primary fluence).  So
    normalise by NSTAT (from the fixture's beam.dat) and the log-bin width to get
    a per-primary spectral density comparable to the OSH curve.
    """
    path = ref_dir / "spectrum.dat"
    if not path.exists():
        return None
    d = load_numeric(path)  # page bin value
    if d.shape[1] < 3:
        return None
    nstat = read_nstat(ref_dir / "beam.dat")
    if not nstat:
        return None
    centers, widths = log_bin_centers_widths(SPEC_LO, SPEC_HI, SPEC_NBINS)
    out = {"ekin": centers}
    pages = sorted(set(d[:, 0].astype(int)))
    keys = ("dphi_all", "dphi_prot", "dphi_deut", "dphi_trit", "dphi_he3", "dphi_alpha")
    if pages == [0, 1, 2]:
        keys = ("dphi_all", "dphi_prot", "dphi_alpha")
    for page, key in enumerate(keys):
        rows = d[d[:, 0].astype(int) == page]
        if rows.size == 0:
            continue
        rows = rows[np.argsort(rows[:, 1])]  # order by bin index
        vals = rows[:, 2]
        n = min(len(centers), len(vals))
        out[key] = vals[:n] / (nstat * widths[:n])
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
        osh_let = depth_let(osh_out, "osh") if osh_out is not None else None
        osh_spec = osh_spectrum(osh_out) if osh_out is not None else None
        sh = sh_let = sh_spec = None
        sh_nstat = None
        if case in SH12A_CASES:
            ref = sh12a_ref_dir(args, case)
            sh_nstat = read_nstat(ref / "beam.dat")
            sh = depth_quantities(ref, "sh12a")
            sh_let = depth_let(ref, "sh12a")
            sh_spec = sh12a_spectrum(ref)
        data[case] = {"osh": osh, "osh_let": osh_let, "osh_err": osh_err, "osh_s": osh_s, "osh_spec": osh_spec,
                      "sh": sh, "sh_let": sh_let, "sh_spec": sh_spec, "sh_nstat": sh_nstat}

    DOSE_SPECIES = [
        ("dose", "all", "C0"), ("dose_prim", "primary p", "C1"),
        ("dose_prot", "all p", "C2"), ("dose_deut", "d", "C3"),
        ("dose_trit", "t", "C4"), ("dose_he3", "He3", "C5"),
        ("dose_alpha", "alphas", "C6"), ("dose_heavy", "heavy rec (Z>=3)", "C7"),
    ]
    FLU_SPECIES = [
        ("flu", "all", "C0"), ("flu_prim", "primary p", "C1"),
        ("flu_prot", "all p", "C2"), ("flu_deut", "d", "C3"),
        ("flu_trit", "t", "C4"), ("flu_he3", "He3", "C5"),
        ("flu_alpha", "alphas", "C6"), ("flu_heavy", "heavy rec (Z>=3)", "C7"),
    ]
    # Only the track-based species read cleanly; alphas/heavy DLET are dominated
    # by rare high-LET transported recoils in low-statistics bins, and heavy-recoil
    # LET is ~0 until point deposits feed LET (score_point energy/dose only today).
    DLET_SPECIES = [
        ("dlet", "all", "C0"), ("dlet_prim", "primary p", "C1"),
        ("dlet_prot", "all p", "C2"), ("dlet_deut", "d", "C3"),
        ("dlet_trit", "t", "C4"), ("dlet_he3", "He3", "C5"),
        ("dlet_alpha", "alphas", "C6"), ("dlet_heavy", "heavy rec (Z>=3)", "C7"),
    ]

    with PdfPages(args.out) as pdf:
        for case in args.cases:
            d = data[case]
            osh, sh = d["osh"], d["sh"]
            if osh is None and sh is None:
                continue
            ref = osh or sh
            zend = max((c["z"].max() for c in (osh, sh) if c), default=ref["z"].max())

            osh_let, sh_let = d["osh_let"], d["sh_let"]
            note = ""
            if osh is None:
                note = f"  [OpenShieldHIT unavailable: {d['osh_err']}]"
            elif case not in SH12A_CASES:
                note = "  [OpenShieldHIT-only NUCRE mode]"
            elif sh is None:
                note = "  [SHIELD-HIT12A fixture missing]"
            elif d["sh_nstat"] is not None:
                note = f"  [SH12A nstat={format_nstat(d['sh_nstat'])}]"

            # Each depth quantity is a main panel with a thin (OSH-SH12A)/SH12A
            # residual strip below it (shared x-axis); the spectrum spans the
            # lower-right.  SH12A underneath (thick, light); OSH on top (thin).
            fig = plt.figure(figsize=(12, 10), constrained_layout=True)
            gs = fig.add_gridspec(4, 2, height_ratios=[3, 1, 3, 1])
            fig.suptitle(f"200 MeV p -> water, {MODEL[case]}   (MSCAT off, STRAGG off; OSH nstat={args.nstat})\n"
                         f"OpenShieldHIT {version} (thin, saturated) vs SHIELD-HIT12A committed fixture "
                         f"(thick, light){note}", fontsize=11)

            def overlay(a, sh_src, osh_src, families):
                for c, code, tag in ((sh_src, "sh", "SH12A"), (osh_src, "osh", "OSH")):
                    if c is None:
                        continue
                    for key, label, color in families:
                        if key in c:
                            a.plot(c["z"], c[key], label=f"{tag} {label}", **code_style(code, color))
                a.set_xlim(0.0, zend)

            # --- Depth dose by species (+ residual) ---
            a = fig.add_subplot(gs[0, 0])
            overlay(a, sh, osh, DOSE_SPECIES)
            a.set(title="Depth dose by species, 1 cm^2 column", ylabel="Dose/primary (MeV/g)")
            a.set_yscale("log")
            a.set_ylim(1e-3, None)  # log so the alpha/heavy-recoil secondary dose is visible
            a.grid(True, which="both", alpha=0.3)
            a.legend(fontsize=6, ncol=2)
            a.tick_params(labelbottom=False)
            # Residual strip: the all/primary/all-proton family, where a % residual
            # is meaningful; the tiny alpha/heavy secondaries are 2-4x apart (off
            # this scale) and are compared directly in the log-y main panel above.
            draw_residual_strip(fig.add_subplot(gs[1, 0], sharex=a), osh, sh, DOSE_SPECIES[:3], zend)

            # --- Depth fluence by species (+ residual) ---
            a = fig.add_subplot(gs[2, 0])
            overlay(a, sh, osh, FLU_SPECIES)
            a.set(title="Depth fluence by species, 1 cm^2 column", ylabel="Fluence/primary (1/cm^2)")
            a.set_yscale("log")
            a.set_ylim(0.5, None)
            a.grid(True, which="both", alpha=0.3)
            a.legend(fontsize=6, ncol=2)
            a.tick_params(labelbottom=False)
            afs = fig.add_subplot(gs[3, 0], sharex=a)
            draw_residual_strip(afs, osh, sh, FLU_SPECIES[:3], zend)
            afs.set_xlabel("Depth (cm)")

            # --- Dose-averaged LET by species (+ residual) ---
            a = fig.add_subplot(gs[0, 1])
            overlay(a, sh_let, osh_let, DLET_SPECIES)
            a.set(title="Dose-averaged LET by species (DLET)", ylabel="DLET (MeV/cm)")
            a.set_yscale("log")
            a.grid(True, which="both", alpha=0.3)
            if osh_let is not None:
                a.legend(fontsize=6, ncol=2)
            else:
                a.text(0.5, 0.5, "no let.dat (run the case with the updated detect.dat)",
                       ha="center", va="center", transform=a.transAxes, fontsize=8, color="gray")
            a.tick_params(labelbottom=False)
            # DLET-all omitted from the residual: OSH's heavy-recoil LET is not yet
            # scored (point deposits feed dose, not LET), so DLET-all is not
            # comparable; primary/all-proton DLET are.
            draw_residual_strip(fig.add_subplot(gs[1, 1], sharex=a), osh_let, sh_let, DLET_SPECIES[1:], zend)

            # --- Plateau secondary spectrum (spans lower-right) ---
            a = fig.add_subplot(gs[2:, 1])
            for spec, code, tag in ((d["sh_spec"], "sh", "SH12A"), (d["osh_spec"], "osh", "OSH")):
                if spec is None:
                    continue
                for key, label, color in SPECTRUM_SPECIES:
                    if key in spec:
                        a.plot(spec["ekin"], spec[key], label=f"{tag} {label}", **code_style(code, color))
            a.set(title="Plateau secondary spectrum (z=9.5-10.5 cm)",
                  xlabel="Ekin (MeV)", ylabel="dPhi/dEkin per primary (1/cm^2/MeV)")
            a.set_xscale("log")
            a.set_yscale("log")
            a.grid(True, which="both", alpha=0.3)
            handles, _ = a.get_legend_handles_labels()
            if handles:
                a.legend(fontsize=7)

            pdf.savefig(fig)
            plt.close(fig)

        # Summary: all-particle dose and plateau proton spectrum across NUCRE modes.
        fig, ax = plt.subplots(1, 2, figsize=(13, 6))
        for i, case in enumerate(args.cases):
            # SH12A first (thick light underlay), OSH second (thin saturated overlay).
            for c, code, tag in ((data[case]["sh"], "sh", "SH12A"), (data[case]["osh"], "osh", "OSH")):
                if c is None:
                    continue
                ax[0].plot(c["z"], c["dose"], label=f"{tag} {MODEL[case]}", **code_style(code, f"C{i}"))
            for spec, code, tag in ((data[case]["sh_spec"], "sh", "SH12A"), (data[case]["osh_spec"], "osh", "OSH")):
                if spec is None:
                    continue
                ax[1].plot(spec["ekin"], spec["dphi_prot"], label=f"{tag} {MODEL[case]}", **code_style(code, f"C{i}"))
        ax[0].set(title="All-particle depth dose vs NUCRE mode", xlabel="Depth (cm)", ylabel="Dose per primary (MeV/g)")
        ax[0].set_yscale("log")
        ax[0].set_ylim(1.0, None)
        ax[0].grid(True, which="both", alpha=0.3)
        ax[1].set(title="Plateau proton spectrum vs NUCRE mode", xlabel="Ekin (MeV)",
                  ylabel="dPhi/dEkin per primary (1/cm^2/MeV)")
        ax[1].set_xscale("log")
        ax[1].set_yscale("log")
        ax[1].grid(True, which="both", alpha=0.3)
        for a in ax:
            handles, _ = a.get_legend_handles_labels()
            if handles:
                a.legend(fontsize=8)
        fig.suptitle(f"OpenShieldHIT {version} (thin, saturated) vs "
                     f"SHIELD-HIT12A committed fixtures (thick, light)", fontsize=11)
        fig.tight_layout(rect=(0, 0, 1, 0.95))
        pdf.savefig(fig)
        plt.close(fig)

    print(f"\nWrote {args.out}")
    print(f"\nOpenShieldHIT wall-clock (nstat={args.nstat}); SHIELD-HIT12A = committed fixture:")
    for case in args.cases:
        sh_label = ""
        if data[case]["sh_nstat"] is not None:
            sh_label = f"; SH12A nstat={format_nstat(data[case]['sh_nstat'])}"
        print(f"  {case:<8}{format_elapsed(data[case]['osh_s'])}{sh_label}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
