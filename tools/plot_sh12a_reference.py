#!/usr/bin/env python3
"""Generate SH12A-vs-OpenShieldHIT overlay plots as a multi-page PDF report."""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/openshieldhit-matplotlib")

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np


PAGE_SUFFIX_RE = re.compile(r"^(?P<stem>.+)_p(?P<page>[1-9][0-9]*)\.dat$")

PAGE_DESCRIPTIONS = {
    "NB_Z_narrow_dose": {
        1: "Depth curve: fluence, all particles",
        2: "Depth curve: dose in transport medium, all particles",
        3: "Depth curve: dose in transport medium, all protons (primary + secondary)",
        4: "Depth curve: fluence, primary protons only",
        5: "Depth curve: fluence, all protons (primary + secondary)",
    },
    "NB_Z_narrow_dose_water": {
        1: "Depth curve: fluence, all particles",
        2: "Depth curve: dose to water, all particles",
        3: "Depth curve: dose to water, all protons (primary + secondary)",
    },
    "NB_Z_narrow_LET": {
        1: "Depth curve: dose-averaged LET, all particles",
        2: "Depth curve: dose-averaged LET, primary protons only",
        3: "Depth curve: dose-averaged LET, all protons (primary + secondary)",
        4: "Depth curve: track-averaged LET, all particles",
        5: "Depth curve: track-averaged LET, primary protons only",
        6: "Depth curve: track-averaged LET, all protons (primary + secondary)",
    },
    "NB_Z_narrow_LET_water": {
        1: "Depth curve: dose-averaged LET to water, all particles",
        2: "Depth curve: dose-averaged LET to water, primary protons only",
        3: "Depth curve: dose-averaged LET to water, all protons (primary + secondary)",
        4: "Depth curve: track-averaged LET to water, all particles",
        5: "Depth curve: track-averaged LET to water, primary protons only",
        6: "Depth curve: track-averaged LET to water, all protons (primary + secondary)",
    },
    "NB_Z_narrow_QEFF": {
        1: "Depth curve: dose-averaged QEFF, all particles",
        2: "Depth curve: dose-averaged QEFF, primary protons only",
        3: "Depth curve: dose-averaged QEFF, all protons (primary + secondary)",
        4: "Depth curve: track-averaged QEFF, all particles",
        5: "Depth curve: track-averaged QEFF, primary protons only",
        6: "Depth curve: track-averaged QEFF, all protons (primary + secondary)",
    },
    "NB_target_diff": {
        1: "Target differential fluence: all particles, scored versus DEDX",
        2: "Target differential fluence: primary protons only, scored versus LET",
        3: "Target differential fluence in Si: all particles, scored versus DEDX",
        4: "Target differential fluence in Si: primary protons only, scored versus DEDX",
    },
    "NB_target_water_diff": {
        1: "Target differential fluence in water: all particles, scored versus DEDX",
        2: "Target differential fluence in water: primary protons only, scored versus LET",
    },
}

PLOT_METADATA = {
    "NB_Z_narrow_dose": {
        1: {"x": "Depth z [cm]", "y": "Fluence [1/cm^2/prim]"},
        2: {"x": "Depth z [cm]", "y": "Dose in transport medium [MeV/g/prim]"},
        3: {"x": "Depth z [cm]", "y": "Dose in transport medium [MeV/g/prim]"},
        4: {"x": "Depth z [cm]", "y": "Fluence [1/cm^2/prim]"},
        5: {"x": "Depth z [cm]", "y": "Fluence [1/cm^2/prim]"},
    },
    "NB_Z_narrow_dose_water": {
        1: {"x": "Depth z [cm]", "y": "Fluence [1/cm^2/prim]"},
        2: {"x": "Depth z [cm]", "y": "Dose to water [MeV/g/prim]"},
        3: {"x": "Depth z [cm]", "y": "Dose to water [MeV/g/prim]"},
    },
    "NB_Z_narrow_LET": {
        1: {"x": "Depth z [cm]", "y": "LET [MeV/cm]"},
        2: {"x": "Depth z [cm]", "y": "LET [MeV/cm]"},
        3: {"x": "Depth z [cm]", "y": "LET [MeV/cm]"},
        4: {"x": "Depth z [cm]", "y": "LET [MeV/cm]"},
        5: {"x": "Depth z [cm]", "y": "LET [MeV/cm]"},
        6: {"x": "Depth z [cm]", "y": "LET [MeV/cm]"},
    },
    "NB_Z_narrow_LET_water": {
        1: {"x": "Depth z [cm]", "y": "LET to water [MeV/cm]"},
        2: {"x": "Depth z [cm]", "y": "LET to water [MeV/cm]"},
        3: {"x": "Depth z [cm]", "y": "LET to water [MeV/cm]"},
        4: {"x": "Depth z [cm]", "y": "LET to water [MeV/cm]"},
        5: {"x": "Depth z [cm]", "y": "LET to water [MeV/cm]"},
        6: {"x": "Depth z [cm]", "y": "LET to water [MeV/cm]"},
    },
    "NB_Z_narrow_QEFF": {
        1: {"x": "Depth z [cm]", "y": "QEFF [1]"},
        2: {"x": "Depth z [cm]", "y": "QEFF [1]"},
        3: {"x": "Depth z [cm]", "y": "QEFF [1]"},
        4: {"x": "Depth z [cm]", "y": "QEFF [1]"},
        5: {"x": "Depth z [cm]", "y": "QEFF [1]"},
        6: {"x": "Depth z [cm]", "y": "QEFF [1]"},
    },
    "NB_target_diff": {
        1: {"x": "DEDX in transport medium [MeV/cm]", "y": "Differential fluence [1/cm^2/prim]"},
        2: {"x": "LET in transport medium [MeV/cm]", "y": "Differential fluence [1/cm^2/prim]"},
        3: {"x": "DEDX in Si [MeV/cm]", "y": "Differential fluence [1/cm^2/prim]"},
        4: {"x": "DEDX in Si [MeV/cm]", "y": "Differential fluence [1/cm^2/prim]"},
    },
    "NB_target_water_diff": {
        1: {"x": "DEDX in water [MeV/cm]", "y": "Differential fluence [1/cm^2/prim]"},
        2: {"x": "LET in water [MeV/cm]", "y": "Differential fluence [1/cm^2/prim]"},
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference_dir", type=Path, help="directory with imported SH12A *.dat fixtures")
    parser.add_argument("run_dir", type=Path, help="directory with OpenShieldHIT outputs")
    parser.add_argument("filenames", nargs="*", help="specific files to overlay; defaults to all reference *.dat files")
    parser.add_argument(
        "--normalize",
        choices=("none", "max"),
        default="none",
        help="optional y-normalization per curve",
    )
    parser.add_argument(
        "--logx",
        action="store_true",
        help="compatibility flag; differential plots are shown on log-log axes by default",
    )
    parser.add_argument("--save", type=Path, help="write a multi-page PDF report; defaults to <run_dir>/<reference_dir.name>_overlay.pdf")
    parser.add_argument("--show", action="store_true", help="also show the figures interactively")
    return parser.parse_args()


def load_table(path: Path) -> np.ndarray:
    data = np.loadtxt(path, comments="#")
    data = np.atleast_2d(data)
    if data.shape[1] < 2:
        raise ValueError(f"{path} does not contain at least two numeric columns")
    return data


def load_reference_curve(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
    data = load_table(path)
    x = data[:, 0]
    y = data[:, 1]
    err = data[:, 2] if data.shape[1] >= 3 else None
    return x, y, err


def resolve_run_path(run_dir: Path, requested_name: str) -> tuple[Path, str | None]:
    path = run_dir / requested_name
    match = PAGE_SUFFIX_RE.match(requested_name)

    if path.exists():
        return path, None
    if match and match.group("page") == "1":
        fallback = run_dir / f"{match.group('stem')}.dat"
        if fallback.exists():
            return fallback, f"using {fallback.name} because convertmc produced a single differential table"
    raise FileNotFoundError(f"missing OpenShieldHIT output: {path}")


def load_run_curve(path: Path, ref_x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    data = load_table(path)
    if data.shape[0] == 1:
        values = data[0]
        if values.shape[0] != ref_x.shape[0]:
            raise ValueError(f"{path} has {values.shape[0]} columns, expected {ref_x.shape[0]} bins")
        return ref_x, values
    x = data[:, -2]
    y = data[:, -1]
    return x, y


def normalize(values: np.ndarray, mode: str) -> np.ndarray:
    if mode == "none":
        return values
    scale = np.max(np.abs(values))
    if scale <= 0.0:
        return values
    return values / scale


def describe_reference_file(filename: str) -> str:
    match = PAGE_SUFFIX_RE.match(filename)
    if not match:
        return filename
    stem = match.group("stem")
    page = int(match.group("page"))
    return PAGE_DESCRIPTIONS.get(stem, {}).get(page, filename)


def choose_filenames(reference_dir: Path, requested: list[str]) -> list[str]:
    if requested:
        return requested
    return sorted(path.name for path in reference_dir.glob("*.dat"))


def default_report_path(reference_dir: Path, run_dir: Path) -> Path:
    return run_dir / f"{reference_dir.name}_overlay.pdf"


def axis_label(filename: str) -> str:
    meta = plot_metadata(filename)
    return meta.get("x", "Axis")


def is_differential_plot(filename: str) -> bool:
    return "diff" in filename


def plot_metadata(filename: str) -> dict[str, str]:
    match = PAGE_SUFFIX_RE.match(filename)
    if not match:
        return {}
    stem = match.group("stem")
    page = int(match.group("page"))
    return PLOT_METADATA.get(stem, {}).get(page, {})


def positive_mask(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    return (x > 0.0) & (y > 0.0)


def plot_one(filename: str,
             ref_path: Path,
             run_dir: Path,
             normalize_mode: str,
             logx: bool) -> tuple[plt.Figure, list[str]]:
    warnings: list[str] = []
    ref_x, ref_y, ref_err = load_reference_curve(ref_path)
    title = describe_reference_file(filename)

    fig, ax = plt.subplots(figsize=(8.2, 5.0))
    fig.suptitle(title, fontsize=13, y=0.98)
    ax.set_title(filename, fontsize=10)

    ref_y_plot = normalize(ref_y, normalize_mode)
    ref_err_plot = ref_err
    if ref_err is not None and normalize_mode == "max":
        scale = np.max(np.abs(ref_y))
        if scale > 0.0:
            ref_err_plot = ref_err / scale

    if is_differential_plot(filename):
        ref_mask = positive_mask(ref_x, ref_y_plot)
        ref_x_plot = ref_x[ref_mask]
        ref_y_plot2 = ref_y_plot[ref_mask]
        if ref_err_plot is not None:
            ref_err_plot = ref_err_plot[ref_mask]
        if ref_err_plot is not None:
            ax.errorbar(ref_x_plot, ref_y_plot2, yerr=ref_err_plot, fmt=".", ms=2.5, lw=0.8, alpha=0.8, label="SH12A")
        else:
            ax.plot(ref_x_plot, ref_y_plot2, ".", ms=3.0, alpha=0.8, label="SH12A")
    else:
        if ref_err_plot is not None:
            ax.errorbar(ref_x, ref_y_plot, yerr=ref_err_plot, fmt=".", ms=2.5, lw=0.8, alpha=0.8, label="SH12A")
        else:
            ax.plot(ref_x, ref_y_plot, ".", ms=3.0, alpha=0.8, label="SH12A")

    try:
        run_path, fallback_note = resolve_run_path(run_dir, filename)
        if fallback_note:
            warnings.append(f"{filename}: {fallback_note}")
        run_x, run_y = load_run_curve(run_path, ref_x)
        run_y_plot = normalize(run_y, normalize_mode)
        if is_differential_plot(filename):
            run_mask = positive_mask(run_x, run_y_plot)
            ax.plot(run_x[run_mask], run_y_plot[run_mask], "-", lw=1.5, label="OpenShieldHIT")
        else:
            ax.plot(run_x, run_y_plot, "-", lw=1.5, label="OpenShieldHIT")
    except (FileNotFoundError, ValueError) as exc:
        warnings.append(f"{filename}: {exc}")
        ax.text(0.5,
                0.08,
                "OpenShieldHIT output missing or not directly plottable\nsee terminal warnings for details",
                transform=ax.transAxes,
                ha="center",
                va="bottom",
                fontsize=9,
                bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "0.7"})

    if logx or is_differential_plot(filename):
        ax.set_xscale("log")
    if is_differential_plot(filename):
        ax.set_yscale("log")
    ax.set_xlabel(axis_label(filename))
    ylabel = plot_metadata(filename).get("y", "Value")
    if normalize_mode != "none":
        ylabel = f"Normalized {ylabel}"
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    return fig, warnings


def main() -> int:
    args = parse_args()
    filenames = choose_filenames(args.reference_dir, args.filenames)
    report_path = args.save if args.save is not None else default_report_path(args.reference_dir, args.run_dir)
    all_warnings: list[str] = []

    report_path.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(report_path) as pdf:
        for filename in filenames:
            ref_path = args.reference_dir / filename
            if not ref_path.exists():
                raise FileNotFoundError(f"missing reference fixture: {ref_path}")
            fig, warnings = plot_one(filename, ref_path, args.run_dir, args.normalize, args.logx)
            all_warnings.extend(warnings)
            pdf.savefig(fig)
            if args.show:
                fig.show()
            plt.close(fig)

    print(f"Wrote multi-page PDF report to {report_path}")
    if all_warnings:
        print("Warnings:")
        for warning in all_warnings:
            print(f"  - {warning}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
