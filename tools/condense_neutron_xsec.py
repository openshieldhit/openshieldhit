#!/usr/bin/env python3
"""
tools/condense_neutron_xsec.py

Reads JEFF-4.0 PENDF0K files and produces a C header with condensed neutron
cross-section tables for use in osh_neutron_xsec_data.h.

Nuclide list, energy grid, and paths are configured in condense_neutron_xsec.toml
(same directory as this script).

Usage:
    python tools/condense_neutron_xsec.py [--jeff-dir DIR] [--config FILE]

MTs extracted per nuclide (missing MT → zero array):
    MT=1   sigma_total
    MT=2   sigma_elastic
    MT=4   sigma(n,n')
    MT=16  sigma(n,2n)
    MT=102 sigma(n,gamma)  capture
    MT=103 sigma(n,p)
    MT=107 sigma(n,alpha)

Units in JEFF files: E in eV, sigma in barns.
Output units:        E in MeV, sigma in mb.
"""

import argparse
import math
import re
from datetime import datetime
from pathlib import Path

import tomllib  # stdlib since Python 3.11

# ---------------------------------------------------------------------------
# Element symbol → atomic number
# ---------------------------------------------------------------------------

_ELEMENT_Z = {
    "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8,
    "F": 9, "Ne": 10, "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15,
    "S": 16, "Cl": 17, "Ar": 18, "K": 19, "Ca": 20, "Sc": 21, "Ti": 22,
    "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29,
    "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36,
    "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43,
    "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50,
    "Sb": 51, "Te": 52, "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57,
    "Ce": 58, "Pr": 59, "Nd": 60, "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64,
    "Tb": 65, "Dy": 66, "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70, "Lu": 71,
    "Hf": 72, "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78,
    "Au": 79, "Hg": 80, "Tl": 81, "Pb": 82, "Bi": 83,
    "Po": 84, "At": 85, "Rn": 86, "Fr": 87, "Ra": 88,
    "Ac": 89, "Th": 90, "Pa": 91, "U": 92, "Np": 93,
    "Pu": 94, "Am": 95, "Cm": 96, "Bk": 97, "Cf": 98,
}

_MTS = [1, 2, 4, 16, 102, 103, 107]
_MT_NAMES = {
    1: "tot", 2: "el", 4: "nn", 16: "n2n", 102: "ng", 103: "np", 107: "na",
}
_MT_COMMENTS = {
    1:   "MT=1   sigma_total [mb]",
    2:   "MT=2   sigma_elastic [mb]",
    4:   "MT=4   sigma(n,n') [mb]",
    16:  "MT=16  sigma(n,2n) [mb]",
    102: "MT=102 sigma(n,gamma) capture [mb]",
    103: "MT=103 sigma(n,p) [mb]",
    107: "MT=107 sigma(n,alpha) [mb]",
}

# ---------------------------------------------------------------------------
# ENDF-6 / PENDF parsing
# ---------------------------------------------------------------------------

_ENDF_FLOAT_RE = re.compile(r"([0-9])([\+\-])([0-9])")


def _endf_float(s: str) -> float:
    s = s.strip()
    if not s:
        return 0.0
    return float(_ENDF_FLOAT_RE.sub(r"\1e\2\3", s))


def _parse_line_fields(line: str):
    """Return (fields[0..5], mat, mf, mt) from a fixed-width ENDF line."""
    if len(line) < 75:
        return None
    fields = [line[i * 11:(i + 1) * 11] for i in range(6)]
    try:
        mat = int(line[66:70])
        mf  = int(line[70:72])
        mt  = int(line[72:75])
    except ValueError:
        return None
    return fields, mat, mf, mt


def read_mf3_section(lines: list[str], mat: int, mt: int) -> list[tuple[float, float]]:
    """
    Read (E_eV, sigma_barns) pairs from MF=3, MT=mt of the given MAT.
    Returns an empty list if the section is absent.

    MF=3 TAB1 record layout (line numbers relative to section start):
      line 1: CONT record  (ZA, AWR, 0, LR, 0, 0)          — skip
      line 2: TAB1 header  (0, 0, 0, 0, NR, NP)             — parse NR, NP
      line 3..2+ceil(NR/3): interpolation pairs              — skip (NR=1 → 1 line)
      remaining: data pairs  E1 s1 E2 s2 E3 s3 ...
    """
    # State machine: LOOK → HEADER → INTERP → DATA
    # In _LOOK, the first matching line is the CONT record; `continue` discards it
    # and the next iteration enters _HEADER.
    _LOOK, _HEADER, _INTERP, _DATA = range(4)
    state = _LOOK
    interp_lines_left = 0
    np_total = 0
    collected: list[float] = []

    for line in lines:
        parsed = _parse_line_fields(line)
        if parsed is None:
            continue
        fields, lmat, lmf, lmt = parsed

        in_our_section = (lmat == mat and lmf == 3 and lmt == mt)

        if state == _LOOK:
            if in_our_section:
                state = _HEADER  # CONT record consumed by this continue; next is header
            continue

        if not in_our_section:
            break               # left the section

        if state == _HEADER:
            # TAB1 header: field[4]=NR, field[5]=NP
            try:
                nr = int(fields[4].strip()) if fields[4].strip() else 1
                np_total = int(fields[5].strip()) if fields[5].strip() else 0
            except ValueError:
                nr, np_total = 1, 0
            interp_lines_left = max(1, (nr + 2) // 3)  # ceil(NR/3), at least 1
            state = _INTERP
            continue

        if state == _INTERP:
            interp_lines_left -= 1
            if interp_lines_left == 0:
                state = _DATA
            continue

        # state == _DATA: parse up to 3 (E, sigma) pairs per line
        for i in range(0, 6, 2):
            es = fields[i].strip()
            ss = fields[i + 1].strip()
            if es:
                collected.append(_endf_float(es))
            if ss:
                collected.append(_endf_float(ss))
        if np_total > 0 and len(collected) >= np_total * 2:
            break

    result: list[tuple[float, float]] = []
    for i in range(0, len(collected) - 1, 2):
        result.append((collected[i], collected[i + 1]))
    return result


# ---------------------------------------------------------------------------
# Nuclide → file lookup
# ---------------------------------------------------------------------------

def nuclide_to_file(jeff_dir: Path, nuclide: str) -> tuple[int, int, Path]:
    """
    Parse 'Symbol-A' (e.g. 'O-16') and find the matching PENDF file.
    Returns (Z, A, path).
    """
    sym, a_str = nuclide.split("-", 1)
    a = int(a_str)
    z = _ELEMENT_Z.get(sym)
    if z is None:
        raise ValueError(f"Unknown element symbol '{sym}' in nuclide '{nuclide}'")
    fname = f"0k-{z}-{sym}-{a}g_p.asc"
    path = jeff_dir / fname
    if not path.exists():
        raise FileNotFoundError(
            f"JEFF file not found for {nuclide}: expected {path}"
        )
    return z, a, path


# ---------------------------------------------------------------------------
# Log-log interpolation onto target grid
# ---------------------------------------------------------------------------

def loglog_interp(e_target_mev: list[float],
                  pairs: list[tuple[float, float]]) -> list[float]:
    """
    Interpolate (E_eV, sigma_barns) pairs onto e_target_mev [MeV].
    Returns sigma in mb (barns → mb: × 1000).
    Uses log-log interpolation; holds edge values outside the table range.
    Zero cross sections are clamped to a tiny positive value for log safety.
    """
    if not pairs:
        return [0.0] * len(e_target_mev)

    # Convert to MeV, mb; filter zero-sigma entries for log interpolation
    ev_list = [p[0] * 1e-6 for p in pairs]          # eV → MeV
    mb_list = [max(p[1] * 1e3, 1e-30) for p in pairs]  # barns → mb

    log_e = [math.log(e) for e in ev_list]
    log_s = [math.log(s) for s in mb_list]

    result = []
    for e_mev in e_target_mev:
        if e_mev <= ev_list[0]:
            result.append(pairs[0][1] * 1e3)   # hold first value
            continue
        if e_mev >= ev_list[-1]:
            result.append(pairs[-1][1] * 1e3)  # hold last value
            continue
        # Binary search for bracketing interval
        le = math.log(e_mev)
        lo, hi = 0, len(log_e) - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if log_e[mid] <= le:
                lo = mid
            else:
                hi = mid
        # Linear interpolation in log-log space
        t = (le - log_e[lo]) / (log_e[hi] - log_e[lo])
        log_s_interp = log_s[lo] + t * (log_s[hi] - log_s[lo])
        s_mb = math.exp(log_s_interp)
        # If original sigma was zero at endpoints, zero out
        if pairs[lo][1] == 0.0 and pairs[hi][1] == 0.0:
            s_mb = 0.0
        result.append(s_mb)
    return result


# ---------------------------------------------------------------------------
# C array formatting
# ---------------------------------------------------------------------------

def _c_array(values: list[float], n_per_row: int = 5) -> str:
    rows = []
    for i in range(0, len(values), n_per_row):
        chunk = values[i:i + n_per_row]
        rows.append("    " + ", ".join(f"{v:.6e}f" for v in chunk))
    return ",\n".join(rows)


def _c_ident(nuclide: str, mt_name: str) -> str:
    """e.g. 'O-16', 'el' → 'osh_neutron_xsec_o16_el_mb'"""
    tag = nuclide.lower().replace("-", "")
    return f"osh_neutron_xsec_{tag}_{mt_name}_mb"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    script_dir = Path(__file__).parent

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", default=str(script_dir / "condense_neutron_xsec.toml"),
                        help="Path to TOML config (default: same dir as script)")
    parser.add_argument("--jeff-dir", default=None,
                        help="Override jeff.dir from config")
    args = parser.parse_args()

    config_path = Path(args.config)
    with open(config_path, "rb") as f:
        cfg = tomllib.load(f)

    jeff_dir = Path(args.jeff_dir or cfg["jeff"]["dir"])
    if not jeff_dir.is_absolute():
        jeff_dir = Path.cwd() / jeff_dir

    grid_cfg = cfg["grid"]
    if "egrid_mev" in grid_cfg:
        egrid = list(grid_cfg["egrid_mev"])
        grid_desc = f"user-supplied irregular grid ({len(egrid)} points)"
    else:
        n = int(grid_cfg.get("npoints", 30))
        emin = float(grid_cfg["emin_mev"])
        emax = float(grid_cfg["emax_mev"])
        log_min = math.log10(emin)
        log_max = math.log10(emax)
        egrid = [10 ** (log_min + i * (log_max - log_min) / (n - 1)) for i in range(n)]
        grid_desc = f"{n} log-spaced points {emin}..{emax} MeV"

    out_path = Path(cfg["output"]["file"])
    if not out_path.is_absolute():
        out_path = Path.cwd() / out_path

    nuclides: list[str] = cfg["nuclides"]["list"]

    out_path.parent.mkdir(parents=True, exist_ok=True)

    now = datetime.now().astimezone().strftime("%Y-%m-%d %H:%M %Z")
    npoints = len(egrid)

    lines_out: list[str] = []
    lines_out.append(f"/* Auto-generated by tools/condense_neutron_xsec.py — do not edit.")
    lines_out.append(f" * Source  : JEFF-4.0 PENDF0K @ 0 K")
    lines_out.append(" * DOI     : https://doi.org/10.82555/wgw94-qcx30")
    lines_out.append(f" * Grid    : {grid_desc}")
    lines_out.append(f" * Units   : energy [MeV], cross sections [mb].")
    lines_out.append(f" * Missing MT → zero array (marked with comment).")
    lines_out.append(f" * Generated: {now}")
    lines_out.append(f" */")
    lines_out.append(f"")
    lines_out.append(f"#ifndef OSH_NEUTRON_XSEC_DATA_H")
    lines_out.append(f"#define OSH_NEUTRON_XSEC_DATA_H")
    lines_out.append(f"")
    lines_out.append(f"#define OSH_NEUTRON_XSEC_NPOINTS {npoints}")
    lines_out.append(f"")

    # Energy grid
    lines_out.append(f"static float const osh_neutron_xsec_egrid_mev[OSH_NEUTRON_XSEC_NPOINTS] = {{")
    lines_out.append(_c_array(egrid))
    lines_out.append(f"}};")
    lines_out.append(f"")

    for nuclide in nuclides:
        print(f"Processing {nuclide} ...", end=" ", flush=True)
        try:
            z, a, fpath = nuclide_to_file(jeff_dir, nuclide)
        except FileNotFoundError as e:
            print(f"SKIP ({e})")
            continue

        file_lines = fpath.read_text(encoding="ascii", errors="replace").splitlines()

        # Determine MAT number — skip the TPID tape-ID record (mat≥1, mf=0)
        # by requiring mf > 0; real section lines always have mf > 0.
        mat = None
        for fl in file_lines:
            parsed = _parse_line_fields(fl)
            if parsed is not None:
                _, lmat, lmf, _ = parsed
                if lmat > 0 and lmf > 0:
                    mat = lmat
                    break
        if mat is None:
            print(f"SKIP (could not determine MAT)")
            continue

        tag = nuclide.lower().replace("-", "")
        lines_out.append(f"/* {nuclide}: Z={z} A={a} (MAT={mat}) */")
        lines_out.append(f"#define OSH_NEUTRON_XSEC_{tag.upper()}_Z {z}")
        lines_out.append(f"#define OSH_NEUTRON_XSEC_{tag.upper()}_A {a}")

        found_mts: list[int] = []
        missing_mts: list[int] = []

        for mt in _MTS:
            pairs = read_mf3_section(file_lines, mat, mt)
            ident = _c_ident(nuclide, _MT_NAMES[mt])
            comment = _MT_COMMENTS[mt]

            if not pairs:
                missing_mts.append(mt)
                lines_out.append(
                    f"/* {comment} — MT={mt} not present in JEFF for {nuclide} */"
                )
                values = [0.0] * npoints
            else:
                found_mts.append(mt)
                values = loglog_interp(egrid, pairs)

            lines_out.append(
                f"static float const {ident}[OSH_NEUTRON_XSEC_NPOINTS] = {{"
            )
            lines_out.append(f"    /* {comment} */")
            lines_out.append(_c_array(values))
            lines_out.append(f"}};")

        lines_out.append(f"")
        present = ", ".join(f"MT={m}" for m in found_mts)
        absent  = ", ".join(f"MT={m}" for m in missing_mts)
        print(f"OK  (present: {present or 'none'}"
              + (f"; absent: {absent}" if absent else "") + ")")

    lines_out.append(f"#endif /* OSH_NEUTRON_XSEC_DATA_H */")

    out_path.write_text("\n".join(lines_out) + "\n", encoding="utf-8")
    print(f"\nWrote {out_path} ({npoints} energy points, {len(nuclides)} nuclides)")


if __name__ == "__main__":
    main()
