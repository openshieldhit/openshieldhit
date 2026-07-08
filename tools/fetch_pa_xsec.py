#!/usr/bin/env python3
"""
tools/fetch_pa_xsec.py

Fetch and condense proton-nucleus cross-section reference data (issue #277)
into provenance-headed tables committed under tests/reference/xsec/.

Sources (all retrieved from IAEA Nuclear Data Services):
  * EXFOR Web-API (nds.iaea.org/exfor): measured nonelastic sigma_R points
    ((P,NON),SIG) and integrated nuclear elastic points ((P,EL),SIG), for the
    target isotope and the natural element.
  * ENDF/B-VIII.0 proton sublibrary (LA150h lineage — the evaluation behind
    the ICRU 63 nonelastic tables): MF3 MT5 (nonelastic).
  * TENDL-2023 proton sublibrary: same (adopts LA150 for C-12/O-16).

ENDF-6 parsing is reused from tools/condense_neutron_xsec.py.

Output table format (whitespace-separated, '#' comments):
    source  E_MeV  sigma_mb  err_mb
where source is one of  exfor:<DatasetID>, exfor-el:<DatasetID>,
endfb8-nonel, tendl23-nonel.

Usage:
    python tools/fetch_pa_xsec.py [--out-dir tests/reference/xsec]
                                  [--cache-dir ~/.cache/openshieldhit-xsec]
"""

from __future__ import annotations

import argparse
import csv
import io
import sys
import urllib.request
import zipfile
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from condense_neutron_xsec import read_mf3_section  # noqa: E402

NDS = "https://www-nds.iaea.org/public/download-endf"
EXFOR_API = "https://nds.iaea.org/exfor"

# MT section extracted from the proton sublibraries: MT5 (all nonelastic,
# lumped).  MF3/MT2 is deliberately NOT extracted: for incident charged
# particles the integrated "elastic" entry is a Coulomb-interference
# bookkeeping convention of the evaluation (LA150 carries the *identical*
# curve for C-12 and O-16, e.g. 61.1 mb at 30 MeV, 9.5 mb at 150 MeV), not a
# transport-usable nuclear elastic cross section.  Elastic calibration uses
# the measured integrated points (below), the SH12A-implied removal, and
# optical-model literature instead.
MT_NONELASTIC = 5

E_MIN_MEV = 1.0
E_MAX_MEV = 300.0
E_MIN_ELASTIC_MEV = 20.0

# Plausibility floor for (P,EL),SIG points: a genuine integrated nuclear
# elastic cross section on C/O above 20 MeV is O(100 mb); EXFOR also returns
# state-specific/partial quantities under the same reaction code (e.g. the
# few-mb Karban p+O-16 excitation-function set), which this floor rejects.
SIGMA_MIN_ELASTIC_MB = 30.0

TARGETS = {
    "O16": {
        "exfor_targets": ("O-16", "O-0"),
        "endfb8": f"{NDS}/ENDF-B-VIII.0/p/p_0825_8-O-16.zip",
        "tendl23": f"{NDS}/TENDL-2023/p/p_008-O-16_0825.zip",
    },
    "C12": {
        "exfor_targets": ("C-12", "C-0"),
        "endfb8": f"{NDS}/ENDF-B-VIII.0/p/p_0625_6-C-12.zip",
        "tendl23": f"{NDS}/TENDL-2023/p/p_006-C-12_0625.zip",
    },
}


def http_get(url: str) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "openshieldhit-fetch-pa-xsec/1.0"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        return resp.read()


def fetch_cached(url: str, cache_dir: Path) -> bytes:
    cache_dir.mkdir(parents=True, exist_ok=True)
    name = url.rsplit("/", 1)[-1] if "/" in url else url
    if "?" in name:
        name = name.replace("?", "_").replace("&", "_").replace("=", "-").replace(",", "_")
    path = cache_dir / name
    if path.exists():
        return path.read_bytes()
    data = http_get(url)
    path.write_bytes(data)
    return data


# ---------------------------------------------------------------------------
# EXFOR
# ---------------------------------------------------------------------------


def exfor_list_datasets(target: str, reaction: str, cache_dir: Path) -> list[dict]:
    """Dataset metadata rows for (P,<reaction>),SIG on the given EXFOR target."""
    url = f"{EXFOR_API}/x4list?Target={target}&Reaction=p,{reaction}&Quantity=SIG&csv"
    raw = fetch_cached(url, cache_dir / "exfor").decode("utf-8", errors="replace")
    rows = list(csv.DictReader(io.StringIO(raw)))
    return [r for r in rows if r.get("DatasetID")]


def exfor_get_points(
    dataset_id: str, cache_dir: Path, e_min: float, sigma_min_mb: float = 0.0
) -> list[tuple[float, float, float]]:
    """(E_MeV, sigma_mb, err_mb) points of one EXFOR dataset (err 0 if absent)."""
    url = f"{EXFOR_API}/x4get?DatasetID={dataset_id}&op=csv&plus=2"
    raw = fetch_cached(url, cache_dir / "exfor").decode("utf-8", errors="replace")
    points = []
    for row in csv.DictReader(io.StringIO(raw)):
        try:
            e_mev = float(row["x2(eV)"]) / 1.0e6
            sigma_mb = float(row["y"]) * 1000.0
        except (KeyError, TypeError, ValueError):
            continue
        try:
            err_mb = float(row["dy"]) * 1000.0
        except (KeyError, TypeError, ValueError):
            err_mb = 0.0
        if e_min <= e_mev <= E_MAX_MEV and sigma_mb > sigma_min_mb:
            points.append((e_mev, sigma_mb, err_mb))
    return points


# ---------------------------------------------------------------------------
# Evaluated libraries (ENDF-6)
# ---------------------------------------------------------------------------


def endf_lines_from_zip(url: str, cache_dir: Path) -> list[str]:
    blob = fetch_cached(url, cache_dir / "endf")
    with zipfile.ZipFile(io.BytesIO(blob)) as zf:
        names = [n for n in zf.namelist() if n.lower().endswith((".dat", ".endf", ".tendl"))]
        if not names:
            names = zf.namelist()
        return zf.read(names[0]).decode("utf-8", errors="replace").splitlines()


def endf_mat_number(lines: list[str]) -> int:
    """MAT of the first real section (MF > 0), skipping the tape-header line."""
    for line in lines:
        if len(line) >= 75:
            try:
                mat = int(line[66:70])
                mf = int(line[70:72])
            except ValueError:
                continue
            if mat > 0 and mf > 0:
                return mat
    raise ValueError("no MAT number found")


def evaluated_curve(lines: list[str], mt: int) -> list[tuple[float, float, float]]:
    """(E_MeV, sigma_mb, 0.0) from MF3/MT, restricted to the energy window."""
    mat = endf_mat_number(lines)
    pairs = read_mf3_section(lines, mat, mt)
    return [
        (e_ev / 1.0e6, sig_b * 1000.0, 0.0) for e_ev, sig_b in pairs if E_MIN_MEV <= e_ev / 1.0e6 <= E_MAX_MEV and sig_b > 0.0
    ]


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------


def write_table(out_path: Path, label: str, meta: list[dict], blocks: dict) -> None:
    lines = []
    lines.append(f"# p+{label} cross-section reference data (issue #277)")
    lines.append(f"# Retrieved {date.today().isoformat()} by tools/fetch_pa_xsec.py from IAEA NDS.")
    lines.append("#")
    lines.append("# Evaluated curves (MF3/MT5, all nonelastic lumped):")
    lines.append("#   endfb8-nonel  : ENDF/B-VIII.0 proton sublibrary (LA150h lineage; the")
    lines.append("#                   evaluation behind the ICRU 63 nonelastic tables)")
    lines.append("#   tendl23-nonel : TENDL-2023 proton sublibrary (adopts LA150 for C-12/O-16)")
    lines.append("#   No evaluated elastic: MF3/MT2 of charged-particle files is a Coulomb-")
    lines.append("#   interference bookkeeping convention (LA150 carries the identical curve")
    lines.append("#   for C-12 and O-16), not a transport-usable nuclear elastic sigma.")
    lines.append("#")
    lines.append("# EXFOR datasets: exfor:<ID> = (P,NON),SIG;  exfor-el:<ID> = (P,EL),SIG")
    lines.append("#   (integrated nuclear elastic as published, >= 20 MeV only):")
    for m in meta:
        npts = m.get("nPts", "?")
        lines.append(
            f"#   {m['DatasetID']:<10} {m.get('Target', ''):<5} {m.get('year1', ''):<5}"
            f" {m.get('author1', ''):<22} {m.get('reference1', '')}  ({npts} pts as published)"
        )
    lines.append("#")
    lines.append("# columns: source  E_MeV  sigma_mb  err_mb   (err 0.0 = not given)")
    for source in sorted(blocks):
        for e, s, err in blocks[source]:
            lines.append(f"{source:<16} {e:10.4f} {s:10.3f} {err:8.3f}")
    out_path.write_text("\n".join(lines) + "\n")
    print(f"wrote {out_path}  ({sum(len(p) for p in blocks.values())} rows)")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", type=Path, default=Path("tests/reference/xsec"))
    ap.add_argument("--cache-dir", type=Path, default=Path.home() / ".cache" / "openshieldhit-xsec")
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    for label, cfg in TARGETS.items():
        blocks: dict = {}
        meta: list[dict] = []
        for tgt in cfg["exfor_targets"]:
            for reaction, prefix, e_min, s_min in (
                ("non", "exfor", E_MIN_MEV, 0.0),
                ("el", "exfor-el", E_MIN_ELASTIC_MEV, SIGMA_MIN_ELASTIC_MB),
            ):
                for m in exfor_list_datasets(tgt, reaction, args.cache_dir):
                    points = exfor_get_points(m["DatasetID"], args.cache_dir, e_min, s_min)
                    if not points:
                        continue
                    meta.append(m)
                    blocks[f"{prefix}:{m['DatasetID']}"] = sorted(points)
        for lib_key, lib_tag in (("endfb8", "endfb8"), ("tendl23", "tendl23")):
            lines = endf_lines_from_zip(cfg[lib_key], args.cache_dir)
            blocks[f"{lib_tag}-nonel"] = evaluated_curve(lines, MT_NONELASTIC)
        write_table(args.out_dir / f"p_{label}.txt", label, meta, blocks)


if __name__ == "__main__":
    main()
