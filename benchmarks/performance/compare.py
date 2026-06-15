#!/usr/bin/env python3
"""Compare two benchmark result files produced by run_bench.py.

Prints a per-scenario delta table for transport time and throughput, plus the
phase decomposition shift.  Negative transport delta = the candidate is
faster.  Exits non-zero when any scenario regresses by more than --threshold
percent, so it can serve as a CI regression gate.

Usage:
  python3 benchmarks/performance/compare.py benchmarks/performance/baseline.json results.json
  python3 benchmarks/performance/compare.py old.json new.json --threshold 5
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

PHASES = ("fill_s", "zone_ref_s", "distance_s", "step_s", "compact_s")


def load(path: Path) -> dict:
    doc = json.loads(path.read_text())
    return {rec["id"]: rec for rec in doc.get("results", []) if rec.get("status") == "ok"}


def pct_delta(old: float, new: float):
    if not old:
        return None
    return 100.0 * (new - old) / old


def fmt_pct(value) -> str:
    return f"{value:+7.1f}%" if value is not None else "    n/a"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("baseline", type=Path, help="reference results file")
    parser.add_argument("candidate", type=Path, help="results file to compare against it")
    parser.add_argument(
        "--threshold",
        type=float,
        default=None,
        help="fail (exit 1) when transport time regresses by more than this percent",
    )
    args = parser.parse_args(argv)

    base = load(args.baseline)
    cand = load(args.candidate)
    common = [sid for sid in base if sid in cand]
    if not common:
        print("error: no common successful scenarios between the two files", file=sys.stderr)
        return 2

    base_doc = json.loads(args.baseline.read_text())
    cand_doc = json.loads(args.candidate.read_text())
    print(f"baseline : {args.baseline}  (commit {base_doc['git']['commit']}, {base_doc['machine']['cpu']})")
    print(f"candidate: {args.candidate}  (commit {cand_doc['git']['commit']}, {cand_doc['machine']['cpu']})")
    if base_doc["machine"]["cpu"] != cand_doc["machine"]["cpu"]:
        print("WARNING: results come from different CPUs; deltas are not meaningful")

    header = (
        f"{'scenario':<24} {'transp old[s]':>13} {'transp new[s]':>13} {'Δtransp':>8} "
        f"{'Δprim/s':>8} {'Δphase f/z/d/s/c [pp]':>24}"
    )
    print()
    print(header)
    print("-" * len(header))

    regressions = []
    for sid in common:
        old = base[sid]["median"]
        new = cand[sid]["median"]
        d_transport = pct_delta(old["transport_s"], new["transport_s"])
        d_pps = pct_delta(old["prim_per_s"], new["prim_per_s"])
        shifts = []
        for key in PHASES:
            old_share = 100.0 * old["phases"][key] / old["transport_s"] if old["transport_s"] else 0.0
            new_share = 100.0 * new["phases"][key] / new["transport_s"] if new["transport_s"] else 0.0
            shifts.append(f"{new_share - old_share:+.0f}")
        print(
            f"{sid:<24} {old['transport_s']:>13.3f} {new['transport_s']:>13.3f} "
            f"{fmt_pct(d_transport)} {fmt_pct(d_pps)} {'/'.join(shifts):>24}"
        )
        if args.threshold is not None and d_transport is not None and d_transport > args.threshold:
            regressions.append((sid, d_transport))

    missing = sorted(set(base) ^ set(cand))
    if missing:
        print(f"\nnot compared (present in only one file): {', '.join(missing)}")

    if regressions:
        print()
        for sid, delta in regressions:
            print(f"REGRESSION: {sid} transport time +{delta:.1f}% (> {args.threshold}%)")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
