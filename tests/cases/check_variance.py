#!/usr/bin/env python3
"""check_variance.py — validate the standard-error columns of a VARIANCE run.

Usage:
    python3 check_variance.py <variance.dat> <plain.dat> [--coord-cols INT]
                              [--rtol FLOAT] [--label STR]

Given a scored ASCII file produced with the detect.dat VARIANCE card and the
matching file from a plain (no-variance) run of the same case, assert that:

  * the variance file has exactly twice as many scored columns as the plain file
    (each value column is followed by its paired standard-error column);
  * every value column matches the plain run within --rtol after per-column L1
    normalisation (batching only reorders the FP summation, so the physics is the
    same) — reusing the shape check from compare_dat.py;
  * every standard-error column is finite and non-negative, and is zero exactly
    where its value is zero (an empty bin has no defined error);
  * at least one standard-error entry is strictly positive (the run actually had
    >= 2 batches and produced a usable error, not an all-zero column).

Lines beginning with '#' are skipped.  Exit code 0 on success, 1 on any failure.
"""

import sys
import argparse
import math


ABS_FLOOR = 1e-10


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("variance", help="scored file from a VARIANCE run")
    p.add_argument("plain", help="scored file from a plain (no-variance) run")
    p.add_argument("--coord-cols", type=int, default=3, help="leading coordinate columns (default 3)")
    p.add_argument("--rtol", type=float, default=0.02, help="relative tolerance on normalised values")
    p.add_argument("--label", default="", help="short label for error messages")
    return p.parse_args()


def rows(path):
    out = []
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if s and not s.startswith("#"):
                out.append([float(x) for x in s.split()])
    return out


def fail(label, msg):
    where = f"[{label}] " if label else ""
    print(f"{where}check_variance FAIL: {msg}", file=sys.stderr)
    sys.exit(1)


def l1(col):
    return sum(abs(v) for v in col)


def main():
    a = parse_args()
    var = rows(a.variance)
    plain = rows(a.plain)
    cc = a.coord_cols

    if not var or not plain:
        fail(a.label, "one of the inputs has no data rows")
    if len(var) != len(plain):
        fail(a.label, f"row count differs: variance {len(var)} vs plain {len(plain)}")

    nplain_scored = len(plain[0]) - cc
    nvar_scored = len(var[0]) - cc
    if nplain_scored <= 0:
        fail(a.label, f"plain file has no scored columns (ncols={len(plain[0])}, coord-cols={cc})")
    if nvar_scored != 2 * nplain_scored:
        fail(
            a.label,
            f"variance file should have 2x the scored columns: got {nvar_scored}, expected {2 * nplain_scored}",
        )

    # Value columns of the variance file are the even scored columns (0,2,4,...);
    # error columns are the odd ones (1,3,5,...).
    nq = nplain_scored
    any_positive = False
    # Per-quantity L1 norms for the shape comparison of the value columns.
    for q in range(nq):
        vcol_var = [r[cc + 2 * q] for r in var]
        ecol_var = [r[cc + 2 * q + 1] for r in var]
        vcol_plain = [r[cc + q] for r in plain]

        nv = l1(vcol_var)
        npl = l1(vcol_plain)
        for i in range(len(var)):
            # Error column sanity.
            e = ecol_var[i]
            if not math.isfinite(e):
                fail(a.label, f"quantity {q} row {i}: error is not finite ({e})")
            if e < 0.0:
                fail(a.label, f"quantity {q} row {i}: error is negative ({e})")
            if vcol_var[i] == 0.0 and e != 0.0:
                fail(a.label, f"quantity {q} row {i}: non-zero error on a zero-value bin")
            if e > 0.0:
                any_positive = True
            # Value shape match against the plain run.
            av = vcol_var[i] / nv if nv > 0 else 0.0
            ev = vcol_plain[i] / npl if npl > 0 else 0.0
            if abs(av) < ABS_FLOOR and abs(ev) < ABS_FLOOR:
                continue
            if abs(av - ev) > a.rtol * max(abs(av), abs(ev)):
                fail(a.label, f"quantity {q} row {i}: value differs from plain run ({av} vs {ev})")

    if not any_positive:
        fail(a.label, "every standard-error entry is zero — variance produced no usable error (B < 2?)")

    print(f"check_variance OK: {a.variance} ({nq} quantities, error columns present and sane)")


if __name__ == "__main__":
    main()
