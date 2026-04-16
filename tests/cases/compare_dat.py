#!/usr/bin/env python3
"""compare_dat.py — shape comparison of ASCII scored-output files.

Usage:
    python3 compare_dat.py <actual> <expected> [--rtol FLOAT]
                           [--coord-cols INT] [--label STR]

Lines beginning with '#' are skipped (timestamps, version strings, etc.).
The first --coord-cols columns (default 3: X Y Z) are compared directly —
they are deterministic grid positions independent of nstat.  The remaining
columns are scored quantities that scale linearly with nstat; each is
normalised by the column's own L1 norm before comparison so that runs with
different nstat can be compared.  The shape check passes when:

    |a_norm - e_norm| <= rtol * max(|a_norm|, |e_norm|)

for every normalised value pair, with an absolute floor (ABS_FLOOR) so that
near-zero bins are treated as zero.
"""

import sys
import argparse


ABS_FLOOR = 1e-10  # normalised values smaller than this are treated as zero


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("actual", help="actual output file")
    p.add_argument("expected", help="reference file")
    p.add_argument("--rtol", type=float, default=0.02, help="relative tolerance on normalised values (default: 0.02 = 2%%)")
    p.add_argument(
        "--coord-cols", type=int, default=3, help="number of leading coordinate columns, not normalised (default: 3)"
    )
    p.add_argument("--label", default="", help="short label for error messages")
    return p.parse_args()


def data_lines(path):
    """Return stripped non-comment, non-empty lines."""
    with open(path) as fh:
        for line in fh:
            stripped = line.strip()
            if stripped and not stripped.startswith("#"):
                yield stripped


def parse_rows(path):
    """Return list of float rows."""
    rows = []
    for line in data_lines(path):
        try:
            rows.append([float(v) for v in line.split()])
        except ValueError as exc:
            print(f"parse error in {path}: {exc}", file=sys.stderr)
            sys.exit(1)
    return rows


def l1_norms(rows, first_col):
    """L1 norm of each column starting at first_col."""
    if not rows:
        return []
    ncols = len(rows[0])
    return [sum(abs(row[c]) for row in rows) for c in range(first_col, ncols)]


def normalise(rows, first_col, norms):
    """Return rows with score columns divided by norms (in-place copy)."""
    result = []
    for row in rows:
        new_row = list(row[:first_col])
        for j, norm in enumerate(norms):
            val = row[first_col + j]
            new_row.append(val / norm if norm > 1e-30 else val)
        result.append(new_row)
    return result


def close_enough(a, b, rtol):
    scale = max(abs(a), abs(b))
    if scale < ABS_FLOOR:
        return True
    return abs(a - b) <= rtol * scale


def main():
    args = parse_args()
    label = f"[{args.label}] " if args.label else ""

    a_rows = parse_rows(args.actual)
    e_rows = parse_rows(args.expected)

    if len(a_rows) != len(e_rows):
        print(f"FAIL {label}data line count mismatch: actual={len(a_rows)} expected={len(e_rows)}", file=sys.stderr)
        print(f"  actual  : {args.actual}", file=sys.stderr)
        print(f"  expected: {args.expected}", file=sys.stderr)
        sys.exit(1)

    if not a_rows:
        print(f"FAIL {label}no data lines found", file=sys.stderr)
        sys.exit(1)

    cc = args.coord_cols

    # Normalise score columns by each file's own L1 norm
    a_norms = l1_norms(a_rows, cc)
    e_norms = l1_norms(e_rows, cc)
    a_norm_rows = normalise(a_rows, cc, a_norms)
    e_norm_rows = normalise(e_rows, cc, e_norms)

    errors = []
    for idx, (a_row, e_row) in enumerate(zip(a_norm_rows, e_norm_rows), 1):
        if len(a_row) != len(e_row):
            errors.append(f"  data line {idx}: column count mismatch (actual={len(a_row)} expected={len(e_row)})")
            continue
        for col, (a, e) in enumerate(zip(a_row, e_row)):
            if not close_enough(a, e, args.rtol):
                rel = abs(a - e) / max(abs(a), abs(e))
                kind = "coord" if col < cc else "score (normalised)"
                errors.append(
                    f"  data line {idx} col {col + 1} ({kind}): "
                    f"actual={a:.6e} expected={e:.6e} "
                    f"rel_err={rel:.3%} > rtol={args.rtol:.3%}"
                )

    if errors:
        print(f"FAIL {label}{len(errors)} value(s) outside {args.rtol:.0%} tolerance:", file=sys.stderr)
        for msg in errors[:20]:
            print(msg, file=sys.stderr)
        if len(errors) > 20:
            print(f"  ... and {len(errors) - 20} more", file=sys.stderr)
        print(f"  actual  : {args.actual}", file=sys.stderr)
        print(f"  expected: {args.expected}", file=sys.stderr)
        sys.exit(1)

    n = len(a_rows)
    print(
        f"OK {label}{n} data lines, shape within {args.rtol:.0%} tolerance "
        f"(normalised by column L1 norm; coord cols 1-{cc} compared directly)"
    )


if __name__ == "__main__":
    main()
