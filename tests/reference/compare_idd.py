#!/usr/bin/env python3
"""compare_idd.py — absolute comparison of an integral depth-dose curve
against an external reference (SH12A, FLUKA, ...).

Unlike cases/compare_dat.py (which checks the *shape* after L1
normalisation), this tool compares the absolute per-primary energy
deposition so that global dose deficits are caught.

Both files are reduced to (z, energy density) pairs:
  - openshieldhit mesh text output: columns X Y Z E [F...]  ->  (col 3, col 4)
  - generic two-column reference:   columns Z E             ->  (col 1, col 2)
Per-bin energies are divided by the local bin width (inferred from the grid)
so that files with different binning compare correctly; the reference is
linearly interpolated onto the test grid over the overlapping z range.

Pass criteria (all must hold):
  1. integral:  |I_test / I_ref - 1|            <= --tol-integral
  2. per-bin:   max |t(z) - r(z)| / max(r)      <= --tol-bin
  3. peak:      |z_peak_test - z_peak_ref| [mm] <= --tol-peak-mm
     (peak position = centroid of bins >= 95% of the curve maximum)

Exit codes: 0 pass, 1 fail, 2 usage/parse error.
"""

import argparse
import sys


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("actual", help="openshieldhit mesh text output")
    p.add_argument("reference", help="reference IDD (SH12A/FLUKA export or osh format)")
    p.add_argument("--tol-integral", type=float, default=0.05, help="tolerance on the integral ratio (default 0.05)")
    p.add_argument(
        "--tol-bin", type=float, default=0.05, help="per-bin tolerance relative to the reference maximum (default 0.05)"
    )
    p.add_argument("--tol-peak-mm", type=float, default=1.0, help="tolerance on the peak position [mm] (default 1.0)")
    p.add_argument("--label", default="", help="short label for messages")
    return p.parse_args()


def load_curve(path):
    """Return (z_centres, density) lists; density = value / bin width."""
    zs = []
    vs = []
    with open(path) as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            cols = [float(c) for c in stripped.split()]
            if len(cols) >= 4:
                zs.append(cols[2])
                vs.append(cols[3])
            elif len(cols) >= 2:
                zs.append(cols[0])
                vs.append(cols[1])
            else:
                raise ValueError(f"{path}: cannot parse line: {stripped!r}")
    if len(zs) < 3:
        raise ValueError(f"{path}: need at least 3 data points, got {len(zs)}")
    pairs = sorted(zip(zs, vs))
    zs = [p[0] for p in pairs]
    vs = [p[1] for p in pairs]
    dz = (zs[-1] - zs[0]) / (len(zs) - 1)
    if dz <= 0.0:
        raise ValueError(f"{path}: non-increasing z grid")
    dens = [v / dz for v in vs]
    return zs, dens


def interp(zs, vs, z):
    """Linear interpolation; assumes zs ascending and z within range."""
    lo = 0
    hi = len(zs) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if zs[mid] <= z:
            lo = mid
        else:
            hi = mid
    frac = (z - zs[lo]) / (zs[hi] - zs[lo])
    return vs[lo] * (1.0 - frac) + vs[hi] * frac


def peak_position(zs, vs):
    """Centroid of bins at >= 95% of the maximum (robust against MC noise)."""
    vmax = max(vs)
    if vmax <= 0.0:
        return zs[0]
    wsum = 0.0
    zsum = 0.0
    for z, v in zip(zs, vs):
        if v >= 0.95 * vmax:
            wsum += v
            zsum += v * z
    return zsum / wsum


def main():
    args = parse_args()
    label = f"[{args.label}] " if args.label else ""

    try:
        tz, tv = load_curve(args.actual)
        rz, rv = load_curve(args.reference)
    except (OSError, ValueError) as exc:
        print(f"ERROR {label}{exc}", file=sys.stderr)
        return 2

    z_lo = max(tz[0], rz[0])
    z_hi = min(tz[-1], rz[-1])
    if z_hi <= z_lo:
        print(f"ERROR {label}no overlapping z range (test {tz[0]}..{tz[-1]}, ref {rz[0]}..{rz[-1]})", file=sys.stderr)
        return 2

    common_z = [z for z in tz if z_lo <= z <= z_hi]
    test_c = [interp(tz, tv, z) for z in common_z]
    ref_c = [interp(rz, rv, z) for z in common_z]

    dz = (common_z[-1] - common_z[0]) / (len(common_z) - 1)
    integral_test = sum(test_c) * dz
    integral_ref = sum(ref_c) * dz
    if integral_ref <= 0.0:
        print(f"ERROR {label}reference integral is not positive", file=sys.stderr)
        return 2

    integral_dev = integral_test / integral_ref - 1.0

    ref_max = max(ref_c)
    bin_dev = max(abs(t - r) for t, r in zip(test_c, ref_c)) / ref_max

    peak_test = peak_position(tz, tv)
    peak_ref = peak_position(rz, rv)
    peak_dev_mm = abs(peak_test - peak_ref) * 10.0  # cm -> mm

    failures = []
    if abs(integral_dev) > args.tol_integral:
        failures.append(f"integral ratio off by {integral_dev:+.2%} (tol {args.tol_integral:.2%})")
    if bin_dev > args.tol_bin:
        failures.append(f"max per-bin deviation {bin_dev:.2%} of ref max (tol {args.tol_bin:.2%})")
    if peak_dev_mm > args.tol_peak_mm:
        failures.append(f"peak position off by {peak_dev_mm:.2f} mm (tol {args.tol_peak_mm:.2f} mm)")

    print(
        f"{label}integral dev {integral_dev:+.2%}; "
        f"max bin dev {bin_dev:.2%} of ref max; "
        f"peak {peak_test:.3f} vs {peak_ref:.3f} cm (|d| = {peak_dev_mm:.2f} mm); "
        f"{len(common_z)} bins over z = {z_lo:.2f}..{z_hi:.2f} cm"
    )

    if failures:
        for msg in failures:
            print(f"FAIL {label}{msg}", file=sys.stderr)
        print(f"  actual   : {args.actual}", file=sys.stderr)
        print(f"  reference: {args.reference}", file=sys.stderr)
        return 1

    print(f"OK {label}IDD within tolerance of reference")
    return 0


if __name__ == "__main__":
    sys.exit(main())
