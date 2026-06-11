#!/usr/bin/env python3
"""Generate the C8 GEMCA stress-test case: a lattice of SPH bodies in water.

Emits a complete, self-contained case directory (geo.dat / mat.dat / beam.dat /
detect.dat) whose zone count is parameterized.  The geometry is a water RCC
target inside the usual vacuum-gap / blackhole-shell pair, with a 3-D cubic
lattice of PMMA spheres carved out of the water matrix:

  - each sphere is its own zone (``+sNNNN``),
  - the water matrix zone is ``+1 -s0001 -s0002 ...`` (one negative term per
    sphere), which is exactly the long boolean expression that stresses
    ``eval_distance`` and zone lookup,
  - vacuum and blackhole complete the universe.

Total zones = nspheres + 3.  The requested zone count is rounded down to the
nearest count realizable with a full n x n x n lattice plus the 3 shell zones
(so ``--zones 1000`` gives a 9x9x9 lattice = 732 zones; pass the exact
realizable count if you need to hit a number precisely).

Usage:
  python3 tools/bench/gen_gemca_stress.py --zones 1000 --out /tmp/c8_z1000
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Water target: RCC along +z, radius 10 cm, z in [0, 20].  The lattice is
# inscribed in a box well inside the cylinder so spheres never poke out.
TARGET_RADIUS = 10.0
TARGET_LENGTH = 20.0
LATTICE_HALF_XY = 6.0
LATTICE_Z_LO = 2.0
LATTICE_Z_HI = 18.0

BEAM_DAT = """\
RNDSEED        89736501            # Fixed seed: benchmark runs must be reproducible
PRIMARY        1        1          # Proton
TMAX0          150.0    0.0        # Incident energy (MeV/nucleon)
BEAMSIGMA      -6.0     -6.0       # 12x12 cm square field: histories sample the lattice
BEAMPOS        0.0      0.0  -2.0  # Beam XYZ start position [cm]
NSTAT          2000    -1          # Overridden by the harness via -n
DELTAE         0.005               # Max relative energy loss per step
DEMIN          0.025               # Minimum energy step (MeV/nucleon)
STRAGG         1                   # Gaussian straggling
MSCAT          2                   # Moliere multiple scattering
NUCRE          0                   # Geometry stress only: keep physics minimal
"""

MAT_DAT = """\
Material Water
    ICRU 276

Material PMMA
    ICRU 223
"""

DETECT_DAT = """\
# C8: GEMCA stress - depth-dose through the sphere lattice
Geometry Mesh
    Name DepthDose
    X -10.0  10.0    1
    Y -10.0  10.0    1
    Z   0.0  20.0  200

Output
    Filename gemca_stress_dose.bdo
    Geo DepthDose
    Quantity Dose
"""


def lattice_side_for_zones(nzones: int) -> int:
    """Largest n with n^3 + 3 <= nzones (at least 1)."""
    if nzones < 4:
        raise ValueError("need at least 4 zones (1 sphere + matrix + vacuum + blackhole)")
    n = int(round((nzones - 3) ** (1.0 / 3.0)))
    while n**3 + 3 > nzones:
        n -= 1
    while (n + 1) ** 3 + 3 <= nzones:
        n += 1
    return max(n, 1)


def sphere_centers(n: int):
    """Centers of an n x n x n lattice inside the target box."""
    xs = [-LATTICE_HALF_XY + (i + 0.5) * (2.0 * LATTICE_HALF_XY / n) for i in range(n)]
    zs = [LATTICE_Z_LO + (i + 0.5) * ((LATTICE_Z_HI - LATTICE_Z_LO) / n) for i in range(n)]
    for z in zs:
        for y in xs:
            for x in xs:
                yield (x, y, z)


def sphere_radius(n: int) -> float:
    """40% of the smallest lattice pitch, so spheres never touch."""
    pitch_xy = 2.0 * LATTICE_HALF_XY / n
    pitch_z = (LATTICE_Z_HI - LATTICE_Z_LO) / n
    return 0.4 * min(pitch_xy, pitch_z)


def generate_geo(n: int) -> str:
    nspheres = n**3
    r = sphere_radius(n)
    lines = []
    lines.append(f"    0    0  GEMCA stress: {nspheres} PMMA spheres ({n}x{n}x{n} lattice) in water")
    lines.append(f"  RCC    1       0.0       0.0       0.0       0.0       0.0      {TARGET_LENGTH:.1f}")
    lines.append(f"                {TARGET_RADIUS:.1f}")
    lines.append("  RCC    2       0.0       0.0      -5.0       0.0       0.0      25.0")
    lines.append(f"                {TARGET_RADIUS:.1f}")
    lines.append("  RCC    3       0.0       0.0     -10.0       0.0       0.0      35.0")
    lines.append("                20.0")
    for i, (x, y, z) in enumerate(sphere_centers(n), start=1):
        lines.append(f"  SPH    s{i:05d}  {x:.6f} {y:.6f} {z:.6f} {r:.6f}")
    lines.append("  END")
    for i in range(1, nspheres + 1):
        lines.append(f"  z{i:05d}       +s{i:05d}")
    lines.append("  matrix       +1")
    for start in range(1, nspheres + 1, 8):
        terms = " ".join(f"-s{i:05d}" for i in range(start, min(start + 8, nspheres + 1)))
        lines.append(f"               {terms}")
    lines.append("  vac          +2     -1")
    lines.append("  bh           +3     -2")
    lines.append("  END")
    lines.append(f"  ASSIGNMAT PMMA z00001 z{nspheres:05d}")
    lines.append("  ASSIGNMAT Water matrix")
    lines.append("  ASSIGNMAT vacuum vac")
    lines.append("  ASSIGNMAT blackhole bh")
    return "\n".join(lines) + "\n"


def generate_case(nzones: int, out_dir: Path) -> int:
    """Write the case files; returns the actual zone count."""
    n = lattice_side_for_zones(nzones)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "geo.dat").write_text(generate_geo(n))
    (out_dir / "mat.dat").write_text(MAT_DAT)
    (out_dir / "beam.dat").write_text(BEAM_DAT)
    (out_dir / "detect.dat").write_text(DETECT_DAT)
    return n**3 + 3


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--zones", type=int, required=True, help="target zone count (>= 4)")
    parser.add_argument("--out", type=Path, required=True, help="output case directory")
    args = parser.parse_args(argv)
    actual = generate_case(args.zones, args.out)
    print(f"generated {args.out} with {actual} zones (requested {args.zones})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
