"""Benchmark scenario registry for the OpenShieldHIT profiling harness.

Each scenario describes one reproducible openshieldhit run: which case
directory to use, how many primaries, and optional input patches.  Scenario
IDs follow the matrix in issue #138 (C = core case, S = parameter sweep).

The registry is data, not code: the runner (run_bench.py) interprets it.
Scenario fields:

  case        case directory, relative to the repository root, OR the string
              "generate:gemca_stress" for procedurally generated geometry
  nstat       default primary count (chosen so transport dominates wall time
              on a ~2000 primaries/s machine; scale with --nstat-scale)
  beam_patch  {KEYWORD: "value string"} replacements applied to beam.dat
  detect      replacement detect.dat content (None keeps the case's own)
  gen_zones   target zone count for generated geometry scenarios
  pool_capacity  rebuild knob: OSH_TRANSPORT_POOL_CAPACITY for this scenario
                 (None uses the default binary; non-None requires --allow-rebuild)
  tags        grouping labels used by --filter
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Optional, Tuple


@dataclass(frozen=True)
class Scenario:
    id: str
    title: str
    case: str
    nstat: int
    beam_patch: Dict[str, str] = field(default_factory=dict)
    detect: Optional[str] = None
    gen_zones: Optional[int] = None
    pool_capacity: Optional[int] = None
    tags: Tuple[str, ...] = ()


def _s4_detect(nfilters: int) -> str:
    """detect.dat for the S4 filter-count sweep: 1 unfiltered + N filtered pages."""
    parts = [
        "# S4: filter-count scaling on one detector",
        "Geometry Mesh",
        "    Name DepthDose",
        "    X -10.0  10.0    1",
        "    Y -10.0  10.0    1",
        "    Z   0.0  18.0  360",
        "",
    ]
    for i in range(1, nfilters + 1):
        parts += [
            "Filter",
            f"    Name F{i:02d}",
            f"    GEN = {i % 3}",
            "",
        ]
    parts += [
        "Output",
        "    Filename s4_dose.bdo",
        "    Geo DepthDose",
        "    Quantity Dose",
    ]
    for i in range(1, nfilters + 1):
        parts.append(f"    Quantity Dose F{i:02d}")
    return "\n".join(parts) + "\n"


_REPO_ROOT = Path(__file__).parents[2]
_CASES = str(Path(__file__).parent.relative_to(_REPO_ROOT) / "cases")

SCENARIOS = [
    # ---- Core scenarios --------------------------------------------------
    Scenario(
        id="c1_p100_dose",
        title="100 MeV p, water, depth-dose (clean CSDA reference)",
        case=f"{_CASES}/c1_p100_water_dose",
        nstat=20000,
        tags=("core",),
    ),
    Scenario(
        id="c2_p100_dose_fluence",
        title="100 MeV p, water, depth-dose + depth-fluence",
        case=f"{_CASES}/c2_p100_water_dose_fluence",
        nstat=20000,
        tags=("core",),
    ),
    Scenario(
        id="c3_p100_let",
        title="100 MeV p, water, dose + DLET + TLET",
        case=f"{_CASES}/c3_p100_water_let",
        nstat=20000,
        tags=("core",),
    ),
    Scenario(
        id="c4_p150_nucre_filters",
        title="150 MeV p, water, NUCRE on, filtered pages (all / GEN=0 / Z=1)",
        case=f"{_CASES}/c4_p150_nucre_filters",
        nstat=10000,
        tags=("core",),
    ),
    Scenario(
        id="c5_ct_voxel",
        title="126 MeV p, CT voxel head phantom, 3-D dose (requires LFS fixtures)",
        case="tests/cases/05_dicom_simple",
        nstat=5000,
        tags=("core", "dicom"),
    ),
    Scenario(
        id="c6_c12_let",
        title="C-12 400 MeV/u, water, dose + DLET/TLET (heavy-ion physics)",
        case=f"{_CASES}/c6_c12_water_let",
        nstat=20000,
        tags=("core",),
    ),
    Scenario(
        id="c7_p70",
        title="70 MeV p, water, dose (energy-scaling low point)",
        case=f"{_CASES}/c7_p_water_escan",
        nstat=30000,
        beam_patch={"TMAX0": "70.0     0.0"},
        tags=("core",),
    ),
    Scenario(
        id="c7_p200",
        title="200 MeV p, water, dose (energy-scaling high point)",
        case=f"{_CASES}/c7_p_water_escan",
        nstat=10000,
        beam_patch={"TMAX0": "200.0    0.0"},
        tags=("core",),
    ),
    # ---- C8 / S3: GEMCA stress, zone-count scaling -----------------------
    Scenario(
        id="c8_z10",
        title="GEMCA stress, ~10 zones",
        case="generate:gemca_stress",
        nstat=2000,
        gen_zones=10,
        tags=("core", "gemca_stress"),
    ),
    Scenario(
        id="c8_z100",
        title="GEMCA stress, ~100 zones",
        case="generate:gemca_stress",
        nstat=1000,
        gen_zones=100,
        tags=("core", "gemca_stress"),
    ),
    Scenario(
        id="c8_z1000",
        title="GEMCA stress, ~1000 zones",
        case="generate:gemca_stress",
        nstat=300,
        gen_zones=1000,
        tags=("core", "gemca_stress"),
    ),
    Scenario(
        id="c8_z5000",
        title="GEMCA stress, ~5000 zones",
        case="generate:gemca_stress",
        nstat=100,
        gen_zones=5000,
        tags=("core", "gemca_stress"),
    ),
    # ---- S1: physics toggles on the C4 base ------------------------------
    Scenario(
        id="s1_mcs_off",
        title="C4 base with MSCAT off",
        case=f"{_CASES}/c4_p150_nucre_filters",
        nstat=10000,
        beam_patch={"MSCAT": "0"},
        tags=("sweep", "s1"),
    ),
    Scenario(
        id="s1_stragg_off",
        title="C4 base with STRAGG off",
        case=f"{_CASES}/c4_p150_nucre_filters",
        nstat=10000,
        beam_patch={"STRAGG": "0"},
        tags=("sweep", "s1"),
    ),
    Scenario(
        id="s1_nucre_off",
        title="C4 base with NUCRE off",
        case=f"{_CASES}/c4_p150_nucre_filters",
        nstat=10000,
        beam_patch={"NUCRE": "0"},
        tags=("sweep", "s1"),
    ),
    # ---- S2: transport pool capacity (one rebuild per point) -------------
    # Note: the wavefront loop keeps per-slot scratch arrays (~40 B/slot) on
    # the stack, so capacity 65536 needs ~2.5 MiB of stack. Fine on Linux
    # (8 MiB default), would overflow the 1 MiB default stack on Windows.
    Scenario(
        id="s2_pool1",
        title="C1 base, pool capacity 1 (scalar reference)",
        case=f"{_CASES}/c1_p100_water_dose",
        nstat=2000,
        pool_capacity=1,
        tags=("sweep", "s2"),
    ),
    Scenario(
        id="s2_pool256",
        title="C1 base, pool capacity 256",
        case=f"{_CASES}/c1_p100_water_dose",
        nstat=20000,
        pool_capacity=256,
        tags=("sweep", "s2"),
    ),
    Scenario(
        id="s2_pool4096",
        title="C1 base, pool capacity 4096 (default)",
        case=f"{_CASES}/c1_p100_water_dose",
        nstat=20000,
        pool_capacity=4096,
        tags=("sweep", "s2"),
    ),
    Scenario(
        id="s2_pool65536",
        title="C1 base, pool capacity 65536",
        case=f"{_CASES}/c1_p100_water_dose",
        nstat=20000,
        pool_capacity=65536,
        tags=("sweep", "s2"),
    ),
    # ---- S4: filter-count scaling -----------------------------------------
    Scenario(
        id="s4_filters0",
        title="C4 base, 0 filtered pages",
        case=f"{_CASES}/c4_p150_nucre_filters",
        nstat=10000,
        detect=_s4_detect(0),
        tags=("sweep", "s4"),
    ),
    Scenario(
        id="s4_filters3",
        title="C4 base, 3 filtered pages",
        case=f"{_CASES}/c4_p150_nucre_filters",
        nstat=10000,
        detect=_s4_detect(3),
        tags=("sweep", "s4"),
    ),
    Scenario(
        id="s4_filters10",
        title="C4 base, 10 filtered pages",
        case=f"{_CASES}/c4_p150_nucre_filters",
        nstat=10000,
        detect=_s4_detect(10),
        tags=("sweep", "s4"),
    ),
]


def by_id() -> Dict[str, Scenario]:
    return {s.id: s for s in SCENARIOS}
