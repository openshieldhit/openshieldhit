#!/usr/bin/env python3
"""One-command benchmark harness for OpenShieldHIT (issue #138, Phase 0).

Runs the scenario matrix from scenarios.py against a release/prof build,
collects the per-run JSON profiles emitted by ``openshieldhit --profile``,
wraps them with case + machine + commit metadata, and writes one aggregated
machine-readable results file plus a human-readable summary table.

Typical use:

  # build the release preset and run the core matrix, 3 repeats each
  python3 tools/bench/run_bench.py --build

  # everything, including sweeps that need per-point rebuilds (S2)
  python3 tools/bench/run_bench.py --build --filter all --allow-rebuild

  # a single scenario with profiler wrappers (requires prof preset tools)
  python3 tools/bench/run_bench.py --filter c1_p100_dose --preset prof --build \\
      --callgrind --perf

  # freeze a baseline
  python3 tools/bench/run_bench.py --build --filter all --allow-rebuild \\
      --output tools/bench/baseline.json

Methodology guard-rails baked in (see issue #138):
  - refuses to measure debug builds,
  - one un-recorded warm-up run per scenario,
  - >= 3 repeats by default, median + spread reported,
  - setup (parse + compile) and transport time reported separately; a warning
    is printed when transport is below 95% of wall time,
  - every results file embeds machine metadata (CPU, ISA flags, cores, RAM,
    OS, compiler) and the git commit.
"""

from __future__ import annotations

import argparse
import datetime
import json
import math
import platform
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from scenarios import SCENARIOS, Scenario  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]

PRESET_BUILD_DIRS = {
    "release": "build",
    "prof": "build_prof",
    "relwithdebinfo": "build-rel",
}

TRANSPORT_FRACTION_WARN = 0.95
PHASE_SUM_TOLERANCE = 0.02


# ---- Machine / commit metadata ---------------------------------------------


def _read_first_match(path: str, pattern: str):
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fp:
            for line in fp:
                m = re.match(pattern, line)
                if m:
                    return m.group(1).strip()
    except OSError:
        return None
    return None


def machine_metadata() -> dict:
    cpu = _read_first_match("/proc/cpuinfo", r"model name\s*:\s*(.+)") or platform.processor()
    flags_line = _read_first_match("/proc/cpuinfo", r"(?:flags|Features)\s*:\s*(.+)") or ""
    interesting = ("sse2", "sse4_2", "avx", "avx2", "avx512f", "fma", "asimd", "neon", "sve")
    isa = sorted(set(flags_line.split()) & set(interesting))
    mem_kb = _read_first_match("/proc/meminfo", r"MemTotal:\s*(\d+)")
    return {
        "cpu": cpu or "unknown",
        "isa": isa,
        "logical_cores": __import__("os").cpu_count(),
        "ram_gb": round(int(mem_kb) / 1024.0 / 1024.0, 1) if mem_kb else None,
        "os": platform.platform(),
        "python": platform.python_version(),
    }


def git_commit() -> dict:
    def run(*args):
        try:
            return subprocess.run(["git", *args], cwd=REPO_ROOT, capture_output=True, text=True, check=True).stdout.strip()
        except (subprocess.CalledProcessError, FileNotFoundError):
            return None

    commit = run("rev-parse", "--short", "HEAD")
    dirty = bool(run("status", "--porcelain"))
    return {"commit": commit, "dirty": dirty}


# ---- Build management --------------------------------------------------------


def build_dir_for(preset: str, pool_capacity=None) -> Path:
    if pool_capacity is not None:
        return REPO_ROOT / f"build_bench_pool{pool_capacity}"
    return REPO_ROOT / PRESET_BUILD_DIRS[preset]


def binary_path(build_dir: Path) -> Path:
    return build_dir / "bin" / "openshieldhit"


def cmake_build(preset: str, pool_capacity=None) -> Path:
    """Configure + build; returns the build directory."""
    bdir = build_dir_for(preset, pool_capacity)
    cfg = [
        "cmake",
        "-S",
        str(REPO_ROOT),
        "-B",
        str(bdir),
        "-DOSH_BUILD_EXAMPLES=OFF",
        "-DCMAKE_BUILD_TYPE=" + ("RelWithDebInfo" if preset == "prof" else "Release"),
    ]
    flags = []
    if preset == "prof":
        flags.append("-fno-omit-frame-pointer")
    if pool_capacity is not None:
        flags.append(f"-DOSH_TRANSPORT_POOL_CAPACITY={pool_capacity}u")
    if flags:
        cfg.append("-DCMAKE_C_FLAGS=" + " ".join(flags))
    subprocess.run(cfg, check=True, capture_output=True, text=True)
    subprocess.run(["cmake", "--build", str(bdir), "--parallel"], check=True, capture_output=True, text=True)
    return bdir


# ---- Case preparation --------------------------------------------------------


def patch_beam(text: str, patch: dict) -> str:
    """Replace the value of existing beam.dat keywords (append when missing)."""
    lines = text.splitlines()
    seen = set()
    out = []
    for line in lines:
        tokens = line.split()
        key = tokens[0].upper() if tokens else ""
        if key in patch:
            out.append(f"{key:<14} {patch[key]}")
            seen.add(key)
        else:
            out.append(line)
    for key, value in patch.items():
        if key not in seen:
            out.append(f"{key:<14} {value}")
    return "\n".join(out) + "\n"


def prepare_case(scenario: Scenario, work_root: Path) -> Path:
    """Materialize the scenario's case directory under work_root.

    Unpatched cases are used in place: case inputs may reference sibling
    files (e.g. tests/cases/05_dicom_simple points at ../../fixtures) that a
    copy would orphan, and inputs are only read, never written.
    """
    case_dir = work_root / scenario.id
    if scenario.case.startswith("generate:"):
        from gen_gemca_stress import generate_case

        generate_case(scenario.gen_zones, case_dir)
    elif not scenario.beam_patch and scenario.detect is None:
        return REPO_ROOT / scenario.case
    else:
        shutil.copytree(REPO_ROOT / scenario.case, case_dir)
    if scenario.beam_patch:
        beam = case_dir / "beam.dat"
        beam.write_text(patch_beam(beam.read_text(), scenario.beam_patch))
    if scenario.detect is not None:
        (case_dir / "detect.dat").write_text(scenario.detect)
    return case_dir


# ---- Running -------------------------------------------------------------------


def run_once(binary: Path, case_dir: Path, nstat: int, out_dir: Path, extra_argv=None) -> dict:
    """One openshieldhit run; returns the parsed profile JSON."""
    out_dir.mkdir(parents=True, exist_ok=True)
    profile_path = out_dir / "profile.json"
    argv = list(extra_argv or []) + [
        str(binary),
        "--workdir",
        str(case_dir),
        "-n",
        str(nstat),
        "-o",
        str(out_dir),
        "--profile",
        str(profile_path),
    ]
    proc = subprocess.run(argv, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"run failed (rc={proc.returncode}): {' '.join(argv)}\n{proc.stderr.strip()[-2000:]}")
    return json.loads(profile_path.read_text())


def aggregate(runs: list) -> dict:
    """Median + spread over repeated runs of one scenario."""

    def med(key_fn):
        return statistics.median(key_fn(r) for r in runs)

    transport = [r["transport_s"] for r in runs]
    median_transport = statistics.median(transport)
    spread_pct = 100.0 * (max(transport) - min(transport)) / median_transport if median_transport else 0.0
    wall = [r["setup_parse_s"] + r["setup_compile_s"] + r["run_s"] + r["save_s"] for r in runs]
    agg = {
        "transport_s": median_transport,
        "transport_spread_pct": spread_pct,
        "prim_per_s": med(lambda r: r["prim_per_s"]),
        "setup_s": med(lambda r: r["setup_parse_s"] + r["setup_compile_s"]),
        "setup_parse_s": med(lambda r: r["setup_parse_s"]),
        "setup_compile_s": med(lambda r: r["setup_compile_s"]),
        "save_s": med(lambda r: r["save_s"]),
        "wall_s": statistics.median(wall),
        "steps_per_primary": med(lambda r: r["counters"]["steps_per_primary"]),
        "phases": {
            k: med(lambda r, k=k: r["phases"][k])
            for k in ("fill_s", "zone_ref_s", "distance_s", "step_s", "compact_s", "sum_s")
        },
        "counters": {
            k: med(lambda r, k=k: r["counters"][k])
            for k in ("steps", "nuclear_events", "secondaries", "neutrons_banked", "fragments_banked")
        },
    }
    agg["transport_fraction"] = agg["transport_s"] / agg["wall_s"] if agg["wall_s"] else 0.0
    return agg


# ---- Profiler wrappers -----------------------------------------------------------


def run_callgrind(binary: Path, case_dir: Path, nstat: int, out_dir: Path, sid: str) -> str:
    """callgrind run + annotated flat profile; returns a status string."""
    if not shutil.which("valgrind"):
        return "skipped (valgrind not found)"
    cg_out = out_dir / f"callgrind.{sid}.out"
    try:
        run_once(
            binary,
            case_dir,
            nstat,
            out_dir / "callgrind_run",
            extra_argv=[
                "valgrind",
                "--tool=callgrind",
                f"--callgrind-out-file={cg_out}",
            ],
        )
    except RuntimeError as exc:
        return f"failed: {exc}"
    annotate = shutil.which("callgrind_annotate")
    if not annotate:
        return f"raw output only at {cg_out} (callgrind_annotate not found)"
    proc = subprocess.run([annotate, "--threshold=99", str(cg_out)], capture_output=True, text=True)
    table = out_dir / f"callgrind.{sid}.txt"
    table.write_text(proc.stdout)
    return f"flat profile at {table}"


def run_perf(binary: Path, case_dir: Path, nstat: int, out_dir: Path, sid: str) -> str:
    """perf stat + perf report flat profile; returns a status string."""
    if not shutil.which("perf"):
        return "skipped (perf not found)"
    stat_file = out_dir / f"perf_stat.{sid}.txt"
    try:
        run_once(
            binary,
            case_dir,
            nstat,
            out_dir / "perf_stat_run",
            extra_argv=["perf", "stat", "-o", str(stat_file)],
        )
    except RuntimeError as exc:
        return f"failed: {exc}"
    perf_data = out_dir / f"perf.{sid}.data"
    try:
        run_once(
            binary,
            case_dir,
            nstat,
            out_dir / "perf_record_run",
            extra_argv=["perf", "record", "-o", str(perf_data), "--"],
        )
    except RuntimeError as exc:
        return f"stat only at {stat_file}; record failed: {exc}"
    proc = subprocess.run(
        ["perf", "report", "--stdio", "--percent-limit", "0.5", "-i", str(perf_data)],
        capture_output=True,
        text=True,
    )
    report = out_dir / f"perf_report.{sid}.txt"
    report.write_text(proc.stdout)
    return f"stat at {stat_file}, report at {report}"


# ---- Reporting -----------------------------------------------------------------


def print_summary(records: list) -> None:
    header = (
        f"{'scenario':<24} {'nstat':>7} {'setup[s]':>9} {'transp[s]':>10} {'±%':>5} "
        f"{'prim/s':>9} {'steps/prim':>10} {'phase split f/z/d/s/c [%]':>28}"
    )
    print()
    print(header)
    print("-" * len(header))
    for rec in records:
        if rec["status"] != "ok":
            print(f"{rec['id']:<24} {rec['status']}")
            continue
        agg = rec["median"]
        ph = agg["phases"]
        total = agg["transport_s"] or math.inf
        split = "/".join(f"{100.0 * ph[k] / total:.0f}" for k in ("fill_s", "zone_ref_s", "distance_s", "step_s", "compact_s"))
        print(
            f"{rec['id']:<24} {rec['nstat']:>7} {agg['setup_s']:>9.3f} "
            f"{agg['transport_s']:>10.3f} {agg['transport_spread_pct']:>5.1f} "
            f"{agg['prim_per_s']:>9.0f} {agg['steps_per_primary']:>10.1f} {split:>28}"
        )
    print()


def check_quality(rec: dict) -> list:
    """Methodology warnings for one finished scenario record."""
    warnings = []
    agg = rec["median"]
    if agg["transport_fraction"] < TRANSPORT_FRACTION_WARN:
        warnings.append(
            f"{rec['id']}: transport is only {100 * agg['transport_fraction']:.0f}% of wall "
            f"time (target >= {100 * TRANSPORT_FRACTION_WARN:.0f}%); raise nstat"
        )
    if agg["transport_s"] > 0:
        ratio = agg["phases"]["sum_s"] / agg["transport_s"]
        if abs(1.0 - ratio) > PHASE_SUM_TOLERANCE:
            warnings.append(
                f"{rec['id']}: phase sum is {100 * ratio:.1f}% of transport wall time "
                f"(expected within {100 * PHASE_SUM_TOLERANCE:.0f}%)"
            )
    return warnings


# ---- Main ----------------------------------------------------------------------


def select_scenarios(filter_expr: str) -> list:
    wanted = [tok.strip() for tok in filter_expr.split(",") if tok.strip()]
    if "all" in wanted:
        return list(SCENARIOS)
    selected = []
    for s in SCENARIOS:
        if s.id in wanted or any(tag in wanted for tag in s.tags):
            selected.append(s)
    return selected


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--filter",
        default="core",
        help="comma-separated scenario ids or tags (core, sweep, s1, s2, s4, gemca_stress, dicom, all); default: core",
    )
    parser.add_argument("--list", action="store_true", help="list scenarios and exit")
    parser.add_argument("--repeats", type=int, default=3, help="recorded repeats (default 3)")
    parser.add_argument("--no-warmup", action="store_true", help="skip the warm-up run")
    parser.add_argument("--nstat-scale", type=float, default=1.0, help="scale every scenario's nstat")
    parser.add_argument(
        "--preset",
        default="release",
        choices=sorted(PRESET_BUILD_DIRS),
        help="CMake preset whose binary is measured (default: release)",
    )
    parser.add_argument("--build", action="store_true", help="configure+build before running")
    parser.add_argument(
        "--allow-rebuild",
        action="store_true",
        help="permit scenarios that need their own build (S2 pool-capacity sweep)",
    )
    parser.add_argument("--bin", type=Path, default=None, help="explicit binary path override")
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO_ROOT / "tools/bench/results/results.json",
        help="aggregated results file (default tools/bench/results/results.json)",
    )
    parser.add_argument("--callgrind", action="store_true", help="also produce a callgrind flat profile per scenario")
    parser.add_argument("--perf", action="store_true", help="also produce perf stat + perf report per scenario")
    parser.add_argument(
        "--prof-nstat-divisor",
        type=int,
        default=10,
        help="divide nstat by this for callgrind/perf wrapper runs (default 10)",
    )
    args = parser.parse_args(argv)

    scenarios = select_scenarios(args.filter)
    if args.list:
        for s in SCENARIOS:
            mark = "*" if s in scenarios else " "
            print(f" {mark} {s.id:<24} nstat={s.nstat:<7} tags={','.join(s.tags):<20} {s.title}")
        return 0
    if not scenarios:
        print(f"error: filter '{args.filter}' selects no scenarios", file=sys.stderr)
        return 2

    if args.build and not args.bin:
        print(f"building preset '{args.preset}' ...")
        cmake_build(args.preset)

    default_binary = args.bin or binary_path(build_dir_for(args.preset))
    if not default_binary.exists():
        print(f"error: binary not found at {default_binary}; pass --build or --bin", file=sys.stderr)
        return 2

    results_dir = args.output.parent
    results_dir.mkdir(parents=True, exist_ok=True)
    records = []
    all_warnings = []

    with tempfile.TemporaryDirectory(prefix="oshbench_") as tmp:
        work_root = Path(tmp)
        for scenario in scenarios:
            nstat = max(1, int(scenario.nstat * args.nstat_scale))
            rec = {
                "id": scenario.id,
                "title": scenario.title,
                "case": scenario.case,
                "nstat": nstat,
                "pool_capacity": scenario.pool_capacity,
                "status": "ok",
            }
            print(f"[{scenario.id}] nstat={nstat} ...", flush=True)
            try:
                binary = default_binary
                if scenario.pool_capacity is not None:
                    if not args.allow_rebuild:
                        rec["status"] = "skipped (needs --allow-rebuild)"
                        records.append(rec)
                        continue
                    print(f"  rebuilding with pool capacity {scenario.pool_capacity} ...")
                    binary = binary_path(cmake_build(args.preset, scenario.pool_capacity))

                case_dir = prepare_case(scenario, work_root)
                run_dir = work_root / f"{scenario.id}_out"

                if not args.no_warmup:
                    run_once(binary, case_dir, nstat, run_dir / "warmup")
                runs = [run_once(binary, case_dir, nstat, run_dir / f"rep{i}") for i in range(args.repeats)]
                rec["runs"] = runs
                rec["median"] = aggregate(runs)
                all_warnings.extend(check_quality(rec))

                prof_nstat = max(1, nstat // args.prof_nstat_divisor)
                if args.callgrind:
                    status = run_callgrind(binary, case_dir, prof_nstat, results_dir, scenario.id)
                    rec["callgrind"] = status
                    print(f"  callgrind: {status}")
                if args.perf:
                    status = run_perf(binary, case_dir, prof_nstat, results_dir, scenario.id)
                    rec["perf"] = status
                    print(f"  perf: {status}")
            except (RuntimeError, subprocess.CalledProcessError) as exc:
                rec["status"] = f"error: {exc}"
            records.append(rec)

    document = {
        "schema": 1,
        "timestamp": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "preset": args.preset,
        "repeats": args.repeats,
        "nstat_scale": args.nstat_scale,
        "git": git_commit(),
        "machine": machine_metadata(),
        "results": records,
    }
    args.output.write_text(json.dumps(document, indent=2) + "\n")
    print(f"results written to {args.output}")

    print_summary(records)
    for warning in all_warnings:
        print(f"WARNING: {warning}")
    return 0 if all(r["status"] == "ok" or r["status"].startswith("skipped") for r in records) else 1


if __name__ == "__main__":
    sys.exit(main())
