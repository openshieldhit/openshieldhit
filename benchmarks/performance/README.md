# OpenShieldHIT benchmark harness

Reusable profiling/benchmark infrastructure for the parallelization work
(issue #138, Phase 0 of `docs/dev/parallelization-plan.md`). One command runs
the scenario matrix, emits machine-readable JSON plus a human-readable summary
table, and (optionally) callgrind / perf flat profiles per scenario.

```
benchmarks/performance/
├── run_bench.py          # the harness: build, run matrix, aggregate, summarize
├── scenarios.py          # scenario registry (C1–C8 core cases, S1–S4 sweeps)
├── compare.py            # diff two result files; optional CI regression gate
├── gen_gemca_stress.py   # C8 generator: parameterized sphere-lattice geometry
├── cases/                # self-contained benchmark case directories
├── profiles/             # committed callgrind flat tables (deterministic)
└── baseline.json         # frozen pre-parallelization reference, committed
                          # once measured on dedicated hardware (see below)
```

## Quick start

```bash
# build the release preset and run the core scenarios (C1–C8), 3 repeats each
python3 benchmarks/performance/run_bench.py --build

# list scenarios; * marks what the current --filter selects
python3 benchmarks/performance/run_bench.py --list

# everything including sweeps; S2 needs one rebuild per pool capacity
python3 benchmarks/performance/run_bench.py --build --filter all --allow-rebuild

# quick sanity pass (scaled-down statistics, not for committed numbers)
python3 benchmarks/performance/run_bench.py --filter core --nstat-scale 0.05

# compare a run against the frozen baseline; fail on >5% transport regression
python3 benchmarks/performance/compare.py benchmarks/performance/baseline.json \
    benchmarks/performance/results/results.json --threshold 5
```

Results land in `benchmarks/performance/results/` (gitignored) unless `--output` says
otherwise.

## How a measurement works

Each scenario is one `openshieldhit` run with `--profile <file>`: the binary
itself reports setup (parse + compile), transport, and save wall times, the
five transport-phase timers (pool fill / zone-ref batch / distance batch /
step loop / compact), and event counters as a one-line JSON record. The
harness wraps that record with scenario, git-commit, and machine metadata
(CPU model, ISA flags, cores, RAM, OS, compiler) and aggregates repeats as
median + spread.

The phase timers live in the wavefront loop (`src/transport/osh_transport_ion.c`)
behind a NULL-pointer gate: when `--profile` is absent the only cost is one
pointer test per phase, and profiled runs are bit-identical to unprofiled
ones (the timers read the monotonic clock only — never the RNG streams).

## Scenario matrix

Core scenarios (tag `core`, the default filter):

| ID | What it stresses |
|----|------------------|
| `c1_p100_dose` | clean CSDA reference: analytic geometry, minimal scoring |
| `c2_p100_dose_fluence` | voxel-crossing raytrace, multi-page scoring |
| `c3_p100_let` | LET table lookups, dose/track-averaged accumulation |
| `c4_p150_nucre_filters` | per-step filter rules, secondary pool, Tripathi σ |
| `c5_ct_voxel` | CT voxel traversal + HU→ρ LUT (needs Git LFS fixtures) |
| `c6_c12_let` | heavy-ion physics (`z_eff`), LET scoring on ions |
| `c7_p70` / `c7_p200` | steps-per-primary scaling, setup-vs-transport ratio |
| `c8_z10…z5000` | GEMCA `eval_distance`/zone-lookup scaling with zone count |

Parameter sweeps:

| Tag | Knob |
|-----|------|
| `s1` | MSCAT / STRAGG / NUCRE individually off on the C4 base |
| `s2` | `OSH_TRANSPORT_POOL_CAPACITY` ∈ {1, 256, 4096, 65536} (rebuild per point, `--allow-rebuild`) |
| `s3` | zone-count scaling — this is the `c8_*` family, tag `gemca_stress` |
| `s4` | 0 / 3 / 10 filtered pages on one detector |

The C8 geometries are generated on the fly by `gen_gemca_stress.py`: a cubic
lattice of PMMA spheres in a water cylinder, where the water matrix zone's
boolean expression carries one negative term per sphere. Generate one by hand
with:

```bash
python3 benchmarks/performance/gen_gemca_stress.py --zones 1000 --out /tmp/c8
```

New scoring or physics features should land together with a scenario here
(see the "Extended scenarios" table in issue #138).

## Profiler wrappers

`--callgrind` runs each scenario (at `nstat / --prof-nstat-divisor`) under
callgrind and writes an annotated flat table next to the results file —
deterministic instruction counts, the right tool for regression comparisons.
`--perf` records `perf stat` (IPC, cache/branch misses — classifies hotspots
as compute- vs memory-bound) and a `perf report --stdio` flat profile. Both
skip gracefully when the tool is not installed. Use `--preset prof --build`
for frame-pointer-friendly binaries.

## Methodology rules (binding, from issue #138)

1. Measure only `release`/`prof` builds (the harness refuses anything else).
2. Transport must dominate: the harness warns when transport < 95% of wall
   time — raise `nstat` for that scenario rather than trusting the number.
3. callgrind for regression comparisons; perf for compute- vs memory-bound
   classification.
4. Attribute libm transcendental cost to the calling physics modules before
   drawing SIMD conclusions.
5. Numbers without their machine are meaningless — the results file embeds
   machine metadata; never hand-edit it.
6. Wall-clock hygiene is on you: mains power, `performance` governor, no
   other load. The harness does a warm-up run and ≥3 repeats and reports the
   spread; treat spreads above a few percent as a noisy machine.

## Freezing a baseline

Run on a quiet, mains-powered machine at a tagged commit:

```bash
python3 benchmarks/performance/run_bench.py --build --filter all --allow-rebuild \
    --output benchmarks/performance/baseline.json
git add benchmarks/performance/baseline.json && git commit
```

`compare.py` renders per-scenario transport-time and throughput deltas plus
phase-share shifts between any two result files, and exits non-zero with
`--threshold` on regressions, so it can be wired into CI later if wanted.
