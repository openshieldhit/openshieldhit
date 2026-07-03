# Command line

`openshieldhit` runs one case directory at a time.

## Basic usage

```bash
build/bin/openshieldhit [OPTIONS] [WORKDIR]
```

`WORKDIR` is the case directory. If given, openshieldhit reads:

- `WORKDIR/geo.dat`
- `WORKDIR/beam.dat`
- `WORKDIR/mat.dat`
- `WORKDIR/detect.dat`

If `WORKDIR` is omitted, the current directory is used.

## Common examples

Run a case and write outputs into the case directory:

```bash
build/bin/openshieldhit tests/cases/00_minimal/
```

Run a case and write outputs into a separate directory:

```bash
build/bin/openshieldhit --outdir /tmp/osh-run tests/cases/00_minimal/
```

The directory passed with `--outdir` is created automatically if it does not
already exist.

Parse and compile inputs only:

```bash
build/bin/openshieldhit --dry-run tests/cases/00_minimal/
```

Override the number of primaries from the command line:

```bash
build/bin/openshieldhit --nstat 100000 tests/cases/00_minimal/
```

Run with verbose logging:

```bash
build/bin/openshieldhit --verbose tests/cases/00_minimal/
```

Run for at most a fixed wall-clock time, then stop cleanly and save the partial
result:

```bash
build/bin/openshieldhit --max-time 30m tests/cases/00_minimal/
```

Refine the output files on disk while the run progresses, writing the current
partial result every 10 minutes of wall time:

```bash
build/bin/openshieldhit --dump-every 10m tests/cases/00_minimal/
```

## Options

```text
-h, --help            Show this help message
-V, --version         Print version information
-v, --verbose         Increase verbosity
-n, --nstat <n>       Number of requested primary histories
-N, --seedoffset <n>  Random seed offset override (max 9999)
    --dry-run         Parse/load inputs only, do not run transport
    --workdir <dir>   Working directory for default input/output files
-g, --geo <file>      Override geometry input file
-b, --beam <file>     Override beam input file
-m, --mat <file>      Override material input file
-d, --detect <file>   Override scoring input file
-o, --outdir <dir>    Override output directory
    --profile <file>  Write a one-line JSON timing/counter profile to <file>
    --max-time <dur>  Wall-time budget; stop cleanly and save the partial result
                      (e.g. 30s, 30m, 1h, or a bare number of seconds)
    --dump-every <dur>            Periodically overwrite the output files with the
                                  current partial result, every <dur> of wall time
    --dump-every-primaries <n>    Same, but every <n> completed primaries
    --score-replicas <n>          Diagnostic: transport the run as <n> sequential
                                  private-accumulator sub-ranges, then merge (n <= nstat)
```

## Stopping a run early

A run stops cleanly — rather than being killed mid-history — in two ways:

- **Wall-time budget.** `--max-time <dur>` accepts a duration with an optional
  unit suffix: `s` (seconds, also the default), `m` (minutes), or `h` (hours).
  It overrides the `MAXTIME` card in `beam.dat`; `MAXTIME` accepts the same
  syntax. `0` (the default) means unlimited.
- **Ctrl-C (SIGINT).** Pressing Ctrl-C once requests the same clean stop. On
  Windows the console Ctrl-C / close handler does the same.

When a clean stop fires, transport stops *injecting* new primaries but lets
every in-flight history finish and drains all the secondary families they
spawned. The saved result is therefore **family-exact** for exactly the
primaries that completed: every output is normalised by that true completed
count (reported on the console and written to the `# PRIMARIES:` header), not by
the originally requested `nstat`. Because each history is an independent
function of its global index, the partial result is an unbiased Monte Carlo
estimate from the histories that finished.

## Periodic partial-result dumps

A long run can refine its output files on disk as it goes, so you can inspect
convergence without waiting for the end or stopping the run.

Throughout, a **cadence** is the interval at which the run pauses at a live
checkpoint to attempt a dump — expressed either by **wall time** (`--dump-every`)
or by **completed primary count** (`--dump-every-primaries`). Two cadences select
how often a dump fires; both overwrite the normal output files with the current
partial result and then let the run continue, untouched.

- **By wall time.** `--dump-every <dur>` dumps roughly every `<dur>` of wall-clock
  time (same duration grammar as `--max-time`: `30s`, `10m`, `1h`, or bare
  seconds). This is the cadence to use in production and on many cores: its
  overhead is bounded per wall-hour regardless of how fast the machine is. It
  overrides the `DUMPEVERY` card in `beam.dat`.
- **By primary count.** `--dump-every-primaries <n>` dumps every `<n>` completed
  primaries. This cadence is deterministic and reproducible (each history is a
  pure function of its global index), which makes it the right choice for tests.
  It overrides the `beam.dat` `NSTAT <n> <step>` save step (`nsave`).

A command-line flag always overrides the corresponding `beam.dat` card. With no
dump flag and no card, no periodic dumps are written (the fastest path).

**Every dump is physically exact.** A dump is taken only at a *family-complete
checkpoint* — a point where every secondary family (neutrons, fragments, …) the
completed primaries banked has been fully transported into scoring. A mid-run
ratio such as `DLET = Σ(LET·dose)/Σdose` is therefore unbiased, not merely
noisier. Each output records this with a `# COMPLETENESS: exact` header (ASCII)
or an `OSHBDO_RT_COMPLETENESS` token (BDO), alongside the `# PRIMARIES:` count it
was normalised by. Taking the snapshot never mutates the live accumulators, so
the run's final result is identical whether or not dumps were taken.

**On-demand dump (POSIX).** Sending `SIGUSR1` to the process
(`kill -USR1 <pid>`) requests a one-off dump at the next checkpoint. Because
checkpoints exist only while a cadence is running, `SIGUSR1` is meaningful **only
alongside** a `--dump-every[-primaries]` cadence: it then services the request at
the next checkpoint (ahead of the cadence's own schedule), so a cadence both
enables on-demand dumps and bounds how soon they land. With no cadence there are
no intermediate checkpoints and the signal has no effect. There is no `SIGUSR1`
on Windows, so on-demand dumps are a no-op there; the scheduled cadences work on
every platform.

**Memory.** When a dump cadence is set, the small extra buffer a snapshot needs
is reserved up front and shown in the `Scoring memory:` line, so a scheduled dump
can never run the process out of memory partway through. Runs without a dump
cadence pay nothing for this.

## Diagnostics

### `--score-replicas <n>` — sequential private-accumulator harness

`--score-replicas <n>` splits the run's `[0, nstat)` histories into `n` contiguous
sub-ranges, transports them **one after another** (no threads), each depositing
into its **own private accumulator set**, then **merges** all `n` into the master
before the normal postprocess and save. It is a **diagnostic and profiling
harness, not a speed-up**: it runs no faster than a plain run (often a touch
slower, from the extra accumulator memory and the final merge).

Its purpose is to exercise — with zero concurrency risk — the per-worker
private-accumulator and merge machinery that the future threaded, MPI, and WASM
Web Worker backends all depend on, and to make the parallel **reproducibility
contract** testable before any threading exists:

| Invocation | Guarantee |
|---|---|
| `--score-replicas 1` vs a plain run | **bit-identical** — one range → private set → merge into an empty master is the same summation order as serial |
| `--score-replicas N` (N>1) vs a plain run | **equal within tolerance** — identical per-history physics (each history's RNG stream is a pure function of its global index); only the cross-partition floating-point summation order differs |
| `--score-replicas N` run **twice** | **bit-identical to each other** — a fixed partition gives a deterministic merge order |
| `--score-replicas N` vs `--score-replicas M` (N≠M) | **not** guaranteed identical — a different partition changes the merge grouping |

`n` must be at least 1 and at most `nstat` (each replica needs at least one
history); `0` or `n > nstat` is rejected. With the flag absent, the run takes the
ordinary shared-master path, byte-for-byte unchanged.

## Notes

- `WORKDIR` and `--workdir` serve the same purpose: selecting the default
  location of `geo.dat`, `beam.dat`, `mat.dat`, and `detect.dat`.
- The file override options replace only one input file each; they do not
  affect the others.
- If `--outdir` is not given, outputs are written into the working directory.
