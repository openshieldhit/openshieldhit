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

## Notes

- `WORKDIR` and `--workdir` serve the same purpose: selecting the default
  location of `geo.dat`, `beam.dat`, `mat.dat`, and `detect.dat`.
- The file override options replace only one input file each; they do not
  affect the others.
- If `--outdir` is not given, outputs are written into the working directory.
