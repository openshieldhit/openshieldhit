# Getting started

## Requirements

- CMake ≥ 3.16
- C11 compiler (GCC, Clang, or MSVC)
- Optional: zlib (compressed `.bdz` output files), libSDL2 (interactive geometry viewers in `examples/`, not needed for the main application)
- DICOM: no external library needed — openshieldhit includes a minimal self-contained reader/writer

## Build

```bash
git clone https://github.com/openshieldhit/openshieldhit.git
cd openshieldhit
cmake --preset release
cmake --build --preset release --parallel
```

The `openshieldhit` binary is written to `build/bin/`.

## Install system-wide

Install to the default prefix (`/usr/local/bin/openshieldhit`):

```bash
sudo cmake --install build
```

Install to `/usr` (`/usr/bin/openshieldhit`):

```bash
sudo cmake --install build --prefix /usr
```

Install into your home directory — no sudo needed (`~/.local/bin/openshieldhit`):

```bash
cmake --install build --prefix ~/.local
```

## Run a case

```bash
build/bin/openshieldhit path/to/case/
```

The case directory must contain `beam.dat`, `geo.dat`, `mat.dat`, and
`detect.dat`.  Output files are written into the same directory.

A few common variations:

Parse and compile inputs only, without running transport:

```bash
build/bin/openshieldhit --dry-run path/to/case/
```

Override the history count:

```bash
build/bin/openshieldhit -n 50000 path/to/case/
```

Enable verbose logging:

```bash
build/bin/openshieldhit -v path/to/case/
```

Write outputs somewhere else (the directory is created automatically):

```bash
build/bin/openshieldhit --outdir /tmp/osh-run path/to/case/
```

Run for at most a fixed wall-clock time, then stop cleanly and save the partial
result (you can also press Ctrl-C at any time for the same clean stop):

```bash
build/bin/openshieldhit --max-time 30m path/to/case/
```

The partial result is normalised by the number of primaries that actually
finished, so it stays a valid (just noisier) estimate — see
[stopping a run early](command-line.md#stopping-a-run-early).

## Controlling when a run ends

A run can end for two reasons, and **whichever limit is reached first wins**:

1. **Primary count** — the number of source histories, `nstat`. This is always
   set (a run needs a target), and the run stops as soon as that many primaries
   have been simulated. Set it with `NSTAT` in `beam.dat` or with `-n` / `--nstat`
   on the command line.
2. **Wall-time budget** — an optional real-time limit. Set it with `MAXTIME` in
   `beam.dat` or with `--max-time` on the command line; `0` (the default) means
   *no* time limit. When the budget elapses the run stops **cleanly**: in-flight
   histories finish, all secondaries drain, and the partial result is saved,
   normalised by the primaries that actually completed.

So with both set, the run stops at the **first** of "all `nstat` primaries done"
or "time budget elapsed". Pressing **Ctrl-C** requests the same clean stop at any
moment. (Getting a *preview while the run keeps going* is a different feature —
see [periodic dumps](command-line.md#periodic-partial-result-dumps).)

### Precedence: command line overrides beam.dat

Each setting can be given in `beam.dat` **and** on the command line. The rule is
simple: **a command-line flag always overrides the matching `beam.dat` card** for
that one setting. The two files are not "merged" beyond this per-setting override.

| Setting | `beam.dat` card | Command-line flag (wins) |
|---|---|---|
| Primary count | `NSTAT <n>` | `-n <n>` / `--nstat <n>` |
| Wall-time budget | `MAXTIME <dur>` | `--max-time <dur>` |
| Periodic dump, by time | `DUMPEVERY <dur>` | `--dump-every <dur>` |
| Periodic dump, by primaries | `NSTAT <n> <step>` (the `<step>` field) | `--dump-every-primaries <n>` |

Durations accept `s` / `m` / `h` suffixes (or a bare number of seconds), e.g.
`30s`, `10m`, `1h`.

### Examples

Stop after one million primaries (no time limit):

```bash
build/bin/openshieldhit -n 1000000 path/to/case/
```

Run up to one million primaries **or** two hours, whichever comes first:

```bash
build/bin/openshieldhit -n 1000000 --max-time 2h path/to/case/
```

`beam.dat` says `NSTAT 1000000` but you want a quick smoke run — the flag wins,
so this simulates 50 000:

```bash
build/bin/openshieldhit -n 50000 path/to/case/
```

`beam.dat` says `MAXTIME 1h`, but tonight you can spare only 20 minutes — the flag
wins:

```bash
build/bin/openshieldhit --max-time 20m path/to/case/
```

Long run that also refines its output files every 10 minutes **without stopping**
(a partial-result dump, not a stop), capped at a 6-hour budget:

```bash
build/bin/openshieldhit -n 100000000 --dump-every 10m --max-time 6h path/to/case/
```

## Try the minimal example

```bash
build/bin/openshieldhit tests/cases/00_minimal/
```

Produces a Bragg-peak depth-dose profile for 200 MeV protons in water.

## Next steps

- [command-line reference](command-line.md) — CLI options, output directories, dry-run, file overrides
- [beam.dat reference](beam.dat.md) — particle, energy, beam geometry, physics switches
- [geo.dat reference](geo.dat.md) — geometry bodies and zones
- [detect.dat reference](detect.dat.md) — scoring detectors
