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
```

## Notes

- `WORKDIR` and `--workdir` serve the same purpose: selecting the default
  location of `geo.dat`, `beam.dat`, `mat.dat`, and `detect.dat`.
- The file override options replace only one input file each; they do not
  affect the others.
- If `--outdir` is not given, outputs are written into the working directory.
