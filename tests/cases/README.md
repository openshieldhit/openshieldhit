# Integration Cases

This directory contains self-contained OpenShieldHIT case directories used for
end-to-end and regression testing.

Each case is expected to live in its own subdirectory, for example:

- `01_simple_detect/`
- `02_sobp/`

A case directory may contain:

- `geo.dat`, `beam.dat`, `mat.dat`, `detect.dat`
- optional `README`
- optional `args.cmake` to override default test runner arguments
- optional `expected/` with reference `stdout`, `stderr`, exit code, or `*.dat` output files
- reference output files such as `ref_*.bdo`

These directories should be runnable as realistic user-style inputs, not just
minimal parser snippets.

Some heavier cases intentionally keep the default `--dry-run` CTest behavior and
are meant for manual transport runs plus plotting against external reference
fixtures rather than automated output comparison.
