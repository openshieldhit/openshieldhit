# dicominspect

`dicominspect` is a tiny command-line inspection tool built on top of the
internal DICOM module.

## Purpose

It exists mainly as:

- a manual smoke-test for the DICOM reader
- a simple developer utility when checking CT or RTDOSE metadata
- a small example of how to wire `osh_diag_sink` into a standalone executable

It is intentionally minimal. It prints a short summary of:

- CT series geometry and HU scaling
- RTDOSE geometry, frame offsets, and peak dose

## Usage

```sh
dicominspect ct <directory>
dicominspect rtdose <file>
```

## Design notes

- Diagnostics go through a local sink callback.
- `INFO`/`DEBUG` go to `stdout`.
- `WARN` and above go to `stderr`.
- The tool inspects data only; it is not a general converter or validator.
