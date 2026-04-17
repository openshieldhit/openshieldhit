# CLI Parser (`src/cli`)

This module implements the internal command-line option parser used by the
`openshieldhit` executable.

## Why this exists

It is a deliberate replacement for `getopt/getopt_long` to keep behavior
portable across compilers and platforms, including MSVC where `<getopt.h>` is
not part of the standard C runtime.

## Scope

- Parses CLI flags into `struct osh_cli_options`
- Provides `osh_cli_print_help()`
- Does not parse input files (that lives in `src/apps/osh`)
- Is not part of the public installed API (`include/openshieldhit/...`)

## Main entry points

- `osh_cli_parse(...)`
- `osh_cli_print_help(...)`
