# Testing

## Test categories

| Label | Location | What it covers |
|-------|----------|----------------|
| `unit::*` | `tests/unit/` | Isolated unit tests — single functions or modules, fast |
| `cases::*` | `tests/cases/` | End-to-end cases — full parse + transport + score, slower |

Run the full suite:

```bash
cd build && ctest --output-on-failure -j4
```

Run only the unit tests:

```bash
cd build && ctest -L unit --output-on-failure
```

Run only the end-to-end cases:

```bash
cd build && ctest -L cases --output-on-failure
```

## Adding a unit test

Create `tests/unit/test_<module>.c`, add it to `tests/unit/CMakeLists.txt`
following the existing pattern, and implement a `main()` that returns 0 on
success and non-zero on failure.

## Adding a case

Create a subdirectory under `tests/cases/` with the four input files and an
expected-output reference.  Register it in `tests/cases/CMakeLists.txt`.
See `tests/cases/README.md` for the expected-output format.
