# Testing

## Test categories

| Label | Location | What it covers |
|-------|----------|----------------|
| `unit::*` | `tests/unit/` | Isolated unit tests — single functions or modules, fast |
| `cases::*` | `tests/cases/` | End-to-end cases — full parse + transport + score, slower |
| `dicom` | `tests/unit/`, `tests/cases/` | Tests that need the downloaded DICOM dataset (see below) |

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

## DICOM test fixtures (downloaded on demand)

A few tests exercise the DICOM reader against a real CT dataset
(`DCPT_headphantom`, 182 `.dcm` files). That data is **not** kept in the git
repository — it is too large and used to exhaust the Git LFS bandwidth budget.
Instead it is published as a [GitHub release asset][release] and downloaded
automatically the first time you run a test that needs it:

- `unit::test_osh_dicom`
- `unit::test_osh_geometry_dcm`
- `cases::05_dicom_simple`

These three carry the `dicom` label.

### When the download happens

The fetch is wired into CTest as a setup fixture (`DICOM_FIXTURES`), so it runs
**only at test time, on demand** — never while configuring or building:

| Action | Downloads the fixtures? |
|--------|-------------------------|
| `cmake -S . -B build` (configure) | No |
| `cmake --build build` (build the binary to run it) | No |
| `ctest -LE dicom` (everything except the DICOM tests) | No |
| Running any of the three DICOM tests above | Yes — once, beforehand |

The archive is **~26 MB compressed (~267 MB unpacked)** on disk. The download
runs once: a stamp file records the completed fetch, so later `ctest` runs
re-check the stamp and skip it. Run only the DICOM tests (fetching if needed)
with:

```bash
cd build && ctest -L dicom --output-on-failure
```

### How it is fetched (and why it works on Windows)

The download is done by a small CMake script, [`cmake/fetch_fixtures.cmake`][script],
using only CMake built-ins — `file(DOWNLOAD)` for the (SHA-256-verified)
transfer and `cmake -E tar` for extraction. It needs **no `curl`, `wget`,
`tar` or shell**, so it behaves identically on Linux, macOS and Windows. This
matters on Windows in particular, where tools like `wget` are usually not
available; because the logic lives inside CMake (which every build already
uses) there is nothing extra to install. `cmake -E tar` auto-detects the
archive format, so it handles `.tar.gz` and `.zip` assets alike.

You can also fetch the dataset manually at any time — for example to pre-warm
an offline checkout — by running, from the repository root:

```bash
cmake -P cmake/fetch_fixtures.cmake
```

See [`tests/fixtures/dicom/README.md`][fixtures-readme] for the dataset
provenance notes and the exact release coordinates.

[release]: https://github.com/openshieldhit/openshieldhit/releases/tag/test-fixtures-v1
[script]: https://github.com/openshieldhit/openshieldhit/blob/main/cmake/fetch_fixtures.cmake
[fixtures-readme]: https://github.com/openshieldhit/openshieldhit/blob/main/tests/fixtures/dicom/README.md

## Adding a unit test

Create `tests/unit/test_<module>.c`, add it to `tests/unit/CMakeLists.txt`
following the existing pattern, and implement a `main()` that returns 0 on
success and non-zero on failure.

## Adding a case

Create a subdirectory under `tests/cases/` with the four input files and an
expected-output reference.  Register it in `tests/cases/CMakeLists.txt`.
See `tests/cases/README.md` for the expected-output format.
