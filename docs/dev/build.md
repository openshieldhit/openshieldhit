# Build system

openshieldhit uses CMake.  See also `DEVELOPER.md` in the repository root for the
full developer build guide including coverage, sanitizers, and packaging.

## Common build configurations

```bash
# release build (optimised, no debug symbols)
cmake -B build_rel -DCMAKE_BUILD_TYPE=Release .
cmake --build build_rel --parallel

# debug build
cmake -B build_dbg -DCMAKE_BUILD_TYPE=Debug .
cmake --build build_dbg --parallel

# coverage build
cmake -B build_cov -DCMAKE_BUILD_TYPE=Debug -DCOVERAGE=ON .
cmake --build build_cov --parallel
```

## Running tests

```bash
cd build && ctest --output-on-failure -j4
```

Tests are split into `unit::*` (fast, isolated) and `cases::*` (end-to-end,
slower).  The CI runs both; see `.github/workflows/test.yml`.

## Packaging

See `PACKAGING.md` for DEB / TGZ / ZIP package creation.
