# DICOM test fixtures

Complete datasets with CT, RS, RD and RN files go here.
Keep this lean, < 500 MB. Only add what is absolutely needed for test and development.

## Fetched on demand (not stored in git)

The `DCPT_headphantom` dataset (182 `.dcm` files, ~267 MB uncompressed) used by
the DICOM tests is **not** committed to the repository. It used to live in Git
LFS, but the project exceeded its LFS bandwidth budget (see
[#187](https://github.com/openshieldhit/openshieldhit/issues/187)), so the data
now ships as an asset on a GitHub release and is downloaded on demand.

- Release: [`test-fixtures-v1`](https://github.com/openshieldhit/openshieldhit/releases/tag/test-fixtures-v1)
- Asset: `DCPT_headphantom.tar.gz`
- SHA-256: `5ab7ac365ccbbf4e5d3eb915b30f8a123fce3d460b05c7a17ce4ce9d64707095`

Git LFS is therefore no longer required to clone or build the project.

### One-step fetch (Linux, macOS, Windows)

Run from the repository root:

```sh
cmake -P cmake/fetch_fixtures.cmake
```

This downloads the archive, verifies its SHA-256, and unpacks it into
`tests/fixtures/dicom/DCPT_headphantom/`. It uses only CMake built-ins
(`file(DOWNLOAD)` and `cmake -E tar`), so it behaves identically on every
platform — no `curl`, `wget`, `tar` or shell required.

The script is idempotent: a stamp file
(`DCPT_headphantom/.fixtures-test-fixtures-v1.stamp`) marks a completed fetch, so
re-running it is a no-op.

### Automatic fetch at configure time

When you configure the project with testing enabled (the default), CMake runs
the fetch script for you. To skip it (e.g. when building offline or when you
only run the non-DICOM tests):

```sh
cmake -S . -B build -DOSH_FETCH_DICOM_FIXTURES=OFF
```

CI fetches the fixtures the same way, during the configure step, before the
DICOM tests run.

## Provenance

> **TODO (tracked in [#187](https://github.com/openshieldhit/openshieldhit/issues/187)):**
> document the origin, anonymization status and redistribution license of the
> `DCPT_headphantom` dataset. It is a head-phantom CT with an associated proton
> RTPLAN/RTDOSE (Brain_fin2). Confirm patient-identifying DICOM tags are scrubbed
> and that the data may be redistributed in a public repository before treating
> this release as permanent.
