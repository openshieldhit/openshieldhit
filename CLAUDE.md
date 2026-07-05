# CLAUDE.md

Operating brief for AI coding agents (Claude Code and others) working **inside**
this repository. It is intentionally short: it points at the authoritative docs
and captures the handful of repo-specific rules and commands that trip agents up.
It is **not** a second copy of those docs — when this file and a linked document
disagree, the linked document wins.

For *what the project is* and where things live, read [`llms.txt`](llms.txt) first.

## Read before editing

| Document | What it is the source of truth for |
|---|---|
| [`DEVELOPER.md`](DEVELOPER.md) | Coding style, struct layout, C11 rules, module/dependency layering, presets — **the** style authority. Organised as numbered rules (§1–§15); cite them (e.g. "violates §1"), don't paraphrase |
| [`llms.txt`](llms.txt) | Project overview, source-tree map, unit conventions, SH12A/DICOM notes |
| [`CONTRIBUTING.md`](CONTRIBUTING.md) | Build prerequisites, formatting, how to submit a patch |
| [`docs/dev/`](docs/dev/) | Architecture, testing, scoring, coordinates, RNG deep-dives |
| [`docs/user/`](docs/user/) | Input-file (`beam.dat` / `geo.dat` / `mat.dat` / `detect.dat`) and CLI reference |
| `src/<module>/README.md` | What a given module owns and its `parse/` vs `runtime/` split |

## Build, test, format

Day-to-day work uses the **debug** preset — a plain `-Og` build (no sanitizers).
For an AddressSanitizer + UBSan run, use the **asan** preset (`build_asan/`),
which is what the `sanitizers` CI job builds and tests.
One command per block:

```bash
cmake --preset debug
```
```bash
cmake --build --preset debug --parallel
```
```bash
ctest --test-dir build_debug --output-on-failure
```

> **There are no CTest presets.** `ctest --preset debug` fails — CMakePresets.json
> defines only configure/build presets. Always point CTest at the binary dir:
> `ctest --test-dir build_debug` (add `-C Debug` on Windows / multi-config).

Optimised build for benchmarking (binaries land in `build/bin/`, debug in `build_debug/bin/`):

```bash
cmake --preset release
```
```bash
cmake --build --preset release --parallel
```

Format C/H before committing — clang-format (18+) and clang-tidy are enforced in CI:

```bash
./tools/clang-format-all.sh
```

Python helpers under `tools/` are formatted with `ruff` (see `.pre-commit-config.yaml`).

## Hard rules (agents break these most)

These are the repo-specific traps, each tagged with the DEVELOPER.md section that
owns it. The § is the authority — read it before relying on the one-liner here.

- **C11**, but declare all variables at the **top of the block** (K&R); no declarations inside `for`/`if`. — §1
- **No `typedef struct`** — write `struct foo` explicitly everywhere. — §2.1
- **No `//` comments** in production code — block `/* */` only; `//` is reserved for temporary/WIP notes. — §4.1
- `double const *p`, **not** `const double *p`. — §2.2
- **No non-portable APIs** — Windows is a target. Avoid POSIX-only APIs (`strcasecmp`/`<strings.h>`, `mkdtemp`, `getpid`, `<unistd.h>`, `<sys/stat.h>`) and also `<threads.h>` (C11 threads; not implemented in MSVC). See the banned-API table, §11.1.
- **No heap allocation on the hot path** — nothing under `osh_scoring_score_step()` / `osh_transport_step()` may `malloc`/`calloc`/`realloc`/`free`. Pre-allocate at setup; scratch/counters are caller-owned (per-worker). — §10
- **Return conventions don't mix:** predicates return `int` (1 = yes); operations return `enum osh_status` (`OSH_OK = 0` = success). Do not confuse the two. — §5
- Naming: public functions use the `osh_` prefix; file-local helpers are `static` and unprefixed; `static inline` helpers use a leading `_`. — §6.1
- Case-insensitive keywords: lowercase once at parse time with `osh_lower_inplace()`, then `strcmp` — never `strcasecmp`. — §11.2

## Where things go

- Module layout: `src/<module>/` (API + domain types), optional `src/<module>/parse/` (input parsing) and `src/<module>/runtime/` (simulation-ready state). Add these layers only when the module actually has those phases. `runtime/` must not depend on `transport/`. Details in DEVELOPER.md §9.
- Internal headers sit next to their `.c`; `include/openshieldhit/` is public headers only.
- **Adding a unit test:** drop `tests/unit/test_<thing>.c` — it is auto-discovered by glob and registered as `unit::<file>`. Tests can read fixtures from `tests/fixtures/` via the `OSH_TEST_FIXTURES_DIR` macro (see `tests/unit/CMakeLists.txt` / `tests/unit/README.md`).

## Provenance — important

OpenShieldHIT shares **zero source code** with SHIELD-HIT or SHIELD-HIT12A (SH12A).
Do not call it a "port", "reimplementation", or "derivative" of them — those imply
code derivation and are incorrect. The shared ground is the input-file format,
application domain, and conceptual inspiration only. See the SH12A note in `llms.txt`.

## Submitting work

- One logical change per commit; PR description explains **why**, not just what (see CONTRIBUTING.md).
- Match the existing commit/PR title style: a module prefix helps, e.g. `scoring: …`, `transport: …`, `docs: …`.
- CI (Ubuntu/macOS/Windows) must be green: build, `ctest`, clang-format, clang-tidy.
