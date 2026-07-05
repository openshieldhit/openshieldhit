# OpenShieldHIT bug-hunting report (2026-07-05)

Deep audit of `main` (e3c619a) prompted by PR #239 (sanitizer wiring) and issues
#235/#236/#237, with special attention to physics correctness, numerical
stability, parallelization readiness (#161, #165–#169), stddev/variance plans
(#169/#195), and architecture (TODO.md). Findings are numbered per subsystem
and tagged with **severity** (critical / high / medium / low) and **confidence**
(high = verified by running code or by construction; medium = strong code
reading; low = suspicion worth checking). Line numbers refer to `main` at
e3c619a unless a branch is named.

Method: manual code audit of the hot path outward, plus hands-on experiments:
debug + ASan/UBSan instrumented builds, full `ctest` runs, targeted standalone
harnesses. Every finding cites `file:line` evidence.

Status: **in progress** — sections are appended step by step, one commit each.

---

## 0. Baseline, test infrastructure, and working-tree observations

### T-1 (high, high) Integration-case registration turns any stray directory into a phantom failing test

`tests/cases/CMakeLists.txt:20-36` registers a CTest case for **every**
subdirectory via `file(GLOB ... "*/")` with no check that the case's input
files exist. On this machine, `cases::06_minimal_nucre` fails on a clean
`main` build: the case directory is tracked only on the local
`feat/parallel-threads` branch, but after switching back to `main` the
directory survives because it still holds *ignored* leftovers
(`NB_msh.bdo`, `NB_msh.dat` — matched by `.gitignore:14-15`). The glob then
registers a test whose `beam.dat` does not exist → exit code 66, red suite.

- Impact: any developer who ever checks out a branch adding a case (or whose
  tools drop a directory under `tests/cases/`) gets a permanently failing
  suite on `main`; erodes trust in the suite ("red is normal").
- Fix direction: in the glob loop, `continue()` unless `${case_dir}/beam.dat`
  (and geo/mat/detect) exist — or better, drive registration from a tracked
  manifest. Same pattern exists in `tests/reference/CMakeLists.txt:19`.

### T-2 (high, high) The entire physics-validation tier (reference::idd_*) never runs — no reference curves are committed

All 13 `reference::idd_*` tests report **Skipped**: the runner
(`tests/reference/run_idd.cmake:41`) skips when `${CASE_DIR}/reference/` has
no data, and *no* `idd_*/reference/` directory exists in the tree (checked all
of `tests/reference/idd_*/`). So the IDD comparisons against SH12A/FLUKA that
TODO.md cites as validation evidence ("distal 80–20% falloff 0.4522 vs
0.4533 cm", STRAGG 2 section) are not locked in by CI at all — a physics
regression in stopping power, straggling, or scattering would merge silently.

- Impact: no automated guard on physics output; the strongest claims in
  TODO.md are unenforced.
- Fix direction: commit the reference curves (they are small 1-D IDD tables)
  for at least 70/150/200 MeV water, and add a CI job that runs
  `ctest -L reference` (possibly nightly if runtime is a concern; cases have
  TIMEOUT 1800 already).

### T-3 (medium, high) Tests write outputs into the current working directory (source-root pollution)

Confirmed stray files at repo root produced by test binaries run from the
source root: `osh_transport_parallel_test_{0..9}.tmp`,
`osh_transport_parallel_out_{serial_0,2t_2,4t_4}.dat`, `serial_out.dat`
(dated 2026-07-03, version `v0.0.8-7-g5725239-dirty`).

- `feat/parallel-threads:tests/unit/test_osh_transport_parallel.c:27,138,152,166`
  builds all scratch/output filenames as bare relative paths; the two smoke
  tests never `remove()` their outputs, and `test_serial_matches_parallel`
  removes them only on the success path.
- `main:tests/unit/test_osh_run_dump.c:72,121,149` similarly creates
  `osh_rundump_out_<n>` directories relative to CWD.
- Impact: `ctest` executed from the source root (or any manual run) litters
  the repo; leftovers then interact badly with T-1.
- Fix direction: root all test scratch under `OSH_TEST_TMPDIR` supplied by
  CMake (`CMAKE_CURRENT_BINARY_DIR`), and always clean up in a teardown path
  (or register with `FIXTURES_CLEANUP`).

### T-4 (info, high) PR #239 validated hands-on

Reproduced the PR's claims independently on `main` + manual
`-fsanitize=address,undefined -fno-sanitize-recover=all` build (GCC, Debug):
the **only** failures in the whole suite are exactly the two the PR fixes —
`unit::test_osh_scoring_parse_geometry` (stack-use-after-scope) and
`unit::test_osh_scoring_step` (stack OOB read via `nenergy = 0` table), both
aborting under ASan. Everything else is green with `detect_leaks=0`.
Conclusion: PR #239's test fixes and CI job are correct and worth merging
as-is; the sanitizer job would have caught both bugs at introduction time.

### T-5 (low, high) Working-tree / repo hygiene observations

- A file literally named `out_%%d.dat` sits at the repo root (written
  2026-07-03 by a run-control-branch binary): a `Filename` containing a
  printf-style `%d` pattern is used verbatim by the savers
  (`src/scoring/save/osh_scoring_save_ascii.c:70` opens `out->filename`
  as-is). If numbered periodic-dump filenames are ever intended (#193), there
  is no templating support; if not intended, nothing warns the user that `%d`
  is meaningless. Low priority, but worth a decision.
- `.venv/**` (a full Python virtualenv, thousands of files) is currently
  *staged* in the local index on `main`. Not a repo bug, but one `git commit
  -a` away from becoming one; `.gitignore` should cover `.venv/`.

---

*(Sections below are appended as the audit proceeds.)*
