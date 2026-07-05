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

## 1. Transport hot path & electromagnetic physics

Scope: `src/transport/osh_transport_ion_step.c`, `osh_transport_ion.c`,
`src/physics/atomic/*` (Bethe, Highland, Bohr/Vavilov/Landau straggling),
`src/material/runtime/*` (table build + lookup), `src/random/osh_rng.c`
(the samplers the hot path consumes).

### P-1 (high, high) Residual kinetic energy is silently deleted when a particle is killed at the energy cutoff

Every ion ends its life at `ion_step_setup()`:

- `osh_transport_ion_step.c:539-543` — `if (ctx->e0 <= ctx->cutoff) { pool->e[slot] = 0.0; ctx->done = 1; return; }` — no scoring call of any kind.
- The exit-energy clamps at `osh_transport_ion_step.c:897-898` and `:934-935` park a stopping particle at exactly `cutoff`, so this kill fires for essentially **every** track that ranges out (rather than escaping).
- `cutoff = max(OSH_MATERIAL_RUNTIME_EMIN, params->tcut, material_rt->emin) × A` (`cutoff_total_energy`, `osh_transport_ion_step.c:1156-1177`; `EMIN = 0.025` MeV/u, `osh_material_runtime.h:88`).

The wavefront loop (`osh_transport_ion.c:294-295`) merely compacts dead
entries; nothing deposits the residual kinetic energy. Contrast with the care
taken for **sub-threshold nuclear fragments**, which do get a point deposit
(`ion_point_deposit`, `osh_transport_ion_step.c:456-458`, issue #179) — the
cutoff kill is the same physical situation without the deposit.

- Impact, default settings: 0.025 MeV/u × A per stopping track (protons:
  25 keV ≈ 0.02% of a 150 MeV beam) — below statistical noise, but a
  *systematic* energy-balance deficit concentrated exactly at track ends
  (Bragg peak), visible in any energy-conservation audit.
- Impact, raised TCUT0: scales linearly. `TCUT0 = 2 MeV` deletes ~1.3% of a
  150 MeV proton beam's energy, all of it from the distal dose. SH12A/FLUKA
  deposit the residual energy locally at cutoff kill.
- Related silent kills on the same pattern: unknown species
  (`osh_transport_ion_step.c:532-535`, kill with no diagnostic at all —
  deserves at least a once-per-run WARN), and the defensive
  `residual_range <= 0` kills (`:772-775`, `:880-884`; provably unreachable
  given `preclip_step_len ≤ ds_csda`, but if ever reached they would drop the
  *entire* remaining energy — worth a diag).
- Secondaries born *between* `EMIN` and the effective cutoff (possible when
  `TCUT0 > EMIN`: the fragment-injection threshold at
  `osh_transport_ion_step.c:456` checks `EMIN`, not the run cutoff) are
  injected, then killed on their first step with the full energy dropped —
  doubly inconsistent.
- Fix direction: at the cutoff kill, point-deposit `e0` at the particle
  position (reuse `ion_point_deposit`); align the fragment-injection threshold
  with the effective cutoff.

### P-2 (high, high) Non-default isotopes are transported with the representative isotope's range table and rest mass — deuteron range wrong by ≈2×

Mechanism, three coupled pieces:

1. Runtime tables have **one projectile column per Z**, built for the
   *default isotope*: `osh_material_compile.c:698-716` (`Z=1→A=1`, `Z=2→A=4`,
   …), and the CSDA range integral bakes that A in:
   `integrate_range(..., proj.a, ...)` at `osh_material_compile.c:161-186`
   (`R(E/u) = ∫ A·d(E/u)/SP` — range is proportional to A at fixed E/u).
2. Transport maps **any isotope** of charge Z onto that column:
   `find_projectile_index`, `osh_transport_ion_step.c:1258-1261` explicitly
   ignores an A mismatch ("differing isotopes share the representative
   projectile for that Z for now") — with no runtime warning.
3. The per-step kinematics then mixes the *particle's* A with the *column's*
   mass: `e0/a_proj` uses the true A for the E/u lookup (good), but
   `ctx->proj_mass_mev = material_rt->projectile_mass_mev[projectile_idx]`
   (`osh_transport_ion_step.c:586`) is the representative isotope's mass.

Consequences (mass stopping power at fixed E/u is isotope-independent, range
is not: `R ∝ A/z²`):

- **Deuterons** (Z=1) use the proton range column: true `R_d(E/u) =
  2·R_p(E/u)`, so a deuteron's sampled energy loss per cm is 2× too high and
  its range is **half** the physical value. **Tritons**: 3×. **³He** (Z=2
  column is ⁴He): range overestimated by 4/3.
- β²/γ in the straggling parameters (`ξ`, `E_max`,
  `osh_transport_ion_step.c:911-914`) and in Highland `pv`
  (`osh_physics_scat_highland.c:28-33` via `proj_mass_mev`) are computed from
  total kinetic energy with the wrong rest mass — a deuteron is treated as a
  proton with twice the kinetic energy (γ−1 doubled).
- Who is affected: d/t/³He *secondaries* from Fermi break-up (small dose
  fraction, but their spatial distribution is systematically wrong), and —
  much worse — any **non-default-isotope primary beam** (deuteron or ³He
  beams are legitimate use cases; the run proceeds with no warning).
- Fix direction: within a Z column, rescale range lookups by
  `A_actual/A_rep` (exact up to the tiny Tmax term in the Bethe log) and use
  `part->mass` for kinematics — both cheap; or add isotope columns for Z=1,2.
  At minimum, WARN once per (Z,A) when the mismatch branch at
  `osh_transport_ion_step.c:1258` is taken.

### P-3 (medium, high) Bohr Gaussian straggling omits the relativistic variance factor — σ ~11% low at 200 MeV; discontinuity at the STRAGG 2 κ=10 dispatch

`osh_physics_strag_sigma()` (`osh_physics_strag_gauss.c:25-37`) implements the
**classical** Bohr variance `σ² = 0.1569·z_eff²·(Z/A)·d`, with no β
dependence. The relativistic form (PDG "Passage of particles", thick-absorber
Gaussian limit; equivalently `σ² = ξ·E_max·(1−β²/2)`) multiplies this by
`(1−β²/2)/(1−β²) = γ²(1−β²/2)`:

| proton energy | γ²(1−β²/2) | σ error |
|---|---|---|
| 70 MeV | 1.078 | −3.8% |
| 200 MeV | 1.236 | −10.5% |
| 1 GeV | ~1.9 | −27% |

Consequences:

- **STRAGG 1** (pure Gaussian) understates straggling width at therapeutic
  and especially at high energies → distal falloff too sharp.
- **STRAGG 2** dispatches Gaussian for κ ≥ 10 (`osh_transport_ion_step.c:916-918`).
  The Vavilov branch's sampled distribution carries the correct
  `ξ·E_max·(1−β²/2)` variance, so at the κ=10 boundary the model's variance
  jumps by the factor above — a discontinuity in the physics as a function of
  step thickness/energy. (In practice the Gaussian regime is entered at low
  energy where γ→1, which is why the SH12A distal-falloff validation in
  TODO.md didn't see it — but the branch is also taken for thick steps at
  high energy, e.g. coarse DELTAE in thick slabs.)
- Fix: multiply the variance by `γ²(1−β²/2)`; β², γ are already computed at
  the STRAGG 2 call site and cheap to add at the STRAGG 1 site. The comment
  in `osh_physics_strag_gauss.c:17-19` ("kept as a literal so the Gaussian
  result is unchanged") suggests the omission is legacy-compat, but nothing
  documents it as a deliberate physics approximation.

### P-4 (medium, high) Nuclear reaction vertex is placed at the step endpoint — secondaries systematically displaced downstream, possibly into the next zone

The interaction is sampled over the whole step (`ion_step_nuclear`,
`osh_transport_ion_step.c:976-1028`: survival probability over `ds_gcm2`),
but the vertex is effectively the step **endpoint**:

- `ion_step_commit()` runs first, moving the particle to the bent-path exit
  `q` and — for boundary-limited steps — **nudging it past the boundary**
  (`osh_transport_ion_step.c:1122-1149`).
- Only then are secondaries injected at `pool->x[slot]`
  (`:351-353` neutrons, `:373-375` ions, `:467-469` fragments), and
  sub-threshold recoils point-deposited there (`:458`).

So all nuclear secondaries/recoil deposits are displaced downstream by up to
one full step (mean ≈ step/2 ≈ mm scale at DELTAE=0.03 in water at
200 MeV), and for boundary-limited steps they are born marginally inside the
**next** zone/material — a reaction in water can deposit its recoil in the
adjacent bone voxel's zone. The primary's ionization `de` over the full step
is also scored even when ABSORB fired mid-step (`st.de` untouched,
`osh_transport_ion_step.c:1088-1092`), i.e. the primary ionizes through the
point where it ceased to exist.

- Impact: second-order for broad-beam dose, but a *systematic* downstream
  shift of the secondary-production field and of (n, fragment) source terms;
  affects the secondary-dose validation done in #212/#222 at the mm level and
  any thin-slab / interface geometry.
- Fix direction: sample the interaction depth from the survival law within
  the step (`s = -λ·ln(1-u·(1-e^{-ds/λ}))`), split the step there: score
  ionization for `[0,s]`, place the vertex and secondaries at `s`, and kill or
  redirect the primary from `s`.

### P-5 (low, high) Straggling clamp to [cutoff, e0] biases the sampled loss distribution for thin steps

`osh_transport_ion_step.c:931-935` clamps the post-straggling exit energy to
`[cutoff, e0]`. The upper clamp (loss ≥ 0) one-sidedly truncates the Gaussian:
for a step with mean loss μ and width σ the mean bias is
`σ·[φ(a) − a·(1−Φ(a))]` with `a = μ/σ`. For physics-limited steps
(`a ≈ 17` at DELTAE=0.03) this is nil, but **boundary-limited thin steps**
have `a ∝ √ds → 0`, where the bias approaches `σ·φ(0) ≈ 0.4σ` per step —
extra energy loss deposited in every thin zone crossing. In many-slab
geometries (ripple filters, stacked detectors) this accumulates as a
systematic overestimate of energy loss. The Landau tail clamp at `cutoff`
acts in the opposite direction for rare events.

- Fix direction: for boundary-limited steps either suppress straggling below
  a σ/μ threshold (classic condensed-history practice), or resample/fold the
  truncated tail so the mean is preserved.

### P-6 (low, medium) RNG child streams are seeded with the parent's *future* draws

`osh_rng_split()` (`osh_rng.c:90-120`) derives a child's (seed, stream) from
the ordinal-k window `[2k, 2k+1]` of a *copy* of the parent's stream — but the
parent keeps living: the very next physics draw the parent makes **is** the
u64 that became child-0's seed. This is structural seed reuse between a
secondary's stream initialisation and its parent's subsequent physics
sequence. PCG32's initialisation permutation makes an observable correlation
unlikely, but it is cheap to make impossible: derive child seeds by hashing
`(parent_state, ordinal)` with splitmix64 instead of consuming the parent's
own output window. (Everything else about the split design — drop-proof,
order-independent, non-advancing — is sound and well documented.)

### P-7 (low, high) Verified-OK notes for this subsystem

- Step-length logic: `preclip = min(boundary, CSDA, θ-limit)` with
  `e1_target` never below cutoff; the CSDA inversion via
  `energy_from_residual_range` is self-consistent (range-difference method,
  no midpoint-SP bias). Defensive `residual_range ≤ 0` kills are unreachable
  because `preclip_step_len ≤ ds_csda` by construction.
- Random-hinge MCS (Fippel/Soukup) correctly samples the hinge uniformly,
  applies full-step θ₀ at the hinge, re-queries the boundary in the scattered
  direction, and documents the O(θ₀²·L) chord error
  (`osh_transport_ion_step.c:722-847`).
- Highland θ₀ with macroscopic-path log correction (Gottschalk-style additive
  variance) is correctly implemented, including the substep/path-scale split
  (`osh_physics_scat_highland.c:8-44`); ξ, E_max, κ, λ̄ match the standard
  Vavilov formulary (`osh_physics_strag.c`).
- The Vavilov/Landau inverse-CDF evaluators are domain-clamped on κ, β², u
  with documented behaviour at the edges and a NaN guard on κ ≤ 0
  (`osh_physics_strag_vavilov.c:60-84`); `logf` precision use is justified in
  a comment and is sound.
- `osh_rng_double` yields 53-bit uniforms in [0,1); it **can return exactly
  0.0** — safe at current transport call sites (hinge, inverse-CDF clamps),
  but any future `-log(u)` free-path sampler must use `1-u` or reject 0
  (checked in neutron transport in a later section).
- Per-history seeding (`seed = f(rndseed, rndoffset + global index)`) is
  cursor-free at the fill site (`osh_transport_ion.c:192-197`), matching the
  #165 design; profile counters are per-worker (`wctx->profile`).
- Wavefront loop structure, pool compaction, progress reporting, and the
  clean-stop drain semantics are correct as documented; the step budget is
  global (`nstat × 1e6`) and a budget trip fails the whole run with
  OSH_ESTATE (`osh_transport_ion.c:233-251`) — acceptable, though a
  kill-the-particle-and-warn policy would be more forgiving late in a long
  run.

### P-8 (info) Float range tables: precision is adequate but close to the edge for DEMIN-limited steps

`range_csda` is stored as float (`osh_material_runtime.h:42`); transport
differences two lookups (`osh_transport_ion_step.c:678-683`). At R ≈ 30 g/cm²
(230 MeV proton in water) a float ULP is ≈ 2·10⁻⁶ g/cm²; the smallest
physics-limited steps (`demin`, default per-nucleon threshold × A) are
~10⁻⁴ g/cm², so the difference carries ~2 good decimal digits of headroom.
Fine today; a future finer DEMIN, or 1 GeV/u ions (R ~ 100 g/cm², ULP
~8·10⁻⁶) erode it. Worth a comment in the header at minimum.

---

*(Sections below are appended as the audit proceeds.)*
