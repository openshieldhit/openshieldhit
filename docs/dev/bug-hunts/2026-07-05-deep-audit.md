# Bug-hunting audit — 2026-07-05

!!! info "About this report"
    This is a periodic **bug hunt**: a deep, point-in-time audit that looks
    beyond what CI and normal code review catch. See
    [Bug hunt reports](index.md) for the methodology, the severity/confidence
    rating scheme shared by all reports, and the full list of reports.
    Findings below are AI-assisted and code-reading/experiment-based, but
    **not yet independently confirmed** by a maintainer — treat them as
    triage input, not settled fact.

| | |
|---|---|
| **Commit audited** | [`e3c619a`][e3c619a] on `main` (2026-07-05) |
| **Source** | [PR #248][#248] |
| **Scope** | Test infra · transport hot path & EM physics · scoring · nuclear/neutron transport · RNG & parallel readiness · GEMCA geometry · parsers/IO · architecture |
| **Follow-up** | 11 suggested issues (5 filed — [#255][#255], fixed by [#257][#257]; [#267][#267], fixed by [`a0f60a1`][a0f60a1]; [#279][#279], fixed by [#285][#285]; [#280][#280], fixed by [#287][#287]; [#299][#299], fixed by [#301]) — see [§8](#sec-8) |


Deep audit of `main` (e3c619a) prompted by PR [#239][#239] (sanitizer wiring) and issues
[#235][#235]/[#236][#236]/[#237][#237], with special attention to physics correctness, numerical
stability, parallelization readiness ([#161][#161], [#165][#165]–[#169][#169]), stddev/variance plans
([#169][#169]/[#195][#195]), and architecture ([TODO.md][todo-md]). Findings are numbered per subsystem
and tagged with **severity** (critical / high / medium / low / info) and **confidence**
(high = verified by running code or by construction; medium = strong code
reading; low = suspicion worth checking). Line numbers refer to `main` at
e3c619a unless a branch is named.

Method: manual code audit of the hot path outward, plus hands-on experiments:
debug + ASan/UBSan instrumented builds, full `ctest` runs, targeted standalone
harnesses. Every finding cites `file:line` evidence.

Status: **complete** (7 sections + executive summary; one commit per audit step).

---

## 0. Baseline, test infrastructure, and working-tree observations {: #sec-0 }

### T-1 (high, high) Integration-case registration turns any stray directory into a phantom failing test {: #t-1 }

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

### T-2 (high, high) The entire physics-validation tier (reference::idd_*) never runs — no reference curves are committed {: #t-2 }

All 13 `reference::idd_*` tests report **Skipped**: the runner
(`tests/reference/run_idd.cmake:41`) skips when `${CASE_DIR}/reference/` has
no data, and *no* `idd_*/reference/` directory exists in the tree (checked all
of `tests/reference/idd_*/`). So the IDD comparisons against SH12A/FLUKA that
[TODO.md][todo-md] cites as validation evidence ("distal 80–20% falloff 0.4522 vs
0.4533 cm", STRAGG 2 section) are not locked in by CI at all — a physics
regression in stopping power, straggling, or scattering would merge silently.

- Impact: no automated guard on physics output; the strongest claims in
  [TODO.md][todo-md] are unenforced.
- Fix direction: commit the reference curves (they are small 1-D IDD tables)
  for at least 70/150/200 MeV water, and add a CI job that runs
  `ctest -L reference` (possibly nightly if runtime is a concern; cases have
  TIMEOUT 1800 already).

### T-3 (medium, high) Tests write outputs into the current working directory (source-root pollution) {: #t-3 }

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
  the repo; leftovers then interact badly with [T-1](#t-1).
- Fix direction: root all test scratch under `OSH_TEST_TMPDIR` supplied by
  CMake (`CMAKE_CURRENT_BINARY_DIR`), and always clean up in a teardown path
  (or register with `FIXTURES_CLEANUP`).

### T-4 (info, high) PR #239 validated hands-on {: #t-4 }

Reproduced the PR's claims independently on `main` + manual
`-fsanitize=address,undefined -fno-sanitize-recover=all` build (GCC, Debug):
the **only** failures in the whole suite are exactly the two the PR fixes —
`unit::test_osh_scoring_parse_geometry` (stack-use-after-scope) and
`unit::test_osh_scoring_step` (stack OOB read via `nenergy = 0` table), both
aborting under ASan. Everything else is green with `detect_leaks=0`.
Conclusion: PR [#239][#239]'s test fixes and CI job are correct and worth merging
as-is; the sanitizer job would have caught both bugs at introduction time.

### T-5 (low, high) Working-tree / repo hygiene observations {: #t-5 }

- A file literally named `out_%%d.dat` sits at the repo root (written
  2026-07-03 by a run-control-branch binary): a `Filename` containing a
  printf-style `%d` pattern is used verbatim by the savers
  (`src/scoring/save/osh_scoring_save_ascii.c:70` opens `out->filename`
  as-is). If numbered periodic-dump filenames are ever intended ([#193][#193]), there
  is no templating support; if not intended, nothing warns the user that `%d`
  is meaningless. Low priority, but worth a decision.
- `.venv/**` (a full Python virtualenv, thousands of files) is currently
  *staged* in the local index on `main`. Not a repo bug, but one `git commit
  -a` away from becoming one; `.gitignore` should cover `.venv/`.

---

## 1. Transport hot path & electromagnetic physics {: #sec-1 }

Scope: `src/transport/osh_transport_ion_step.c`, `osh_transport_ion.c`,
`src/physics/atomic/*` (Bethe, Highland, Bohr/Vavilov/Landau straggling),
`src/material/runtime/*` (table build + lookup), `src/random/osh_rng.c`
(the samplers the hot path consumes).

### P-1 (high, high) Residual kinetic energy is silently deleted when a particle is killed at the energy cutoff {: #p-1 }

!!! success "Resolved"
    Filed as [#279][#279] and fixed by [#285][#285].  An ion or neutron killed
    at the transport energy cutoff (default or a raised `TCUT0`/`NEUTRLCUT`)
    now point-deposits its residual kinetic energy at the kill site instead of
    discarding it, gated on positive-density material and attributed to the
    particle's own generation, with an energy-conservation regression test
    (`tests/unit/test_osh_cutoff_deposit.c`) for the ion side and a companion
    test for the neutron side.  The related unknown-species kill noted below
    is also fixed: it now increments a dedicated `n_unknown_dropped` counter
    and WARNs once per run instead of failing silently.  The `residual_range
    <= 0` kills and the EMIN/effective-cutoff fragment-injection mismatch
    (both noted below) are untouched by #285 and remain open.

Every ion ends its life at `ion_step_setup()`:

- `osh_transport_ion_step.c:539-543` — `if (ctx->e0 <= ctx->cutoff) { pool->e[slot] = 0.0; ctx->done = 1; return; }` — no scoring call of any kind.
- The exit-energy clamps at `osh_transport_ion_step.c:897-898` and `:934-935` park a stopping particle at exactly `cutoff`, so this kill fires for essentially **every** track that ranges out (rather than escaping).
- `cutoff = max(OSH_MATERIAL_RUNTIME_EMIN, params->tcut, material_rt->emin) × A` (`cutoff_total_energy`, `osh_transport_ion_step.c:1156-1177`; `EMIN = 0.025` MeV/u, `osh_material_runtime.h:88`).

The wavefront loop (`osh_transport_ion.c:294-295`) merely compacts dead
entries; nothing deposits the residual kinetic energy. Contrast with the care
taken for **sub-threshold nuclear fragments**, which do get a point deposit
(`ion_point_deposit`, `osh_transport_ion_step.c:456-458`, issue [#179][#179]) — the
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

### P-2 (high, high) Non-default isotopes are transported with the representative isotope's range table and rest mass — deuteron range wrong by ≈2× {: #p-2 }

!!! success "Resolved"
    Filed as [#267][#267] and fixed by [`a0f60a1`][a0f60a1].  Transport still
    shares stopping-power columns by `Z`, but CSDA range is now scaled by
    `A_actual/A_table`, isotope-specific kinematics use `part->mass`, and a
    regression unit covers non-default H/He isotopes in a non-unit-density
    material.

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

### P-3 (medium, high) Bohr Gaussian straggling omits the relativistic variance factor — σ ~11% low at 200 MeV; discontinuity at the STRAGG 2 κ=10 dispatch {: #p-3 }

!!! success "Resolved"
    Filed as [#280][#280] and fixed by [#287][#287].  `osh_physics_strag_sigma()`
    now takes `beta2` and applies the `γ²(1−β²/2)` relativistic correction to
    the Bohr variance, and both call sites in `osh_transport_ion_step.c`
    (STRAGG 1 and the STRAGG 2 Gaussian branch) pass `beta2` through, so the
    κ=10 dispatch no longer has a variance discontinuity.

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
  [TODO.md][todo-md] didn't see it — but the branch is also taken for thick steps at
  high energy, e.g. coarse DELTAE in thick slabs.)
- Fix: multiply the variance by `γ²(1−β²/2)`; β², γ are already computed at
  the STRAGG 2 call site and cheap to add at the STRAGG 1 site. The comment
  in `osh_physics_strag_gauss.c:17-19` ("kept as a literal so the Gaussian
  result is unchanged") suggests the omission is legacy-compat, but nothing
  documents it as a deliberate physics approximation.

### P-4 (medium, high) Nuclear reaction vertex is placed at the step endpoint — secondaries systematically displaced downstream, possibly into the next zone {: #p-4 }

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
  affects the secondary-dose validation done in [#212][#212]/[#222][#222] at the mm level and
  any thin-slab / interface geometry.
- Fix direction: sample the interaction depth from the survival law within
  the step (`s = -λ·ln(1-u·(1-e^{-ds/λ}))`), split the step there: score
  ionization for `[0,s]`, place the vertex and secondaries at `s`, and kill or
  redirect the primary from `s`.

### P-5 (low, high) Straggling clamp to [cutoff, e0] biases the sampled loss distribution for thin steps {: #p-5 }

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

### P-6 (low, medium) RNG child streams are seeded with the parent's *future* draws {: #p-6 }

!!! success "Resolved"
    Filed as [#299][#299] and fixed by [#301].  `osh_rng_split()`
    now derives the child's `(seed, stream)` by hashing the parent's *internal
    state* words — which the engine's output function permutes before ever
    emitting — keyed by the ordinal, through the existing `rng_mix_stream`
    SplitMix64 machinery.  A child seed therefore can never coincide with a
    value the parent will itself draw next, closing the structural reuse.  The
    split stays non-advancing, drop-proof, and order-independent (and is now
    O(1) in the ordinal rather than O(ordinal)); regression tests
    `test_split_not_seeded_from_parent_output` and
    `test_split_xoshiro_properties` lock the non-reuse and per-engine
    behaviour.

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

### P-7 (low, high) Verified-OK notes for this subsystem {: #p-7 }

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
  [#165][#165] design; profile counters are per-worker (`wctx->profile`).
- Wavefront loop structure, pool compaction, progress reporting, and the
  clean-stop drain semantics are correct as documented; the step budget is
  global (`nstat × 1e6`) and a budget trip fails the whole run with
  OSH_ESTATE (`osh_transport_ion.c:233-251`) — acceptable, though a
  kill-the-particle-and-warn policy would be more forgiving late in a long
  run.

### P-8 (info) Float range tables: precision is adequate but close to the edge for DEMIN-limited steps {: #p-8 }

`range_csda` is stored as float (`osh_material_runtime.h:42`); transport
differences two lookups (`osh_transport_ion_step.c:678-683`). At R ≈ 30 g/cm²
(230 MeV proton in water) a float ULP is ≈ 2·10⁻⁶ g/cm²; the smallest
physics-limited steps (`demin`, default per-nucleon threshold × A) are
~10⁻⁴ g/cm², so the difference carries ~2 good decimal digits of headroom.
Fine today; a future finer DEMIN, or 1 GeV/u ions (R ~ 100 g/cm², ULP
~8·10⁻⁶) erode it. Worth a comment in the header at minimum.

---

## 2. Scoring: estimators, binning, normalization, variance/parallel readiness {: #sec-2 }

Scope: `src/scoring/runtime/` (step scorer, accumulator, postprocess, shadow,
snapshot), `src/scoring/save/` (ASCII/BDO), plus the deposit call sites in
transport.

### S-1 (high, high) Particle weight `st->wt` is ignored by every deposit site {: #s-1 }

The weight is carefully plumbed end-to-end — `pool->wt` at fill, copied to
secondaries (`osh_transport_ion_step.c:358,380,474`), into the neutron pool
and back (`osh_transport_neutron.c:103,146`), and into every scoring step
(`step_from_pool`, `osh_transport_ion_step.c:1311`; point deposits
`osh_transport_ion_step.c:237`) — and then **no scorer reads it**:
`grep -n 'st->wt' src/scoring/runtime/*.c` matches nothing. `score_group_energy`
deposits `st->de·frac` (`osh_scoring_step.c:706-711`), dose deposits
`path_len·vol_inv·dose_scale` (`:891-896`), fluence `path_len·vol_inv`
(`:778-782`), DLET numerator/denominator `let·w`, `w = de·path/score_len`
(`:986-990`) — all weight-free. Same in `osh_scoring_point.c`.

Today this is invisible: the beam always sets `pool->wt[slot] = 1.0` (SOBP
spot weights are absorbed into the CDF *sampling*, `osh_beam_runtime.c:231`),
and neutron transport has no implicit capture yet. But:

- **MCPL phase-space import is planned** ([TODO.md][todo-md] Beam, [#41][#41]) — MCPL records
  carry per-particle weights; the moment that lands, every scorer silently
  produces wrong results with no test to catch it.
- Any future variance-reduction (implicit capture in the neutron family is
  the standard technique for exactly this code's regime) hits the same wall.
- The stddev work ([#169][#169]) must also decide *now* whether batch `weight` means
  "history count" or "sum of statistical weights" — the current accumulator
  bookkeeping (`acc->weight`) assumes unit weights.
- Fix direction: multiply `st->wt` into the deposited value at each
  `score_group_*` (numerator *and* denominator for the two-pass LET/Qeff
  pages), add a unit test transporting a `wt = 2` primary and asserting
  doubled deposits; decide the weighted-variance contract in [#169][#169] at the same
  time.

### S-2 (medium, high) Nuclear-kill steps evaluate LET/Qeff/differential axes at half the true midpoint energy {: #s-2 }

For a primary destroyed by an inelastic nuclear event, `ion_step_commit()`
signals death via `st.q[3] = 0.0` (`osh_transport_ion_step.c:1088-1092`).
But every scorer that needs a midpoint energy computes it as
`0.5·(st->p[3] + st->q[3])` (`osh_scoring_step.c:441` `compute_step_let_medium`,
`:548` `compute_step_qeff` path, `:566` `diff_step_val`): for the absorption
step this is `e0/2` instead of `(e0 + exit)/2`. LET at half the energy is
substantially larger (for 100 MeV protons, S(50 MeV) ≈ 1.8× S(100 MeV)), so:

- DLET/TLET maps acquire a *systematic high bias* in exactly the regions
  where nuclear reactions matter, weighted by the (non-small) ionization
  `de` of those steps;
- differential (spectrum) scorers mis-bin those steps' deposits toward lower
  energy / higher LET.
- Elastic events have the milder variant (`q[3]` = post-elastic energy,
  `osh_transport_ion_step.c:1096`, instead of the pre-elastic CSDA exit).
- Fix direction: carry the ionization exit energy separately in `struct step`
  (e.g. keep `q[3] = exit_energy` and signal death via an explicit flag or
  via `pool->e` only), so scoring midpoints are physics-consistent.

### S-3 (medium, high) Dose uses the transport voxel's density at step start for all scoring-mesh crossings {: #s-3 }

`score_group_dose` computes one `base_scale = st->de/(score_len·st->rho)`
per step (`osh_scoring_step.c:838-841`), where `st->rho` is the density of the
transport zone/voxel where the step *began*. Every scoring-mesh crossing of
that step — which may span several CT voxels of very different density —
gets that single ρ. This is the known TODO item ("DOSE scoring with per-voxel
density weighting for CT geometries"), confirmed in code; the error is largest
exactly at lung/soft-tissue/bone interfaces where dose accuracy is most
scrutinised. (Also note the step is scored on the chord `p→q` while the
transport step may have been hinge-bent — documented O(θ₀²·L), fine.)

### S-4 (low, high) Differential-axis edge semantics: both edges exclusive, out-of-range silently dropped, ASCII spectra now divided by bin width {: #s-4 }

`diff_axis_bin` (`osh_scoring_step.c:407-421`) treats the axis as the *open*
interval `(lo, hi)`: a value exactly at `lo` (or `hi`) is discarded via the
out-of-range sentinel, and there are no under/overflow bins anywhere, so the
discarded tally is invisible. The original audit also found that differential
outputs were written as per-bin counts; issue #215 moved bin-width
normalisation into scoring postprocess, so ASCII and BDO2019 writers now receive
differential quantities.

### S-5 (info, high) Verified-OK: the accumulator/merge layer is sound and correctly shaped for #169 {: #s-5 }

- `osh_scoring_accumulator_merge` implements the pairwise Welford/West M2
  combine with the Schubert–Gertz cross-term, empty-batch identities, ordered
  correctly before the additive fold, with consistency guards
  (`osh_scoring_accumulator.c:100-210`). This is the *right* variance
  representation (batch means / M2), not the precision-hazardous naive Σx²
  — the [#169][#169] decision can simply ratify what is scaffolded.
- `data_var`/`data2_var` are pure scaffold today: nothing allocates them
  (only `alloc(data, data2)` exists), so merge/zero paths on them are no-ops.
  One inconsistency to fix when wiring them: `osh_scoring_accumulator_alloc`
  cannot create them, and `osh_scoring_accumulator_rescale`
  (`osh_scoring_accumulator.c:67-75`) rescales only `data` — currently dead
  code (zero call sites); if ever used on a two-pass (LET) page before the
  divide, it would corrupt the ratio. Delete it or make it scale all arrays
  consistently.
- Normalization chain verified: BDO stores raw sums + `OSHBDO_RT_NSTAT` (+
  the [#195][#195] completeness label) for post-hoc/merge use
  (`osh_scoring_save_bdo2019.c:9-23,119-131`); ASCII divides NORM/SUM pages by
  `completed_nstat` and correctly leaves AVER (DLET/TLET) undivided
  (`osh_scoring_save_ascii.c:329-341`); `completed_nstat` is the exact drained
  count, including clean-stop paths (`osh_simulation.c:645-650`,
  `osh_transport_ion.c:329-333`).
- The deposit seam `osh_score_deposit` (issue [#158][#158]) is in place at all sites
  (plain `+=`, `osh_scoring_accumulator.h:111-113`), so per-worker
  redirection ([#166][#166]) is mechanical.
- Crossing scratch is caller-owned with capacity checked before use
  (`osh_scoring_step.c:135-137,155-157` return OSH_ESTATE rather than
  overflow); capacity `Σ nbins` per geometry is a safe DDA bound.

---

## 3. Nuclear physics & neutron transport {: #sec-3 }

Scope: `src/physics/nuclear/` (handler, Tripathi, pp/pA elastic, abrasion,
Fermi break-up), `src/physics/neutron/` (JEFF Tier-1 tables, reaction
sampling), `src/transport/osh_transport_neutron.c`, `osh_neutron_pool.c`.

### N-1 (high, high) Neutron interactions deposit **zero** energy — and n–p elastic recoil energy vanishes entirely {: #n-1 }

`apply_event()` (`osh_transport_neutron.c:166-202`) discards
`ev->local_deposit_mev` on **every** channel (comments say "deposit when
scoring wired"), and drops all non-neutron secondaries. The n–H elastic case
is actively wrong rather than just unwired: `do_elastic()` builds the recoil
proton secondary and **zeroes the local deposit** on the assumption the
proton carries the energy (`osh_neutron_reaction.c:91-108`,
`ev->local_deposit_mev = 0.0; /* proton secondary carries the recoil */`) —
then `apply_event`'s ELASTIC branch ignores `ev->secondaries` entirely
(`osh_transport_neutron.c:171-178`). Net effect:

- The dominant neutron dose mechanism in tissue (n–p recoil) contributes
  **nothing** — the energy is neither deposited nor transported; energy
  conservation is broken at every n–H collision.
- Capture, (n,p)/(n,α), compound-residual channels likewise deposit nothing.
- Neutron transport today only attenuates and redirects neutrons (fluence
  scorers work); every energy/dose scorer sees zero from the neutron family.
- The neutron **cutoff kill** (`osh_transport_neutron.c:287-291`, default
  `ncut = 1 keV`, `osh_transport.h:48`) also discards the remaining energy —
  same class as [P-1](#p-1). **Resolved by [#285][#285]** (tracked under
  [P-1](#p-1)/[#279][#279]): the cutoff kill now point-deposits the residual
  in positive-density material. The interaction/reaction deposits below
  (capture, (n,p)/(n,α), compound-residual, n-p elastic recoil) are a
  separate mechanism and remain open.
- **[TODO.md][todo-md] misdescribes this**: the "What is implemented" list says "elastic
  (isotropic CM; n-p recoil proton returned)" and "(currently deposited
  locally)" for charged secondaries. Neither is true in code — only the
  *open-items* list ("Local neutron energy deposits scored") matches reality.
  Anyone reading [TODO.md][todo-md] (or the [#154][#154] merge summary) will assume neutron
  dose exists.
- Fix direction: point-deposit `local_deposit_mev` at the interaction site
  now (the [#179][#179] point-deposit path already exists and is exactly this
  mechanism); feed the recoil proton into the ion pool when ion-feedback
  lands. Until then, correct [TODO.md][todo-md].

### N-2 (high, high) Tier-1 cross sections are clamped at the 20 MeV table edge — fast neutrons interact with σ(20 MeV) {: #n-2 }

`interp_one()` clamps: `e ≥ egrid[last] → arr[last]`
(`osh_neutron_xsec.c:84-87`), and Tier-1 entries always win over Tier-2 for
tabulated nuclides (`osh_neutron_xsec_lookup`, `:264-276`). Abrasion neutrons
from a 130–230 MeV proton beam populate the spectrum up to near beam energy,
far above the 1 meV–20 MeV JEFF grid. For n–p, σ_tot(20 MeV) ≈ 0.48 b vs
σ_tot(100 MeV) ≈ 0.075 b — a ×6 overestimate; O-16 is similar (×~5). Fast
neutrons therefore interact much too often and are attenuated on far too
short a mean free path — their fluence map is systematically compressed
toward the production site. (Today the *dose* consequence is masked by [N-1](#n-1);
fixing [N-1](#n-1) without [N-2](#n-2) will bake this bias into neutron dose.)

- Fix direction: above the Tier-1 grid, dispatch to the Tier-2
  Tripathi σ_R + geometric σ_el path (already implemented for missing
  nuclides), or extend the condensed tables; at minimum warn once when a
  lookup exceeds the grid.

### N-3 (medium, high) Lazily-initialized `static` species cache — explicit single-thread assumption inside the neutron family {: #n-3 }

`osh_neutron_reaction.c:20-33`: `static struct particle s_proton, s_alpha;
static int s_species_ready;` filled by `ensure_species()` on first call,
commented "all transport is single-threaded". The WIP parallel driver
(branch `feat/parallel-threads`) drains per-worker neutron pools
**concurrently** — this becomes a data race (benign-looking but UB, and
TSan-fatal) the day `--threads` ships. Same pattern to audit elsewhere
(`osh_neutron_xsec` warn-once state is per-instance and fine; Molière tables
are built once on the setup thread before transport starts —
`osh_transport_ion.c:562-566` — which is the right pattern).

- Fix: initialise at handler compile time (setup thread), or make them
  `const` compile-time tables.

### N-4 (low, medium) Number densities use the integer mass number, not atomic mass {: #n-4 }

`neutron_sigma_tot_cm()`: `nd_i = mass_fraction·ρ·N_A / A_int`
(`osh_transport_neutron.c:58`); the nuclear handler's rate loops use
`ai = elems[i].a` or the `z*2` fallback (`osh_nuclear_handler.c:305,344`).
Using A instead of the isotope atomic mass in g/mol overestimates hydrogen
number density by ~0.8% (1 vs 1.008) and <0.1% for heavier isotopes —
systematic, though small against the model's other approximations. The
`a=0 → 2Z` fallback deserves a diagnostic if it can ever fire in practice.

### N-5 (info) Physics-approximation register (documented, not bugs — should stay visible) {: #n-5 }

- Elastic scattering is isotropic in CM (P0) at all energies
  (`osh_neutron_reaction.c:60-62`); real n–A angular distributions are
  strongly forward-peaked above ~1 MeV, so energy transfer per collision (and
  hence moderation rate) is overestimated at high energy.
- Non-relativistic elastic kinematics (fine below 20 MeV; not above — ties
  into [N-2](#n-2)).
- No thermal treatment below 1 eV ([#178][#178], known).

### N-6 (info, high) Verified-OK notes {: #n-6 }

- Free-path sampling guards `u = 0` explicitly with the HUGE_VAL fall-through
  (`osh_transport_neutron.c:338-345`) — correct and well-commented.
- The Kopylov N-body phase-space generator mirrors G4FermiPhaseSpaceDecay:
  momentum conservation by construction (back-to-back subsystem CM + boosts),
  fragment 0 closed on-shell (`osh_nuclear_fermi_breakup.c:424-508`);
  `beta_kopylov` rejection envelope is correct.
- pp-elastic channel bookkeeping is consistent: kinetic energy split
  `e1 + e2` with equal-mass lab kinematics, Møller-convention swap keeps the
  primary on the forward proton (`osh_nuclear_handler.c:370-401`).
- Channel competition (`rate_inel + rate_pp + rate_pa`), the struck-element
  selection proportional to per-element hazard, and the RNG-stream-stability
  comment (channel draw only when elastic competes) are all sound.
- Neutron pool: `n_created`/`n_dropped` accounting, per-slot RNG streams via
  `osh_rng_split`, `prim_idx`/`gen` propagation (good for future per-primary
  statistics), wavefront batching bounded by shared scratch capacity.

---

## 4. RNG & parallelization readiness (#161, #165–#169) {: #sec-4 }

Scope: `src/random/`, worker context and history loop (partly covered in §1),
whole-tree shared-mutable-state sweep, and the **WIP pthreads driver on the
local branch `feat/parallel-threads`** (`osh_transport_parallel.c` — reviewed
because it is the concrete next step of [#161][#161]).

### R-1 (info, high) RNG core verified against references {: #r-1 }

- PCG32 matches O'Neill's reference exactly: state advance
  `state·6364136223846793005 + (inc|1)`, XSH-RR output, and the canonical
  `srandom` init sequence (`osh_rng_pcg32.c:15-47`).
- Per-history stream derivation `rng_mix_stream(seed, hist_index, purpose)`
  is a proper hash: golden-ratio Weyl spacing on the index axis, distinct odd
  constant on the purpose axis, Stafford Mix13 finaliser, all with excellent
  in-source rationale (`osh_rng.c:25-88`) and an empirical disjointness test.
  This is genuinely parallel-ready seeding ([#165][#165]'s contract holds).
- `osh_rng_double` = 53-bit in [0,1); can return exactly 0 — the one `-log(u)`
  consumer guards it explicitly (`osh_transport_neutron.c:338-345`).
- Marsaglia-polar gaussian with per-stream spare cache — safe because each
  pool slot owns its stream; no cross-thread hazard as long as streams stay
  slot-owned.
- Known-and-accepted: PCG32 "streams" (different `inc`) are not provably
  independent sequences; with hashed 64-bit stream ids and per-history
  reseeding this is a non-issue in practice, but worth one line in the RNG
  doc.

### R-2 (low, medium) `osh_rng_split` seeds children from the parent's *future* output window {: #r-2 }

!!! success "Resolved"
    Same fix as [P-6](#p-6) — filed as [#299][#299], fixed by [#301].

Detailed as [P-6](#p-6) in §1: child ordinal k consumes the parent-copy's draws
[2k, 2k+1] as (seed, stream) — the same u64s the parent itself will consume
for its next physics decisions (`osh_rng.c:90-120`). Structural seed reuse;
practically masked by PCG's init permutation; cheap to eliminate by hashing
`(parent_state, ordinal)` with the existing `rng_mix_stream` machinery. The
drop-proof/order-independent design itself is excellent.

### R-3 (medium, high) Shared-mutable-state inventory — what actually blocks `--threads N` {: #r-3 }

Complete sweep of file-scope mutable statics in `src/` (excluding apps
parse-dispatch tables, which are logically-const and should just gain
`const`):

| State | Where | Hazard under threads | Needed change |
|---|---|---|---|
| `s_proton`, `s_alpha`, `s_species_ready` | `osh_neutron_reaction.c:23-33` | lazy init inside neutron drain — data race (see [N-3](#n-3)) | init at handler compile (setup thread) |
| `s_n_warned` | `osh_nuclear_compound.c:17` | warn-once counter incremented on hot path — benign race, still UB | per-handler counter or C11 atomic |
| `osh_moliere_inited` + tables | `osh_physics_scat_moliere.c:38` | initialised on the setup thread via `validate_transport_modes` (`osh_transport_ion.c:562-566`) — **safe in the serial driver**, but the WIP parallel driver has each worker call `osh_transport_ion_run_range` → `validate_transport_modes` → `osh_physics_moliere_init` concurrently → first-call race | init once in the parallel driver before spawning workers |
| `g_stop`, `g_dump` | `osh_apps/osh_signals.c:77-85` | `volatile sig_atomic_t` — correct pattern | none |

Beyond statics, the known open items hold as documented: master accumulator
deposits ([#166][#166] — the seam exists), and everything else on the hot path
(pools, scratch, profile, RNG) is already per-worker-ownable. The codebase is
closer to thread-ready than [#161][#161]'s "honest gap" list implies — the gaps are
the table above plus [#166][#166].

### R-4 (high, high) WIP parallel driver (`feat/parallel-threads`) — concrete bugs to fix before it merges {: #r-4 }

`src/transport/osh_transport_parallel.c` at branch head 83b3575:

1. **`pthread_join` on a never-created thread.** If `pthread_create` fails
   the code sets `w->rc = OSH_ESTATE` (`:579-583`) but the join loop
   unconditionally joins `w->thread` for every non-empty slice (`:586-592`)
   — joining an uninitialized `pthread_t` is UB. Track a `spawned` flag.
2. **Error-path leaks in worker setup.** The neutron-pool failure path
   (`:519-531`) and fragment-pool failure path (`:533-546`) free only worker
   *i*'s resources and `free(workers)`, leaking workers `0..i-1`'s pools,
   accumulator sets and scratch (the earlier failure paths do clean up
   predecessors — the last two branches were evidently written later).
3. **Per-worker neutron pools are each sized to the full run.**
   `neutron_capacity = nstat` (`:443-444`) and every worker gets one
   (`:520-521`): memory scales O(nthreads·nstat) — for the [#161][#161] benchmark's
   200-core scenario this is a blow-up. Size by the worker's range
   (`hist_hi - hist_lo`), matching the serial pool-accumulation argument.
4. **Single ion-pass → single neutron-drain per worker.** The serial path
   runs the *family scheduler* (`run_families_over_range`,
   `osh_transport.c`) with has-work discipline; the worker function instead
   hard-codes one ion range + one neutron drain (`:352-381`). Today that is
   equivalent (neutron→ion feedback does not exist — [N-1](#n-1)), but the moment
   ion-feedback lands, the parallel path silently loses the fed-back
   charged secondaries. Workers should call the same family scheduler over
   their range, not a hand-rolled subset of it.
5. Cosmetics with teeth: `parallel_alloc_worker_fragment_pool` "clones"
   nothing (fine only because `struct osh_fragment_pool` is counters-only —
   worth a comment); `worst_rc` keeps the *last* failure, not the worst;
   the smoke tests write outputs to the CWD and don't clean up (see [T-3](#t-3)).

Also verified on the branch: the per-worker accumulator alloc mirrors
`master_acc` shape including `data2` presence (`:92-126`), the merge loop
uses `osh_scoring_accumulator_merge` per page (`:608-613`), and worker slices
`[nstat·i/n, nstat·(i+1)/n)` are exact and disjoint — the core of the design
is right.

### R-5 (info, high) Reproducibility (#168) status {: #r-5 }

- Per-history streams are pure functions of the global index — verified at
  the fill site (`osh_transport_ion.c:192-197`) and by the seeding design
  ([R-1](#r-1)); range splitting/replay reproduces canonical streams.
- Bitwise reproducibility across different `nthreads` is impossible by
  design (FP summation order in merge); for *fixed* worker count the branch
  driver's ordered merge (worker 0..n-1, `:600-620`) is deterministic. The
  [#168][#168] CI gate should therefore pin `nthreads` and assert byte-identical
  output, plus a statistical-tolerance check across thread counts — exactly
  what the branch's `test_serial_matches_parallel` (rel 1e-9) already
  sketches.
- One genuine nondeterminism source to keep out of the physics path:
  time-based decisions. The run-control stop/dump checks read wall time but
  only at wavefront-pass boundaries and only affect *how many* primaries run
  (counted exactly), never a history's own randoms — sound.

---

## 5. Geometry (GEMCA) & shared numerics {: #sec-5 }

Scope: `src/gemca/` (quadric distances, dispatch), boundary nudge policy,
`src/common/raytrace/` (scoring DDA), material runtime lookup edge cases.

### G-1 (critical, high) Missing factor 2 in the linear coefficient — ellipsoid, elliptic-cylinder and cone surface distances are mathematically wrong {: #g-1 }

**Filed:** [#255][#255] (covers this finding and [G-2](#g-2)) · **Fixed by:** [#257][#257].

`_quadratic_solver(a, b, c)` solves `a·t² + b·t + c = 0` with
`d = b² − 4ac`, roots `(−b ± √d)/(2a)` (`osh_gemca2_dist.c:482-510`), i.e. it
expects the **full** linear coefficient. The call sites disagree:

| surface | b passed | correct? |
|---|---|---|
| `_dist_sphere` (`:353-363`) | `2·(cp·p)` | ✔ |
| `_dist_cyl` CYLZ (`:380-390`) | `2·(cpx·px + cpy·py)` | ✔ |
| `_dist_elipcyl` ELLZ (`:406-417`) | `cpx·px/ra² + cpy·py/rb²` | ✗ missing ×2 |
| `_dist_cone` CONE (`:433-446`) | `cpx·px + cpy·py − t·cpz` | ✗ missing ×2 |
| `_dist_ellipsoid` ELLIPSOID (`:459-470`) | `Σ cpᵢ·pᵢ/rᵢ²` | ✗ missing ×2 |

Numerical proof by internal inconsistency: a unit circular cylinder queried
from `p = (0.5, 0, 0)` along `+x` returns **0.5** via `_dist_cyl`
(b = 1: d = 4, roots −1.5, **0.5**) but **0.6514** via `_dist_elipcyl` with
`ra² = rb² = 1` (b = 0.5: d = 3.25, roots −1.151, **0.651**) — the same
surface, 30% apart. Every ray–ELL/ELLZ/CONE intersection with `b ≠ 0` is
wrong by an error that varies along the ray (only center-through rays are
correct).

- Impact: any geometry using `ELL`, elliptic cylinders, or truncated cones
  (`_dist_surface` dispatch, `osh_gemca2_dist.c:158-172`) transports
  particles across the wrong boundary positions — dose systematically
  misplaced, zones effectively deformed. **No test or example in the tree
  references these bodies** (`grep -rl 'ELL|TRC|REC' tests/ examples/` is
  empty), which is why the suite is green.
- Additional suspicion on `_dist_cone` (`:443-444`): the quadratic's `a`
  has `+cpz²/rb²` — a cone's z-term should enter with a **negative** sign
  (`x² + y² − k²(z−z₀)² = 0`); the parametrization via `t = (pz − ra²)/rb²`
  is opaque enough that this needs a dedicated unit test rather than reading.
- Fix direction: pass `2·(…)` at the three sites (or better: make the solver
  take the half-coefficient `h` and use the numerically stable
  `q = −(h + sign(h)·√(h²−ac)); t₁ = q/a; t₂ = c/q` form — see [G-2](#g-2)); add
  per-surface distance unit tests (inside/outside/tangent/grazing) for
  every body type, which would have caught this immediately.

### G-2 (medium, high) Cancellation-prone quadratic root formula near surfaces {: #g-2 }

**Filed:** [#255][#255] (tracked together with [G-1](#g-1)) · **Fixed by:** [#257][#257].

Even the correct call sites use `(−b + √d)/(2a)` (`osh_gemca2_dist.c:507-509`).
When the particle sits just off a surface (`c ≈ 0` — precisely the
post-nudge state every boundary crossing produces), the small root suffers
catastrophic cancellation; a garbage small-positive root re-detects the
surface at distance ~0 and costs an extra nudge/wavefront round (bounded by
the step budget, but wasteful and noisy near grazing incidence). The stable
citardauq/sign-split form costs nothing and removes the failure mode.

### G-3 (info, high) Boundary-nudge policy: absolute ε = 1e-8 cm is sound for the stated domain {: #g-3 }

`OSH_TRANSPORT_BOUNDARY_EPS = 1.0e-8` cm applied along the direction
(`osh_transport_boundary.h:13-27`). At clinic-scale coordinates (|x| ≤ 10³ cm)
double roundoff in the distance solvers is ≤ ~1e-11 cm, so the nudge
dominates it — no stuck-particle hazard from magnitude alone. Two residual
notes: (a) grazing rays get a *normal* displacement of ε·cosθ → repeated
re-nudges are possible but each is caught by the `boundary_ds ≤ EPS` guard
(`osh_transport_ion_step.c:571-578`) and bounded by the step budget; (b) the
ε path segment per crossing is untracked (≤1e-8 cm per boundary — utterly
negligible, worth a one-line comment at most).

### G-4 (info, high) Verified-OK: scoring DDA and index plumbing {: #g-4 }

- The Jacobs alpha-parametric traversal (`osh_raytrace_jacobs_msh.c`) is
  correct and robust: slab clip, ε-inside entry-voxel probe with index clamp,
  tie-tolerant multi-axis advance at corners (`same_crossing`, 64·DBL_EPSILON
  relative), segment count bounded by `Σnᵢ` which matches the caller's
  scratch-capacity check — no negative-index truncation bug (entry index is
  clamped after the cast).
- Zone sentinel constants are consistent (`OSH_ZONE_INDEX_INVALID` ==
  `OSH_GEMCA_ZONE_INDEX_INVALID` == `(size_t)-1`).
- CT density lookup clamps HU to [−1000, 1600]
  (`osh_material_runtime.h:195-204`): correct for the LUT size, but **metal
  implants (HU 2000–3000+) silently become HU 1600** (~dense bone). Worth a
  once-per-load warning when the CT contains HU > 1600 (ties into the
  Schneider LUT domain; clinically relevant for implant patients).
- `src/gemca/osh_gemca2_dist.c` uses `//` comments and `// TODO vectorize
  me` markers in production code — violates DEVELOPER.md §4.1 (style; noted
  because CLAUDE.md calls it a hard rule and clang-tidy evidently does not
  catch it).

---

## 6. Parsers & I/O (beyond the known #235/#236/#237) {: #sec-6 }

Scope: `src/dicom/`, `src/apps/osh/*_parse.c`, `src/common/osh_readline.c`,
save writers. Overall verdict: this layer is in noticeably good shape — the
findings below are minor relative to §1–§5.

### IO-1 (low, high) DICOM walker length guard can wrap on 32-bit `size_t` {: #io-1 }

`osh_dicom_walk()` guards each tag with `if (pos + length > size)`
(`osh_dicom_parse.c:124-131`). `length` is attacker-controlled `uint32_t`; on
an LP32/ILP32 target (`size_t` 32-bit) `pos + length` can wrap and bypass the
guard → OOB read of up to 4 GB offsets in the tag callback. 64-bit builds
are safe. Fix: `if (length > size - pos)` (pos ≤ size holds by loop
invariant).

### IO-2 (low, high) CT slice with rejected Pixel Data proceeds with `pixels == NULL` {: #io-2 }

In `_tag_cb`, if `_slice_pixel_count()` fails or the declared length is short
(`osh_dicom_ct.c:99-107`), the slice is silently kept with `pixels = NULL`;
the volume assembler then dereferences `slices[i].pixels`. A malformed CT
(pixel-data length lied about) is a NULL-deref crash rather than a
diagnostic. The [#235][#235] fix (uniform-dims validation) should also reject
pixel-less slices — same validation loop, one more condition.

### IO-3 (low, medium) RTDOSE `NumberOfFrames` parsed with unchecked `atoi` {: #io-3 }

`rd->n_frames = atoi(tmp)` (`osh_dicom_rtdose.c:89`): garbage/negative values
flow onward (the frame-offset capacity logic happens to clamp ≥ 1, and a
mismatch only warns). Validate `n_frames ≥ 1` and consistency with
rows/cols/frame-offset count before assembling the dose grid.

### IO-4 (low, medium) `osh_readline` silently splits over-long lines {: #io-4 }

`fgets(buff, OSH_MAX_LINE_LENGTH, …)` (`osh_readline.c:47`): an input line
longer than the buffer is processed as two lines with no warning — the
eventual parse error (if any) points at a confusing place. Detect a missing
`\n` and emit a truncation diagnostic with the line number.

### IO-5 (info, high) Verified-OK {: #io-5 }

- All `sscanf` uses in the apps parse layer are width-bounded with the
  `%c`-sentinel trailing-garbage idiom (`osh_material_parse.c:259,493,674,…`)
  — consistent and safe.
- The DICOM walker bounds-checks every tag header and value, handles long/short
  VR forms, and stops cleanly on undefined-length sequences.
- OOM policy is uniform (`osh_abort_oomf`), so no unchecked-malloc NULL
  derefs in the parse layer.
- One environmental note: numeric parsing (`strtod`/`sscanf %lf`) assumes the
  C locale. The CLI never calls `setlocale`, so it is fine standalone — but if
  `libopenshieldhit` is embedded in a host that sets `LC_NUMERIC` (Qt/GUI
  apps do), every `geo.dat`/DICOM DS parse breaks. Worth pinning with
  `strtod_l`/`_configthreadlocale` or documenting the embedding requirement.

---

## 7. Architecture, layering & plans (TODO.md, #161/#169/#195) {: #sec-7 }

### A-1 (medium, high) Layering rule violations: two dead `transport/` includes in lower layers {: #a-1 }

The declared rule (module READMEs, DEVELOPER.md §9; "runtime/ must not depend
on transport/") is violated twice — both **dead** includes (no transport
symbol is used in either file):

- `src/beam/osh_beam_model.h:7` → `#include "transport/osh_transport.h"`
  (and it is a *header*, so every beam consumer inherits the back-edge);
- `src/gemca/osh_gemca2_calc_zone.c:10` → same include.

Trivial to remove; the real finding is that nothing enforces the rule. A
20-line CI script asserting the module dependency DAG (grep `#include`
per module, compare against an allowlist) would keep this from regressing.

### A-2 (info, high) Hot-path allocation rule (§10) verified clean {: #a-2 }

No `malloc/calloc/realloc/free` anywhere under `osh_scoring_score_step()`,
`osh_scoring_score_point()` or `osh_transport_ion_step()`; the neutron
xsec state is allocation-free and per-run. Scratch is caller-owned as the
rule demands. (The only per-batch allocation is `osh_neutron_xsec_compile`
once per neutron drain — setup-scale, acceptable.)

### A-3 (medium, high) TODO.md contradicts the code in load-bearing places {: #a-3 }

- The **neutron section's "What is implemented" list is wrong** about the two
  points that matter most for dose: "n-p recoil proton returned" and charged
  secondaries "(currently deposited locally)" — see [N-1](#n-1); neither deposits nor
  feedback exist. The open-items list is accurate; the implemented list will
  mislead anyone triaging validation discrepancies.
- **STRAGG 2 "validated vs SH12A"** is cited with numbers, but the entire
  `reference::idd_*` tier that would lock those numbers in never runs ([T-2](#t-2)).
- "**Lazy extension of material/projectile runtime tables at batch
  boundaries**" (Material section) conflicts with the parallel plan: under
  `--threads N` the material runtime is a shared read-only structure; lazy
  extension mid-run means either per-worker table copies or a
  stop-the-world barrier — this needs to be reconciled with [#161][#161]'s memory
  policies *before* someone implements it as a simple `realloc`.
- The checkpoint/batch scheduler items marked done are consistent with the
  code (`osh_transport.c` outer loop, K = nstat fast path, family-exact dumps,
  completeness labels in the BDO writer) — verified, good.

### A-4 (info, high) Placement and design notes (no action forced) {: #a-4 }

- `osh_checkpoint_policy` / `osh_run_control` living in `src/transport/` is
  defensible (they gate the transport outer loop); when [#161][#161]'s merge step
  lands they will also coordinate scoring — if they start including scoring
  headers beyond the snapshot API, revisit ownership (a thin
  `simulation/`-level driver is the natural home).
- The [#195][#195] quiescence contract is honestly implemented: dumps fire only at
  family-quiescent checkpoint boundaries; `completed_nstat` is exact on all
  paths including clean stops; the BDO completeness label exists.
- The variance scaffold (see [S-5](#s-5)) already matches the right [#169][#169] answer
  (batch-means M2 + Schubert–Gertz pairwise merge). The remaining [#169][#169]
  decisions are: batch weight semantics under weighted particles (couples to
  [S-1](#s-1)), and who allocates `data_var` (extend
  `osh_scoring_accumulator_alloc`, currently impossible).

### A-5 (info) Test architecture {: #a-5 }

Covered as [T-1](#t-1)..[T-3](#t-3) in §0: glob-registered phantom cases, never-running
reference tier, CWD-writing tests. One addition: `tests/unit` auto-discovery
by glob is fine (a stray unit file fails loudly at compile, unlike a stray
case dir which fails at *someone else's* `ctest`).

---

## 8. Executive summary & suggested issue breakdown {: #sec-8 }

### Highest-impact findings (result-correctness)

| # | Finding | Severity | Where |
|---|---|---|---|
| [G-1](#g-1) | Ellipsoid/elliptic-cyl/cone distances wrong (missing ×2 on linear coefficient); empirically confirmed 0.651 vs 0.500 — filed as [#255][#255], fixed by [#257][#257] | **critical** (for ELL/ELLZ/CONE users) | `osh_gemca2_dist.c:413,443,467` |
| [N-1](#n-1) | Neutron interactions deposit zero energy; n–p recoil energy vanishes; [TODO.md][todo-md] claims otherwise | **high** | `osh_transport_neutron.c:166-202`, `osh_neutron_reaction.c:91-108` |
| [N-2](#n-2) | Neutron σ clamped at 20 MeV table edge — fast neutrons ×5–6 over-attenuated | **high** | `osh_neutron_xsec.c:84-87` |
| [P-2](#p-2) | Isotope conflation: deuteron/triton/³He use wrong range table (×2/×3/×¾) and wrong rest mass — filed as [#267][#267], fixed by [`a0f60a1`][a0f60a1] | **high** (resolved) | `osh_transport_ion_step.c:1258`, `osh_material_compile.c:698-716` |
| [P-1](#p-1) | Residual energy deleted at cutoff kill (scales with TCUT0) — filed as [#279][#279], fixed by [#285][#285] | **high** (resolved) | `osh_transport_ion_step.c:539-543` |
| [S-1](#s-1) | `st->wt` ignored by every scorer — latent for MCPL ([#41][#41])/VR | **high (latent)** | `src/scoring/runtime/*` |
| [P-3](#p-3) | Bohr straggling missing relativistic factor; κ-dispatch variance discontinuity — filed as [#280][#280], fixed by [#287][#287] | medium (resolved) | `osh_physics_strag_gauss.c:25-37` |
| [P-4](#p-4) | Nuclear vertex at step endpoint (even past boundary nudge) | medium | `osh_transport_ion_step.c:307-484` |
| [S-2](#s-2) | LET/spectra midpoint energy uses `q[3]=0` on nuclear-kill steps | medium | `osh_scoring_step.c:441` vs `osh_transport_ion_step.c:1092` |
| [S-3](#s-3) | CT dose uses step-start density for all crossings (known TODO, confirmed) | medium | `osh_scoring_step.c:838-841` |

### Parallelization (#161) blockers found beyond the tracked ones

[R-3](#r-3) statics table (`s_proton/s_alpha` lazy init, `s_n_warned`, Molière init
race in the WIP driver) and [R-4](#r-4)'s four concrete bugs in
`feat/parallel-threads:osh_transport_parallel.c` (join-on-uncreated-thread UB,
two leaky error paths, O(threads·nstat) neutron pools, hand-rolled family
scheduling that will silently drop ion-feedback secondaries once [N-1](#n-1) is
fixed).

### Infrastructure

[T-1](#t-1) phantom test dirs, [T-2](#t-2) physics-validation tier never runs (no reference
curves committed), [T-3](#t-3) CWD-writing tests, [A-1](#a-1) dead layering back-edges,
PR [#239][#239] validated and worth merging as-is ([T-4](#t-4)).

### Suggested issue slicing (one issue each, in priority order)

1. GEMCA quadric linear-coefficient bug + per-surface distance unit tests ([G-1](#g-1), [G-2](#g-2)) — filed as [#255][#255], fixed by [#257][#257].
2. Neutron energy deposition: point-deposit `local_deposit_mev` + fix [TODO.md][todo-md] claims ([N-1](#n-1)); follow-up: recoil-proton ion feedback.
3. Neutron σ above 20 MeV: Tier-2 dispatch above Tier-1 grid ([N-2](#n-2)).
4. Isotope-aware range/mass in transport (A-rescale within Z column + `part->mass` kinematics) ([P-2](#p-2)) — filed as [#267][#267], fixed by [`a0f60a1`][a0f60a1].
5. Deposit residual energy at all cutoff/species kills ([P-1](#p-1), and neutron cutoff from [N-1](#n-1)) — filed as [#279][#279], fixed by [#285][#285]. Broader neutron reaction-deposit scoring in [N-1](#n-1) (capture, (n,p)/(n,α), n-p elastic recoil) remains open.
6. Weight-aware scoring + weighted-variance contract ([S-1](#s-1), feeds [#169][#169] and [#41][#41]).
7. Relativistic Bohr straggling factor ([P-3](#p-3)) — filed as [#280][#280], fixed by [#287][#287].
8. Nuclear vertex sampling within the step ([P-4](#p-4)) + step-midpoint energy for scoring on kill steps ([S-2](#s-2)).
9. Reference-test data: commit IDD curves, wire `ctest -L reference` into CI ([T-2](#t-2)); case-dir input validation ([T-1](#t-1)); test tmpdir hygiene ([T-3](#t-3)).
10. Parallel-driver fixes on `feat/parallel-threads` before merge ([R-4](#r-4) + [R-3](#r-3) statics).
11. Housekeeping: dead transport includes ([A-1](#a-1)), dead `accumulator_rescale` ([S-5](#s-5)), HU>1600 clamp warning ([G-4](#g-4)), DICOM robustness nits ([IO-1](#io-1)..[IO-4](#io-4)).

### Method note

Everything above was verified by direct code reading with exact line
references, plus: full debug + ASan/UBSan builds and test runs (only PR [#239][#239]'s
two known test bugs fail), and a standalone numerical harness for [G-1](#g-1)
(`0.6513878189` vs true `0.5`; `b×2` restores correctness). Severities assume
current usage; "latent" items are correct today but break planned features.

<!-- Issue / PR references (reference-style links; also render as clickable #nnn on GitHub). -->
[#41]: https://github.com/openshieldhit/openshieldhit/issues/41
[#154]: https://github.com/openshieldhit/openshieldhit/issues/154
[#158]: https://github.com/openshieldhit/openshieldhit/issues/158
[#161]: https://github.com/openshieldhit/openshieldhit/issues/161
[#165]: https://github.com/openshieldhit/openshieldhit/issues/165
[#166]: https://github.com/openshieldhit/openshieldhit/issues/166
[#168]: https://github.com/openshieldhit/openshieldhit/issues/168
[#169]: https://github.com/openshieldhit/openshieldhit/issues/169
[#178]: https://github.com/openshieldhit/openshieldhit/issues/178
[#179]: https://github.com/openshieldhit/openshieldhit/issues/179
[#193]: https://github.com/openshieldhit/openshieldhit/issues/193
[#195]: https://github.com/openshieldhit/openshieldhit/issues/195
[#212]: https://github.com/openshieldhit/openshieldhit/issues/212
[#222]: https://github.com/openshieldhit/openshieldhit/pull/222
[#235]: https://github.com/openshieldhit/openshieldhit/issues/235
[#236]: https://github.com/openshieldhit/openshieldhit/issues/236
[#237]: https://github.com/openshieldhit/openshieldhit/issues/237
[#239]: https://github.com/openshieldhit/openshieldhit/pull/239
[#248]: https://github.com/openshieldhit/openshieldhit/pull/248
[#255]: https://github.com/openshieldhit/openshieldhit/issues/255
[#257]: https://github.com/openshieldhit/openshieldhit/pull/257
[#267]: https://github.com/openshieldhit/openshieldhit/issues/267
[#279]: https://github.com/openshieldhit/openshieldhit/issues/279
[#280]: https://github.com/openshieldhit/openshieldhit/issues/280
[#285]: https://github.com/openshieldhit/openshieldhit/pull/285
[#287]: https://github.com/openshieldhit/openshieldhit/pull/287
[#299]: https://github.com/openshieldhit/openshieldhit/issues/299
[#301]: https://github.com/openshieldhit/openshieldhit/pull/301
[todo-md]: https://github.com/openshieldhit/openshieldhit/blob/e3c619a1328e9351bcbc1dc599321ac2770ad622/TODO.md
[e3c619a]: https://github.com/openshieldhit/openshieldhit/commit/e3c619a1328e9351bcbc1dc599321ac2770ad622
[a0f60a1]: https://github.com/openshieldhit/openshieldhit/commit/a0f60a1d201ed147719f26751befbd65488ab536
