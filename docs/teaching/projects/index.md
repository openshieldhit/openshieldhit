# Group projects (problem-based learning)

Unlike the [exercise sets](../courses/index.md), which are designed to be
solved without touching source code, the projects on this page are
multi-week, open-ended group assignments that **do** modify OpenShieldHIT —
its physics, its tooling, or its benchmarks. They are sized for graduate or
post-graduate groups of 2–4, and each lists scope tiers so a group can
succeed at tier 1 and stretch toward the later tiers.

Every project follows the same working pattern, which mirrors how the code
itself is developed:

1. **Proposal** (1 page): the chosen tier, the validation data you will
   compare against, and the observable that defines success — *before*
   writing code.
2. **Benchmark first**: build or adopt the reference comparison early; it is
   the ground truth for everything that follows.
3. **Deliverables**: one or more focused pull requests (one logical change
   each, green CI), a reproducible case directory under `tests/` or
   `examples/`, and a short report that separates *what was measured* from
   *what was concluded*.

## Working with LLM assistants

These projects are explicitly designed for groups working alongside LLM
coding agents. The repository is prepared for that: `llms.txt` gives the
project map, `CLAUDE.md` the agent operating rules, and `DEVELOPER.md` the
numbered style rules that CI enforces. Ground rules that experience in this
repository has validated:

- **Benchmarks are the ground truth, not model confidence.** A model will
  argue eloquently for wrong physics; a committed reference comparison will
  not. Decide observables before generating code.
- **Verify every physics formula against the cited primary literature.**
  Models mis-remember constants and exponents; the repository convention is
  to implement from published equations and cite the paper, never to copy
  from other codes (license and provenance both matter — see the note in
  [the nuclear physics reference](../../physics/nuclear-inelastic.md)).
- **Keep the loop evidence-driven.** A good physics-development cycle is:
  measurement stage → model stage → calibration stage, with each decision
  backed by a numbers table.  This is a better mental model for LLM-assisted
  physics development than any style guide.
- **Declare the collaboration.** Reports must state what the LLM produced,
  what the group verified, and how (academic integrity, and it makes the
  work reproducible).

---

## P1 — Photon transport (physics + code)

**Background.** OpenShieldHIT transports ions and neutrons; photons hit a
17-line stub that raises "photon transport is not implemented"
(`src/transport/osh_transport_photon.c`). Meanwhile the physics produces
photons that are silently dropped or deposited locally: neutron-capture
gammas (e.g. the 2.2 MeV line in the BNCT example, `examples/02_bnct`) and
nuclear de-excitation energy. A partial implementation is genuinely useful —
it does not have to be complete.

**Scope tiers.**

1. Photoelectric absorption + Compton scattering (Klein–Nishina sampling),
   cross sections from NIST XCOM tables or a standard parametrisation;
   photon pool modeled on the existing neutron pool.
2. Pair production above 1.022 MeV with annihilation-photon emission;
   coherent (Rayleigh) scattering as a refinement.
3. Close the sources: route neutron-capture and de-excitation gammas into
   the photon pool instead of local deposit, and quantify the effect on the
   BNCT benchmark.

**Validation.** Attenuation coefficients vs NIST XCOM per material; depth
dose of a Co-60 / 6 MV-like spectrum in water vs published data; energy
conservation per event (the repository's standing invariant — see
`tests/unit/test_osh_nuclear_*.c` for the pattern).

---

## P2 — Mining EXFOR: new nuclear benchmarks (nuclear data)

**Background.** The nuclear models are validated against a small set of
references (`tests/reference/`, SH12A and TOPAS fixtures, and the measured
data sets listed in the
[nuclear physics reference](../../physics/nuclear-inelastic.md)). The EXFOR
library holds far more: light-ion production DDX, neutron production, and
activation channels for p + C/N/O/Ca at therapy energies.

**Scope tiers.**

1. Retrieve and digitise one measured data set (e.g. Tippawan 96 MeV n+O
   light-ion DDX, or Bertrand & Peelle proton-induced spectra), build a
   reference case in the `tests/reference/` layout with a comparison script
   (follow `tests/reference/README.md` and `tools/plot_nucre.py`).
2. Add an **activation benchmark**: C-11 / O-15 production along a proton
   beam in PMMA or water (the PET range-verification observable) — a channel
   the current chain has never been scored against.
3. Improvement study: can a calibration knob (`OSH_PREEQ_GAMMA_SCALE_*`,
   `OSH_ABRASION_SIGMA_PN_MB`, the Tripathi σ_R systematics) be tuned to
   improve your new benchmark *without regressing the committed ones*?
   Propose the change with the full before/after table.

**Validation.** Your own new reference cases, plus the requirement that all
existing `reference::` and `cases::` comparisons stay green.

---

## P3 — Computational science: performance, correctness, parallelism

Three sub-projects, each self-contained. The performance harness under
`benchmarks/performance/`, the parser/runtime split, and the existing reference
tests give concrete anchors.

**P3a — Performance engineering.** Profile the transport loop (`perf`,
flame graphs), produce a roofline analysis, then attack a measured hotspot.
Possible directions include tabulating expensive cross sections on hot paths,
evaluating structure-of-arrays layouts for rate evaluation, or assessing
link-time/profile-guided optimisation.  The hard rule: physics results must
remain bit-identical or statistically indistinguishable, and you must define
which of the two you claim and how you test it.

**P3b — Correctness engineering.** The parsers (`src/*/parse/`) are cleanly
separated from the runtime — ideal fuzzing targets (libFuzzer/AFL++). The
nuclear event generators expose exact conservation invariants — ideal for
property-based testing beyond the fixed-seed unit tests. Grade the test
suite with mutation testing; extend sanitizer coverage; triage and fix what
you find. `docs/dev/bug-hunts/` documents previous audits and their method —
continue the series.

**P3c — Parallelism and reproducibility.** Define and implement
worker-parallel reproducibility: RNG streams that are stable per history,
deterministic or well-characterised reductions, and statistical-equivalence
tests built on batch-means standard errors.  Ambitious groups can study
browser/WebAssembly workers or a GPU feasibility prototype, but the CPU
reproducibility contract should come first.

---

## P4 — Medical physics: cavities, heterogeneities, and TPS comparison

**Background.** Proton dose calculation degrades near low-density cavities
and material interfaces, where multiple scattering and straggling models are
stressed — clinically relevant (sinuses, lung, implants) and a classic
model-comparison problem. Reference fixtures from TOPAS already exist under
`tests/reference/topas/`.

**Scope tiers.**

1. Slab phantoms with air cavities and bone inserts: depth-dose and lateral
   profiles vs a TPS-grade Monte Carlo (TOPAS/Geant4 or MCsquare, both
   free for research); gamma-index (2%/2 mm) pass-rate analysis.
2. Systematics: cavity size/depth scan; identify where OpenShieldHIT's
   Molière/Highland MCS and straggling models drive the deviations
   (see [the MCS reference](../../physics/multiple-scattering.md)).
3. Toward patient-like geometry: a simple voxelised phantom — scoping study
   for the planned CT/voxel revival (`src/gemca/voxel/` is present but
   stale).

**Validation.** Cross-code comparison with quantified statistical and model
uncertainty; a written judgement of *which* model differences explain the
observed deviations.

---

## P5 — Delta-ray production and restricted LET (physics)

**Background.** Ion energy loss is currently deposited on-track using
unrestricted stopping power. For microdosimetry, detector response, and any
future electron transport, the energetic knock-on electrons (δ-rays) must be
separated: sample deltas above a threshold Δ from the free-electron
cross section, and deposit only the restricted stopping power L_Δ on-track.
The maximum-energy-transfer machinery (E_max) already exists in the
straggling code (`src/transport/osh_transport_ion_step.c`).

**Scope tiers.**

1. Scorer-only budget study: what fraction of the local dose is carried by
   deltas above 1 / 10 / 100 keV, as a function of depth and primary LET —
   no transport change, one new scorer.
2. Production: sample δ-rays above Δ, bank them (count, spectrum, angular
   distribution), deposit L_Δ on-track; exact energy bookkeeping between
   track, bank, and deposit.
3. Groundwork electron transport: CSDA-range-based electron stepping for the
   banked deltas (a natural companion to P1's photon pool).

**Validation.** L_Δ vs ICRU 37 / PSTAR restricted stopping powers; total
energy conservation with the bank included; comparison of the δ-threshold
sensitivity against microdosimetry literature.

---

## P6 — Track-structure detector response: quenching scorers

**Background.** Real detectors under-respond to high-LET radiation
(scintillator quenching, alanine, TLD, radiochromic film). OpenShieldHIT
already scores dose- and track-averaged LET and (z_eff/β)² (the
`DLET/TLET/DQEFF/TQEFF` estimators in `src/scoring/runtime/`) — the
infrastructure for **response-weighted dose** scorers is in place.

**Scope tiers.**

1. Birks-type scorer: D_response = Σ dose_i · η(LET_i) with the Birks kB as
   the knob; validate against published proton scintillator quenching data.
2. Amorphous-track (Katz-type) relative-effectiveness tables for alanine or
   TLD, as pluggable η(z, β) tables — connecting the scorer to the
   track-structure literature.
3. Detector-comparison study: predict the response of two detector types in
   the same field (e.g. plateau vs distal edge of a proton beam) and compare
   with published measurements; quantify the impact of recoil and fragment LET
   components that are not yet represented in the scorer.

**Validation.** Published quenching/relative-effectiveness data; internal
consistency (η ≡ 1 must reproduce the plain dose scorer exactly).

---

## Further seeds

Smaller or more speculative directions, suitable as tier-1-only projects or
as extensions:

| Seed | Possible direction |
|---|---|
| Particle trajectory dump + interactive web viewer | Export selected tracks and render them in a lightweight browser viewer |
| Variance reduction: splitting / implicit capture study | Use batch-means uncertainty as the comparison metric |
| Statistical-coverage audit of the MC error bars | Check whether reported confidence intervals have the expected coverage |
| Weisskopf evaporation for A > 16 residues | Extend the heavy-residue de-excitation path beyond the current break-up domain |
| Native plots from scorer output (PNG/SVG) | Generate quick-look plots directly from text or BDO scorer output |
