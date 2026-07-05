# Fast nuclear reaction stage — findings and implementation plan

Status: **design study** for issue
[#221](https://github.com/openshieldhit/openshieldhit/issues/221)
(fast nuclear reaction stage: intranuclear cascade + coalescence +
pre-equilibrium, before Fermi break-up).  Nothing in this document is
implemented yet; it records the analysis, the alternatives that were
considered, and the staged plan.

---

## 1. Problem recap and quantitative targets

OSH's inelastic chain for protons is *abrasion → Fermi break-up (FBU)*
(`osh_nuclear_handler.c`, abrasion at `osh_nuclear_abrasion.c`, FBU at
`osh_nuclear_fermi_breakup.c`).  Abrasion emits only fast **nucleons**; every
light cluster (d, t, ³He, α) comes from the statistical break-up of a residue
whose mean excitation is only E\* ≈ 18 MeV.  The result, measured against
SH12A on the 200 MeV p → water NUCRE decks of issue #212:

| Quantity (alphas)                        | OSH        | SH12A      | Target        |
|------------------------------------------|------------|------------|---------------|
| production yield [α / primary]           | ~0.20      | ~0.3 (est.)| ≈ 0.3         |
| production spectrum mean [MeV]           | 3.2        | (harder)   | ≥ ~6          |
| fluence-weighted plateau mean [MeV]      | ~6         | ~10.5      | ~10           |
| fluence ratio OSH/SH12A                  | 0.66×      | 1          | ≈ 1           |
| dose/fluence (≈ mean LET) ratio          | 1.4–1.7×   | 1          | within ~20%   |
| α dose ratio                             | 0.96×      | 1          | keep ≈ 1      |

The α **dose** already matches because the chain conserves energy — fewer,
softer alphas deposit the same total.  What is missing is a **hard component**
of the light-cluster spectrum (a tail reaching several tens of MeV), which
dominates the fluence-weighted mean and therefore the dose-averaged LET.
Per the issue's scope note, approximate is acceptable: the goal is a
qualitatively right spectrum and DLET moving in the right direction, not an
exact match.

## 2. What the reference chain (SH12A / CEM lineage) actually does

The Dubna cascade–exciton lineage (Toneev & Gudima; Gudima, Mashnik & Toneev;
today's CEM03.03 in MCNP6) runs four stages:

1. **Intranuclear cascade (INC)** — the projectile and struck nucleons
   propagate through the nucleus, colliding on parametrised free hh cross
   sections with Pauli blocking; π/Δ above threshold.  Output: escaped fast
   nucleons *with correlated, forward-peaked momenta*, plus a residue with an
   exciton configuration (p particles, h holes) and E\*.
2. **Coalescence** — escaped cascade nucleons whose *relative momenta* fall
   inside a sphere of radius p_c form d, t, ³He, α.  CEM03.03 uses
   p_c(d) = 90, p_c(t/³He) = 108, p_c(α) = 115 MeV/c below 300 MeV incident
   (150/175/175 MeV/c up to 1 GeV).  This is where the *fast/hard* alphas
   come from.
3. **Pre-equilibrium (Modified Exciton Model)** — the (p, h) residue emits
   n, p, d, t, ³He, α while relaxing toward equilibrium; cluster emission uses
   a condensation ("coalescence of excitons") probability γ_j.  This fills the
   mid-energy part of the cluster spectra (from the Coulomb barrier up to
   roughly E\*).
4. **Equilibrium de-excitation** — FBU for light residues (OSH already has
   this), SMM/evaporation/fission for heavier ones (out of scope here, OSH
   targets are light).

Load-bearing observation: **stages 2 and 3 are the alpha levers, and they sit
*on top of* a cascade that supplies (a) several correlated fast nucleons per
event and (b) residues with a broad E\* distribution reaching tens of MeV.**
OSH's abrasion currently supplies neither.

## 3. Why coalescence on the current abrasion output cannot work

This kills the "first experiment" candidate from the issue (coalescence on the
existing abrasion output), on multiplicity grounds alone:

- Wounded-nucleon mean for 200 MeV p + ¹⁶O:
  ⟨ν⟩ = σ_pN·A/σ_R ≈ 30 mb × 16 / ≈300 mb ≈ **1.6** participants; the typical
  event emits 1–3 nucleons (knockouts + the degraded cascade proton).
- An α needs **≥ 4 escaped nucleons (2p + 2n)** in one event.  With
  ν ~ zero-truncated Poisson(1.6), P(ν ≥ 4) ≈ 10%, and the 2p2n composition
  requirement cuts that further.
- Worse, abrasion's knockouts are sampled with *isotropic CM* angles and tens
  of MeV of kinetic energy; typical pairwise relative momenta are several
  hundred MeV/c.  The coalescence acceptance (p_c ≈ 115 MeV/c for α — a
  4-body momentum-space coincidence) is then vanishingly small.  Both
  effects compound: the expected α coalescence yield is effectively zero.

Coalescence only pays off after a *real* multi-generation cascade exists
(higher multiplicities, momentum-correlated forward emission).  That is the
expensive road (Section 4).  Conclusion: **for a lean OSH, the hard-alpha
deficit must be fixed by direct/pre-equilibrium cluster *emission* channels,
not by coalescence of abrasion nucleons.**

## 4. Reference implementations surveyed (Geant4 tree)

Line counts (`.cc` + `.hh`) of the candidate references, from the Geant4
source attached to this session:

| Model                                   | LOC     | Verdict for OSH |
|-----------------------------------------|---------|-----------------|
| INCL++ (cascade + surface coalescence)  | ~68 500 | far too large; needs per-nucleon position+momentum tracking |
| Bertini cascade (`cascade/cascade`)     | ~58 700 | too large |
| Binary cascade                          | ~9 900  | too large, C++-idiomatic |
| **Exciton pre-equilibrium** (`pre_equilibrium/exciton_model`) | ~6 700 | **tractable**: the physics core (transition rates + emission widths) is a few hundred lines; the rest is class boilerplate |
| G4 Fermi break-up                       | ~4 200  | OSH already has its own |

Useful extractions from the exciton model (all verified in source):

- `G4PreCompoundTransitions.cc` — the CEM-flavoured Δn = +2/0/−2 transition
  rates (Pauli-corrected ⟨σ(v_rel)·v_rel⟩/V_int form) are ~80 lines of math.
- `G4PreCompoundIon.cc` + per-fragment classes — cluster emission widths with
  *closed-form* condensation factors (e.g. γ_α = 4096/A³, plus simple
  combinatorial R_j and factorial factors), Coulomb barriers, and inverse
  cross sections.  No data tables required beyond masses OSH already has.
- INCL++'s clustering (`G4INCLClusteringModelIntercomparison.cc`) confirms
  that a faithful dynamical coalescence needs full phase-space (r, p) nucleon
  tracking — not available in abrasion, and not worth adding for this.

## 5. Candidate approaches

### A. Exciton pre-equilibrium module with cluster emission (mechanistic core)

A single-component exciton model in plain C, slotted **between abrasion and
FBU**:

- State: (p, h, E\*) of the residue; single-particle level density
  g = (6/π²)·a with a ≈ A/13 MeV⁻¹ (one calibration constant).
- Per iteration, compete: emission widths Γ_j(e) for j ∈ {n, p, d, t, ³He, α}
  (CEM/G4-style: e·σ_inv(e)·γ_j·R_j·[exciton phase-space ratio]) against the
  Δn = +2 transition rate λ₊; sample, emit or absorb a 2-exciton pair; stop at
  the equilibrium exciton number n_eq = √(2 g E\*) (≈ 6–7 for E\* ≈ 30 MeV on
  O) or when no channel is open; hand the equilibrated residue to FBU.
- Angular distributions: Kalbach continuum systematics (two-parameter
  forward-peaked form) — ~30 lines, reusable for abrasion knockouts too.
- Allocation-free at runtime (fixed-size state, bounded loop); σ_inv can be
  precomputed at handler-compile time if profiling demands it.

Cost: **~500–700 lines** of new C + tests.  This is the same physics whose
absence the issue diagnoses ("no pre-equilibrium emission"), in its cheapest
credible form.  It produces d/t/³He/α from the Coulomb barrier up to ~E\*, so
its reach depends on the E\* supply (→ approach A′).

### A′. Cascade-nucleon retention: feed the exciton model (prerequisite, tiny)

Abrasion currently lets *every* knockout escape; holes-only residues
(p = 0) have no exciton particles to emit and E\* stays pinned near
ν × 13.3 MeV.  Add a retention step: each knockout nucleon is absorbed into
the residue with probability P_ret(e) (large at low e, → transparency at high
e), booking its kinetic energy into E\* and counting it as a particle exciton.
Export (p, h) alongside E\* in `struct osh_nuclear_fragment`.

This raises and broadens the E\* distribution (mean 18 → ~30+ MeV, with a
tail), which simultaneously (i) opens the pre-equilibrium cluster channels
and (ii) hardens FBU's own output.  Energy bookkeeping stays exact.
Cost: **~50–100 lines** inside `osh_nuclear_abrasion.c` + one struct field.
One knob (P_ret shape/scale).

### B. Quasi-free cluster knockout branch (the hard-tail channel)

A phenomenological direct channel exploiting the α-cluster structure of the
therapy-relevant light targets (¹²C, ¹⁶O ≈ 3–4 preformed alphas): with
probability P_ko per inelastic event on an α-conjugate target, one abrasion
collision is replaced by **quasi-elastic p + α scattering** (existing
`osh_kinematics_elastic_lab` handles unequal masses), paying the α separation
energy (7.16 MeV for ¹⁶O → ¹²C, 7.37 MeV for ¹²C → ⁸Be) and leaving the
(A−4, Z−2) residue with a hole excitation.  Kinematic ceiling: the α can carry
up to 4·m_p·m_α/(m_p+m_α)² ≈ 0.64 of the proton energy — ~120 MeV at 200 MeV,
i.e. exactly the hard tail OSH lacks.  Optional refinement: smear with cluster
Fermi motion.

This is the *effective stand-in* for what coalescence provides in SH12A, with
physics grounding in the knockout term of Kalbach's pre-equilibrium
phenomenology for complex particles and in (p,pα) quasi-free knockout data on
light nuclei.  Cost: **~150–250 lines**, one dominant knob (P_ko, possibly
energy-dependent), calibrated against the SH12A spectra / evaluated DDX.

### C. Momentum-space coalescence (Toneev–Gudima p_c criterion) — defer

Cheap to *run* (pairwise momentum checks over ≤ ~8 event nucleons, p_c radii
from CEM03.03), but per Section 3 it has nothing to chew on until OSH has a
multi-generation cascade with realistic (forward-peaked, Pauli-blocked)
kinematics.  Revisit only if/when the fast stage is upgraded for ion
projectiles or higher energies; then it drops in at ~200–300 lines behind the
same event interface.

### D. Data-driven inclusive sampling (the "completely alternative approach")

Skip reaction mechanics entirely for protons: per (target, incident energy),
sample ejectile multiplicities, energies, and angles for {p, n, d, t, ³He, α}
from **evaluated inclusive spectra**, with per-event energy capping.  This is
exactly what fast proton-therapy MCs do — MCsquare samples inelastic final
states from ICRU 63 differential data and *transports the alphas*; VMCpro
(Fippel & Soukup) and FRED use cruder parametrisations (alphas local).
A modern open data source is the TENDL proton sub-library (TALYS-evaluated,
CC-licensed, emission spectra with Kalbach angular parameters up to 200 MeV);
tables for the handful of tissue/phantom targets are tiny (~tens of kB) and
would be generated offline by a `tools/` script.

Pros: matches evaluated spectra (including the hard tails) by construction;
smallest *physics* risk; ~300 lines of runtime code.
Cons: proton projectiles on tabulated targets only (no path to carbon beams,
which are core SHIELD-HIT domain); inclusive data break OSH's event-level
exact energy conservation (only approximate per-event balance is possible);
adds a nuclear-data-file dependency and generation tooling; needs
extrapolation 200 → 250 MeV.

**Verdict: not as the runtime model, but adopt its data as the calibration
target.**  The same TENDL/ICRU-63 spectra that option D would sample are the
ideal reference curves for tuning the 3–4 knobs of A′/A/B.

### E. Rejected outright

- **Port/reimplement a full INC** (Bertini, INCL, ISABEL): 10–70 kLOC-class,
  months of validation — against the lean mandate.
- **QMD / BUU class models**: heavier still.
- **ML surrogate generators**: no lean C inference path, opaque energy
  conservation, large validation burden — not competitive with a ~1 kLOC
  mechanistic chain at these energies.

## 6. Recommendation

Implement the **mechanistic lean chain A′ → A → B** (retention → exciton
pre-equilibrium with cluster emission → quasi-free α knockout), calibrated
against SH12A decks and TENDL/ICRU-63 spectra (D-as-data).  Rationale:

1. It is the smallest change set that adds the two missing *mechanisms*
   (pre-equilibrium emission, a direct hard-cluster channel) named in the
   issue's diagnosis.
2. It preserves OSH invariants: exact per-event energy/momentum bookkeeping,
   no runtime nuclear-data tables, no hot-path allocation, ~1 kLOC total.
3. It generalises: the exciton module is projectile-agnostic (it de-excites
   *residues*, so it will serve ion projectiles and already-planned SMM work),
   and wiring it through a common de-excitation entry point finally unifies
   the ion path with the neutron path's `osh_nuclear_compound_step()` adapter
   (currently not wired into the ion handler — noted in the issue).
4. Fallback is cheap: if calibration shows the mechanistic chain cannot reach
   the spectral targets, option D can replace the cluster channels for protons
   behind the same `osh_nuclear_event` interface without touching transport.

## 7. Staged plan

Each stage is one PR-sized logical change, independently testable, with a
go/no-go metric before the next.

### Stage 0 — Benchmark harness + E\*-ceiling experiment (no new physics)

- Reduce the #212 NUCRE decks to a repeatable local benchmark: tally α (and
  d/t/³He) production spectra at emission, yields/primary, fluence-weighted
  plateau spectra, DLET; script the OSH-vs-SH12A comparison plots.
- Cheap knob experiment: crank `OSH_ABRASION_EXCITATION_PER_HOLE_MEV` (and/or
  disable the cascade-energy clamp) to measure how far E\* alone can push the
  FBU alpha spectrum.  This quantifies the split between "E\* under-supply"
  and "missing direct channel" **before any new code is written**.
- Deliverable: baseline report + acceptance thresholds for later stages.

### Stage 1 — Nucleon retention and (p, h) bookkeeping in abrasion (A′)

- Retention probability P_ret(e) for knockout nucleons; retained energy → E\*,
  retained nucleon → particle exciton; export (p, h) in the fragment struct
  (2 small fields), keep exact energy balance.
- Re-run Stage-0 harness: expect higher/broader E\*, harder FBU alphas, α dose
  ratio still ≈ 1.  Calibrate P_ret jointly with the per-hole constant so the
  *nucleon* spectra (already validated) do not regress.

### Stage 2 — Exciton pre-equilibrium module (A)

- New `src/physics/nuclear/osh_nuclear_preeq.{c,h}`: single-component exciton
  model as specified in Section 5A; emits into the existing
  `osh_nuclear_event` (whitelist {n, p, d, t, ³He, α} — same as FBU);
  equilibrated residue routed to FBU via the compound adapter
  (`osh_nuclear_compound_step()`), unifying ion and neutron paths.
- Unit tests: energy/momentum conservation per event; width sanity vs
  hand-computed cases; equilibrium-limit behaviour (preeq off ⇒ chain reduces
  to current abrasion → FBU).
- Go/no-go: mid-energy (5–20 MeV) cluster yield appears; α production mean
  moves ≥ halfway to target; α dose ratio stays within ~10% of 1
  (pre-equilibrium consumes E\* that used to feed FBU — the joint calibration
  in Stage 4 owns this trade-off).

### Stage 3 — Quasi-free α-knockout branch (B)

- One extra branch in the abrasion collision loop for α-conjugate targets
  (flagged per element at handler-compile time), p + α quasi-elastic
  kinematics + separation energy + residue update; knob P_ko.
- Go/no-go: hard tail (> 20 MeV) appears with roughly the SH12A slope;
  fluence-weighted plateau mean reaches ~10 MeV.

### Stage 4 — Joint calibration + validation

- Tune {P_ret, g (via a), γ-scale, P_ko} against: #212 SH12A spectra, TENDL
  p+¹²C/¹⁶O emission spectra (generated offline into `tests/fixtures/` or a
  benchmark folder by a small `tools/` script), and EXFOR α DDX on C/O where
  available (e.g. Bertrand–Peelle 30–62 MeV sets for the low-energy anchor).
- Acceptance: table in Section 1.  Also verify: total inelastic multiplicity
  and nucleon spectra unregressed; DLET vs depth moves toward SH12A; runtime
  cost of the nuclear step < a few % of total (inelastic events are rare;
  expected impact ≈ none).

### Later (separate issues)

- Kalbach angular systematics for abrasion knockouts and preeq emission
  (drop-in refinement of the isotropic-CM placeholder; helps plateau fluence).
- π/Δ production ≥ 280 MeV (removes the H-inelastic ABSORB fallback).
- SMM / evaporation for A > 16 residues (#176) — the exciton module's
  equilibrium exit is exactly where they plug in.
- Coalescence (C) + a real multi-generation cascade — only with ion-projectile
  fast-stage work; keep the p_c radii and the event-interface slot in mind.

## 8. Risks and open questions

- **E\* budget double-counting**: retention (Stage 1) and the per-hole
  constant both feed E\*; calibrate jointly, don't stack defaults.
- **Preeq may still undershoot the > 30 MeV tail** (E\*-limited): that is what
  Stage 3 is for; the Stage-0 experiment sizes this gap up front.
- **RNG stream changes**: every stage adds draws; fixed-seed regression
  baselines will shift — plan test updates per stage (statistical asserts,
  not bitwise).
- **SH12A production-spectrum reference**: the 0.3 α/primary SH12A yield in
  Section 1 is inferred (fluence 0.66× + dose 0.96×), not directly tallied;
  Stage 0 should pin it down from the decks.
- **Isospin bookkeeping in preeq**: single-component exciton models need the
  charge ratio R_j fudge (G4's form is fine); a two-component model is not
  worth the size at OSH's accuracy bar.

## 9. References

- V. D. Toneev, K. K. Gudima, *Particle emission in light and heavy-ion
  reactions*, Nucl. Phys. A **400** (1983) 173c — cascade + coalescence
  (p_c criterion).
- K. K. Gudima, S. G. Mashnik, V. D. Toneev, *Cascade-exciton model of nuclear
  reactions*, Nucl. Phys. A **401** (1983) 329 — Modified Exciton Model
  (pre-equilibrium).
- S. G. Mashnik et al., *CEM03.03 and LAQGSM03.03 event generators*,
  LA-UR-08-2931 / arXiv:0805.0751 — coalescence radii, MEM details as used in
  MCNP6.
- C. Kalbach, *Systematics of continuum angular distributions*, Phys. Rev. C
  **37** (1988) 2350 — angular systematics.
- C. Kalbach, *Preequilibrium reactions with complex particle channels*,
  Phys. Rev. C **71** (2005) 034606 — pickup/knockout terms for cluster
  emission (grounding for Stage 3).
- A. Iwamoto, K. Harada, Phys. Rev. C **26** (1982) 1821 — exciton-coalescence
  (condensation) picture for cluster pre-equilibrium emission.
- A. J. Koning, M. C. Duijvestijn, Nucl. Phys. A **744** (2004) 15 —
  two-component exciton model, global parametrisations (TALYS).
- A. Boudard, J. Cugnon, J.-C. David, S. Leray, D. Mancusi, Phys. Rev. C
  **87** (2013) 014606 — INCL4.6 incl. surface coalescence for light clusters
  (the faithful-but-large road).
- M. Fippel, M. Soukup, Med. Phys. **31** (2004) 2263 — simplified nuclear
  interactions for fast proton MC (VMCpro).
- K. Souris, J. A. Lee, E. Sterpin, Med. Phys. **43** (2016) 1700 — MCsquare:
  inelastic final states sampled from ICRU 63 data, alphas transported.
- ICRU Report 63, *Nuclear Data for Neutron and Proton Radiotherapy and for
  Radiation Protection* (2000) — evaluated DDX for p on C/N/O up to 250 MeV.
- A. J. Koning, D. Rochman et al., *TENDL: Complete Nuclear Data Library*,
  Nucl. Data Sheets **155** (2019) 1 — open evaluated spectra (calibration
  data source).
- J. P. Bondorf et al., Phys. Rep. **257** (1995) 133 — SMM (future, A > 16).
- Geant4 references used for sizing/extraction: `pre_equilibrium/exciton_model`
  (transition rates, cluster emission widths, condensation factors),
  `inclxx` (clustering algorithm), `cascade/cascade` (Bertini) — see the
  Geant4 source tree.
