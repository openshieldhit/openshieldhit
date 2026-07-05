# Fast nuclear reaction stage — findings and implementation plan

Status: **design study** for issue
[#221](https://github.com/openshieldhit/openshieldhit/issues/221)
(fast nuclear reaction stage: intranuclear cascade + coalescence +
pre-equilibrium, before Fermi break-up).  Nothing in this document is
implemented yet; it records the analysis, a cross-check against the major
reaction-model lineages (CEM/MCNP6, INCL+ABLA, FLUKA/PEANUT, Geant4,
TALYS) and fast therapy MCs (MCsquare, VMCpro, FRED), the alternatives that
were considered, and the staged plan.  Key numbers and formulas were
verified against the primary sources (Section 10).

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
For context on why this matters clinically: secondary particles contribute
little dose in a proton beam but raise LET_d by ~50% in the plateau and by
far more in the penumbra (Grassberger & Paganetti 2011), so getting the
light-cluster spectrum qualitatively right is what moves DLET the right way.
Per the issue's scope note, approximate is acceptable.

## 2. What the reference chain (CEM / Dubna lineage) actually does — and a correction to the issue's emphasis

The cascade–exciton lineage (Toneev & Gudima; Gudima, Mashnik & Toneev;
today's CEM03.03 in MCNP6 — the model family SH12A descends from
conceptually) runs four stages:

1. **Intranuclear cascade (INC)** — dynamic hh collisions on parametrised
   free cross sections, Pauli blocking, π/Δ above threshold.  Output: escaped
   fast nucleons with correlated forward-peaked momenta, plus a residue with
   a **broad distribution** of exciton configuration (p, h), E\*, momentum and
   angular momentum (CEM03.03 manual, §2: the initial pre-equilibrium
   configuration "differs significantly from that usually postulated in
   exciton models").
2. **Coalescence** — escaped cascade nucleons whose *relative momenta* fall
   inside a sphere of radius p_c form d, t, ³He, α.  Verified values
   (CEM03.03 manual Eq. 19): **p_c(d) = 90, p_c(t) = p_c(³He) = 108,
   p_c(⁴He) = 115 MeV/c** for T₀ < 300 MeV (and ≥ 1 GeV);
   150/175/175 MeV/c for 300 MeV < T₀ ≤ 1 GeV.  Coalesced nucleons are
   removed from the nucleon tallies (no double counting).
3. **Pre-equilibrium (Modified Exciton Model, MEM)** — the (p, h) residue
   emits n, p, d, t, ³He, α while relaxing toward equilibrium; cluster
   emission via a condensation probability γ_j (Section 5A has the verified
   formulas).
4. **Equilibrium de-excitation** — FBU for light residues (OSH already has
   this), evaporation/fission or SMM for heavier ones.

**Correction to the issue's emphasis.**  The issue attributes SH12A's fast
alphas primarily to coalescence.  The CEM03.03 manual says otherwise for the
energy range OSH cares about: for nucleon-induced reactions at ~100 MeV on a
medium nucleus (96 MeV n + ⁵⁶Fe, manual Fig. 17) *"the contribution from
coalescence to the total angle-integrated energy spectra of complex particles
is very low, less than a few percent"*, because the cascade-nucleon
multiplicity is small; coalescence only dominates for heavy-ion and GeV-range
reactions.  The manual states repeatedly that CEM produces d/t/³He/⁴He at
63–480 MeV incident energy **"mainly via preequilibrium emission"** (shown
against Guertin 62.9 MeV p+Pb, Cowley 160 MeV p+Al/Co/Au α spectra, Green &
Korteling 210–480 MeV p+Ag).  Its acknowledged deficiency is the *opposite*
end: the very-forward high-energy tails, where the authors "expect a
contribution from direct processes like pick-up and knock-out, **not
considered so far in our models**".

So the load-bearing physics for OSH's alpha deficit at ≤ 250 MeV is
**pre-equilibrium cluster emission** (+ a small direct knockout component for
the extreme tail) — *not* coalescence.  This is good news: the exciton model
is by far the cheapest of the three mechanisms to implement.

## 3. Why coalescence on the current abrasion output cannot work

Independent of Section 2, the "first experiment" candidate from the issue
(coalescence on the existing abrasion output) fails on multiplicity grounds:

- Wounded-nucleon mean for 200 MeV p + ¹⁶O:
  ⟨ν⟩ = σ_pN·A/σ_R ≈ 30 mb × 16 / ≈300 mb ≈ **1.6** participants; the typical
  event emits 1–3 nucleons (knockouts + the degraded cascade proton).
- An α needs **≥ 4 escaped nucleons (2p + 2n)** in one event.  With
  ν ~ zero-truncated Poisson(1.6), P(ν ≥ 4) ≈ 10%, before the composition cut.
- Abrasion's knockouts are sampled with *isotropic CM* angles and tens of MeV
  of kinetic energy; typical pairwise relative momenta are several hundred
  MeV/c against a coalescence acceptance of 90–115 MeV/c (a 4-body
  momentum-space coincidence for α).  The expected α coalescence yield is
  effectively zero.

This is the same low-multiplicity regime in which even CEM's full cascade
yields only a few percent of its clusters from coalescence (Section 2).
Conclusion unchanged and now doubly supported: **fix the alpha deficit with
pre-equilibrium/direct cluster *emission*, not with coalescence of abrasion
nucleons.**

## 4. Cross-model survey: how each code family makes fast light clusters

| Code / lineage | Fast stage | Fast-cluster mechanism at ≤ 250 MeV | Size / notes |
|---|---|---|---|
| CEM03.03 (MCNP6); SH12A conceptually | Dubna INC | **MEM pre-equilibrium emission** (γ_j·M_j condensation); coalescence few-% at these energies; no pickup/knockout (admitted gap at forward tails) | INC+MEM+GEM ~10⁴ lines Fortran |
| INCL4.6 / INCL++ (+ABLA07); PHITS option; Geant4 option | full INC with r–p phase-space tracking | **surface coalescence**: leading nucleon at the surface drags phase-space-close nucleons (h = 1 fm; h₀ = 424/300/300 fm·MeV/c for A=2/3/4); at 63 MeV p+Pb clusters incl. α come *predominantly from this cascade channel* | ~68 kLOC (C++); needs per-nucleon (r, p) — not lean |
| FLUKA / PEANUT | GINC | coalescence "all along the reaction chain" + exciton pre-equilibrium; coalescence credited with the high-energy fragment tails | closed source; same slot filled by two mechanisms |
| Geant4 default therapy lists (BIC + G4PreCompound + de-excitation) | Binary cascade | **exciton pre-equilibrium cluster emission only** (closed-form condensation factors; *no coalescence stage at all*) + FBU/evaporation | precompound module ~6.7 kLOC, physics core a few hundred lines |
| TALYS (evaluation code, source of TENDL) | none (direct+preeq+compound) | two-component exciton model **plus Kalbach's phenomenological pickup/stripping/knockout** for complex particles ("not included in the exciton model"); known to sometimes overestimate α yields | not a transport code; its spectra are the de-facto evaluated reference |
| MCsquare | none | samples inelastic final states (incl. α, transported) **directly from ICRU 63 data** | proof that inclusive-data sampling suffices for proton therapy dose/LET |
| VMCpro / FRED (GPU) | none | parametrised secondaries; α local deposit | cruder than OSH needs for DLET |

Reading across the table: every mechanistic code fills the "fast cluster"
slot with *some* combination of (i) exciton-model condensation emission and
(ii) a direct surface/knockout-like channel, and the evaluation community
(TALYS/Kalbach) uses exactly (i)+(ii) in phenomenological form.  The
175 MeV quasi-monoenergetic n + O/Si/Fe/Bi measurement (Uppsala TSL) makes
the same point experimentally — the authors single out *"the formation and
emission of composite particles in the pre-equilibrium stage"* as the most
crucial aspect separating INCL4.5-ABLA07, MCNP6, TALYS and PHITS variants.

## 5. Candidate approaches

### A. Exciton pre-equilibrium module with cluster emission (mechanistic core)

A single-component exciton model in plain C, slotted **between abrasion and
FBU**, following the MEM (verified formulas, CEM03.03 manual Eqs. 20–34) and
its Geant4 sibling (`G4PreCompound*`):

- State: (p, h, E\*) of the residue; equidistant single-particle level
  density g, exciton state density (Williams):
  ω(p, h, E) = g·(gE)^{p+h−1} / (p!·h!·(p+h−1)!).
- Emission width for ejectile j (detailed balance):
  λ_c^j(p,h,E,T) = (2s_j+1)/(π²ħ³) · μ_j · ℜ_j(p,h) ·
  [ω(p−p_j, h, E−B_j−T)/ω(p,h,E)] · T·σ_inv(T), integrated from the Coulomb
  barrier V_j to E−B_j; σ_inv from Dostrovsky-style parametrisation.
- Cluster emission: p_j excitons condense with probability
  **γ_j ≃ p_j³·(p_j/A)^{p_j−1}**, *times an empirical scale factor* —
  CEM03.03 multiplies this "rather crude estimate" by fitted coefficients
  M_j(A, Z, T₀); Geant4 uses closed-form equivalents (e.g. γ_α = 4096/A³)
  with its own R_j charge-composition factors.  Lesson for OSH: ship γ_j with
  **one calibration scale per ejectile species**, expect to fit it (Stage 4).
- Transition rates λ₊ (Δn = +2) from the CEM form
  (⟨σ(v_rel)·v_rel⟩/V_int with a Pauli factor — ~80 lines in
  `G4PreCompoundTransitions.cc`); λ₋/λ₀ optional at OSH's accuracy bar
  (Δn = +2 dominates far from equilibrium, which is where emission matters).
- Equilibrium exit at **n_eq = √(2gE)** (CEM03.03 Eq. 34), then hand the
  residue to FBU (A ≤ 16) via the compound adapter.
- Angular distributions: **Kalbach continuum systematics** — this is not a
  refinement but the standard: CEM03.03 itself uses Kalbach's systematics for
  nucleons *and* complex particles below 210 MeV incident because it beats
  their own exciton-based angular model.  ~30 lines, reusable for abrasion
  knockouts.
- Allocation-free at runtime (fixed-size state, bounded loop); σ_inv
  precomputable at handler-compile time if profiling demands.

Cost: **~500–700 lines** of new C + tests.  Per Sections 2 and 4 this is the
*primary* lever, not a supplement: it is the mechanism by which the CEM
lineage — and Geant4's default therapy configuration, which has no
coalescence at all — produces fast light clusters at these energies.

### A′. Cascade-nucleon retention: feed the exciton model (prerequisite, tiny)

Abrasion currently lets *every* knockout escape; holes-only residues (p = 0)
have no exciton particles to emit and E\* stays pinned near ν × 13.3 MeV.
Add a retention step: each knockout nucleon is absorbed into the residue with
probability P_ret(e), booking its kinetic energy into E\* and counting it as
a particle exciton; export (p, h) alongside E\* in
`struct osh_nuclear_fragment`.

Precedent: INCL4.6's "back to spectator" recipe — a cascade nucleon whose
energy falls below the Fermi level plus ξ ≈ 7 MeV (plus ⅔ of the Coulomb
barrier for protons) is re-absorbed into the remnant.  A kinetic-energy
threshold with a smooth turn-on is the lean equivalent.  This raises and
broadens the E\* distribution (mean 18 → ~30+ MeV with a tail), which
simultaneously (i) opens the pre-equilibrium cluster channels and
(ii) hardens FBU's own output.  Energy bookkeeping stays exact.
Cost: **~50–100 lines** inside `osh_nuclear_abrasion.c` + one struct field.
One knob (P_ret shape/scale), calibrated jointly with the per-hole constant.

### B. Quasi-free cluster knockout branch (the hard-tail channel)

A phenomenological direct channel exploiting the α-cluster structure of the
therapy-relevant light targets (¹²C, ¹⁶O): with probability P_ko per
inelastic event on an α-conjugate target, one abrasion collision is replaced
by **quasi-elastic p + α scattering** (existing
`osh_kinematics_elastic_lab` handles unequal masses), paying the α separation
energy (7.16 MeV for ¹⁶O → ¹²C, 7.37 MeV for ¹²C → ⁸Be) and leaving the
(A−4, Z−2) residue with a hole excitation.  Kinematic ceiling: the α can
carry up to 4·m_p·m_α/(m_p+m_α)² ≈ 0.64 of the proton energy — ~120 MeV at
200 MeV, i.e. exactly the extreme tail that even CEM admits to missing
("pick-up and knock-out, not considered so far in our models") and that
TALYS covers with Kalbach's phenomenological knockout term (which treats α
knockout as significant precisely for light cluster-structured targets).
INCL's surface coalescence is the other lineage's answer to the same gap.
Optional refinement: smear with cluster Fermi motion.

Cost: **~150–250 lines**, one dominant knob (P_ko, possibly energy-dependent).
Caution from the TALYS experience: phenomenological cluster channels
overshoot easily — calibrate against spectra, don't stack defaults.

### C. Momentum-space coalescence (Toneev–Gudima p_c criterion) — defer indefinitely for protons

Cheap to *run* (pairwise momentum checks, verified p_c radii above), but per
Sections 2–3 it is a few-percent effect at proton-therapy energies even atop
a full cascade, and OSH's abrasion gives it nothing to chew on.  Revisit only
with an ion-projectile fast stage or > 300 MeV physics, where multiplicities
grow and coalescence becomes the dominant hard-cluster source; it then drops
in behind the same event interface at ~200–300 lines (remember to remove
coalesced nucleons from the nucleon output, as CEM does).

### D. Data-driven inclusive sampling (the "completely alternative approach")

Skip reaction mechanics entirely for protons: per (target, incident energy),
sample ejectile multiplicities, energies, and angles for {p, n, d, t, ³He, α}
from **evaluated inclusive spectra**, with per-event energy capping.  This is
exactly what MCsquare does (inelastic final states from ICRU 63 differential
data, alphas transported) — with clinically validated dose/LET.  A modern
open data source is TENDL's proton sub-library: explicit channels up to
30 MeV, continuum emission spectra (MF6/MT5, Kalbach-Mann representation,
i.e. angular shape comes for free) up to **200 MeV**, for every stable
target; tables for the handful of tissue/phantom nuclides are tens of kB and
would be generated offline by a `tools/` script.

Pros: matches evaluated spectra (including hard tails) by construction;
smallest physics risk; ~300 lines of runtime code.
Cons: proton projectiles on tabulated targets only (no path to carbon beams —
core SHIELD-HIT domain); inclusive data break OSH's event-level exact energy
conservation (CEM03.03 itself moved *to* exact per-event conservation, from
on-average); a nuclear-data-file dependency plus generation tooling; TENDL
tops out at 200 MeV (extrapolation needed to 250); TALYS-based evaluations
inherit TALYS's α-overestimation tendency on some targets.

**Verdict: not as the runtime model, but adopt its data as calibration
targets.**  The same TENDL/ICRU-63 spectra that option D would sample are the
reference curves for tuning the knobs of A′/A/B (together with the measured
DDX sets in Section 8).

### E. Rejected outright

- **Port/reimplement a full INC** (Bertini ~59 kLOC, INCL++ ~68 kLOC,
  ISABEL): months of validation — against the lean mandate; and Section 2
  shows the INC itself is not what makes the alphas at these energies.
- **QMD / BUU class models**: heavier still; PHITS needed a *modified* QMD
  with an added surface-coalescence model to describe the 175 MeV n + O
  composite data at all.
- **ML surrogate generators**: no lean C inference path, opaque energy
  conservation, large validation burden — not competitive with a ~1 kLOC
  mechanistic chain at these energies.

## 6. Recommendation

Implement the **mechanistic lean chain A′ → A → B** (retention → exciton
pre-equilibrium with cluster emission → quasi-free α knockout), calibrated
against SH12A decks, measured n/p + C/O light-ion DDX (Section 8), and
TENDL/ICRU-63 spectra (D-as-data).  Rationale:

1. The literature pass upgraded confidence in exactly this split: the CEM
   lineage makes its fast clusters at ≤ 250 MeV via **pre-equilibrium
   emission** (A), needs an empirical γ-scale fit (built into the plan), uses
   **Kalbach angular systematics** (in the plan), and lacks only the
   **knockout tail** (B).  Geant4's default therapy chain ships (A)-only and
   is the most widely used therapy MC; TALYS ships (A)+(B) phenomenology and
   feeds the evaluated libraries.
2. It preserves OSH invariants: exact per-event energy/momentum bookkeeping
   (CEM03.03 itself converged there), no runtime nuclear-data tables, no
   hot-path allocation, ~1 kLOC total.
3. It generalises: the exciton module de-excites *residues*, so it serves ion
   projectiles and the planned SMM work; wiring it through a common
   de-excitation entry finally unifies the ion path with the neutron path's
   `osh_nuclear_compound_step()` adapter (currently not wired into the ion
   handler).
4. Fallback is cheap: if calibration shows the chain cannot reach the
   spectral targets, option D can replace the cluster channels for protons
   behind the same `osh_nuclear_event` interface without touching transport.

## 7. Staged plan

Each stage is one PR-sized logical change, independently testable, with a
go/no-go metric before the next.

### Stage 0 — Benchmark harness + E\*-ceiling experiment (no new physics)

- Reduce the #212 NUCRE decks to a repeatable local benchmark: tally α (and
  d/t/³He) production spectra at emission, yields/primary, fluence-weighted
  plateau spectra, DLET; script the OSH-vs-SH12A comparison plots.
- Add digitised reference spectra from the measured DDX sets (Section 8) to
  the harness (angle-integrated energy spectra are enough).
- Cheap knob experiment: crank `OSH_ABRASION_EXCITATION_PER_HOLE_MEV` (and/or
  disable the cascade-energy clamp) to measure how far E\* alone can push the
  FBU alpha spectrum.  This quantifies the split between "E\* under-supply"
  and "missing emission mechanism" **before any new code is written**.
- Deliverable: baseline report + acceptance thresholds for later stages.

### Stage 1 — Nucleon retention and (p, h) bookkeeping in abrasion (A′)

- Retention probability P_ret(e) for knockout nucleons (INCL-style
  low-energy re-absorption); retained energy → E\*, retained nucleon →
  particle exciton; export (p, h) in the fragment struct; keep exact energy
  balance.
- Re-run Stage-0 harness: expect higher/broader E\*, harder FBU alphas, α
  dose ratio still ≈ 1.  Calibrate P_ret jointly with the per-hole constant
  so the *nucleon* spectra (already validated) do not regress — watch the
  20–50 MeV nucleon region in particular; INCL reports that adding cluster
  channels "eats" nucleons exactly there.

### Stage 2 — Exciton pre-equilibrium module (A)

- New `src/physics/nuclear/osh_nuclear_preeq.{c,h}`: single-component exciton
  model per Section 5A (Williams densities, detailed-balance widths, γ_j with
  per-species scale knob, λ₊ CEM-form, n_eq = √(2gE) exit, Kalbach angular
  sampling); emits into the existing `osh_nuclear_event` (whitelist
  {n, p, d, t, ³He, α} — same as FBU); equilibrated residue routed to FBU via
  `osh_nuclear_compound_step()`, unifying ion and neutron paths.
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
- Go/no-go: hard tail (> 20 MeV) appears with roughly the reference slope
  (SH12A + measured DDX); fluence-weighted plateau mean reaches ~10 MeV.
- Decision point: if Stage 2 alone already meets the spectral targets
  (possible — CEM describes 160 MeV (p, xα) data with pre-equilibrium alone),
  demote Stage 3 to a backlog item.

### Stage 4 — Joint calibration + validation

- Tune {P_ret, g (via a), γ_j scales, P_ko} against: #212 SH12A spectra, the
  measured n/p + C/O light-ion DDX (Section 8), and TENDL p+¹²C/¹⁶O emission
  spectra generated offline by a small `tools/` script.
- Acceptance: table in Section 1.  Also verify: total inelastic multiplicity
  and nucleon spectra unregressed; neutron-path compound events unchanged;
  DLET vs depth moves toward SH12A; runtime cost of the nuclear step < a few
  % of total (inelastic events are rare; expected impact ≈ none).

### Later (separate issues)

- π/Δ production ≥ 280 MeV (removes the H-inelastic ABSORB fallback).
- SMM / evaporation for A > 16 residues (#176) — the exciton module's
  equilibrium exit is exactly where they plug in.
- Ion-projectile fast stage; only there, revisit coalescence (C) with the
  verified p_c radii.
- Two-component exciton refinement (proton/neutron exciton bookkeeping à la
  Koning–Duijvestijn) if isotopic yields (t vs ³He) ever matter.

## 8. Calibration and validation data

Measured light-ion (p, d, t, ³He, α) double-differential and
energy-differential cross sections on the OSH-relevant light targets — these
anchor the γ_j scales and P_ko independently of SH12A:

| Dataset | Reaction | Energies | Why it matters |
|---|---|---|---|
| Tippawan et al., PRC 73 034611 (2006) | n + ¹⁶O → lcp, 8 angles | 96 MeV | closest measured proxy to therapy p+O; α spectra to ~70 MeV |
| Tippawan et al., PRC 79 064611 (2009) | n + ¹²C → lcp | 96 MeV | same for carbon |
| Bevilacqua et al. (TSL Medley), arXiv:1303.4637 | n + O/Si/Fe/Bi → lcp | 175 MeV | highest-energy O data; paper *is* a cross-model shoot-out (INCL4.5-ABLA07 / MCNP6 / TALYS / PHITS±surface-coalescence) |
| Benck et al., At. Data Nucl. Data Tables | n + ¹⁶O, n + ¹²C → lcp | 25–75 MeV | low-energy anchor for evaporation/FBU region |
| Bertrand & Peelle, PRC 8 1045 (1973) | p + A (A = 12–209) → H/He spectra | 29–62 MeV | proton-projectile anchor incl. carbon |
| Cowley et al. (as used in CEM03.03 validation) | p + Al/Co/Au → α | 160 MeV | right energy, heavier targets — shape/systematics check |
| ICRU Report 63 (2000) | p + C/N/O/Ca evaluated DDX | to 250 MeV | evaluated reference incl. the 200–250 MeV range TENDL lacks |
| TENDL proton sub-library | p + anything, MF6 spectra | to 200 MeV | machine-readable; free; Kalbach-Mann angular shapes included |

Charge symmetry makes the 96/175 MeV *neutron* data a legitimate stand-in for
proton-induced spectra at OSH's accuracy bar (mirror the Coulomb-barrier
shift at the low-energy end).

## 9. Risks and open questions

- **E\* budget double-counting**: retention (Stage 1) and the per-hole
  constant both feed E\*; calibrate jointly, don't stack defaults.
- **Cluster channels eat nucleons**: both CEM (coalescence bookkeeping) and
  INCL (20–50 MeV nucleon dip) warn that adding cluster production perturbs
  the already-validated nucleon spectra; the Stage-0 harness must track
  nucleon observables, not just alphas.
- **γ_j needs fitting, by design**: CEM03.03's own γ_j is "a rather crude
  estimate" corrected by fitted M_j(A, Z, T₀); do not expect the closed form
  to land on data unaided, and beware the TALYS-style α overshoot.
- **Preeq may still undershoot the > 30 MeV tail** (E\*-limited): that is
  what Stage 3 is for; the Stage-0 experiment sizes this gap up front, and
  the Stage-3 decision point allows dropping it if unneeded.
- **RNG stream changes**: every stage adds draws; fixed-seed regression
  baselines will shift — plan test updates per stage (statistical asserts,
  not bitwise).
- **SH12A production-spectrum reference**: the 0.3 α/primary SH12A yield in
  Section 1 is inferred (fluence 0.66× + dose 0.96×), not directly tallied;
  Stage 0 should pin it down from the decks.
- **Isospin bookkeeping in preeq**: single-component exciton models need the
  charge-composition factor ℜ_j/R_j; a two-component model is not worth the
  size at OSH's accuracy bar (see "Later").

## 10. References

Primary sources verified during this study are marked ✓.

- V. D. Toneev, K. K. Gudima, *Particle emission in light and heavy-ion
  reactions*, Nucl. Phys. A **400** (1983) 173c — cascade + coalescence
  (p_c criterion).
- K. K. Gudima, S. G. Mashnik, V. D. Toneev, *Cascade-exciton model of
  nuclear reactions*, Nucl. Phys. A **401** (1983) 329 — Modified Exciton
  Model.
- ✓ S. G. Mashnik, K. K. Gudima, R. E. Prael, A. J. Sierk, M. I. Baznat,
  N. V. Mokhov, *CEM03.03 and LAQGSM03.03 Event Generators for the MCNP6,
  MCNPX, and MARS15 Transport Codes*, LA-UR-08-2931, arXiv:0805.0751 —
  verified: p_c values (Eq. 19), coalescence-is-minor for nucleon-induced
  reactions (Fig. 17 & §4), MEM formulas (Eqs. 20–31), γ_j and fitted M_j,
  n_eq (Eq. 34), Kalbach angular systematics used below 210 MeV, broad
  initial (p, h, E\*) from cascade, exact per-event conservation.
- ✓ A. Boudard, J. Cugnon, J.-C. David, S. Leray, D. Mancusi, *New
  potentialities of the Liège intranuclear cascade model for reactions
  induced by nucleons and light charged particles*, Phys. Rev. C **87**
  (2013) 014606, arXiv:1210.3498 — verified: surface-coalescence parameters
  (Table I), cascade-dominance of lcp production at 63 MeV, "back to
  spectator" re-absorption (ξ ≈ 7 MeV + ⅔ V_C), cluster channels eating
  20–50 MeV nucleons.
- C. Kalbach, *Systematics of continuum angular distributions*, Phys. Rev. C
  **37** (1988) 2350 — angular systematics (used by CEM and the evaluated
  libraries).
- C. Kalbach, *Preequilibrium reactions with complex particle channels*,
  Phys. Rev. C **71** (2005) 034606 — pickup/stripping/knockout terms for
  cluster emission (grounding for Stage 3; used by TALYS).
- A. Iwamoto, K. Harada, Phys. Rev. C **26** (1982) 1821 —
  exciton-coalescence (condensation) picture for cluster pre-equilibrium
  emission.
- A. J. Koning, M. C. Duijvestijn, Nucl. Phys. A **744** (2004) 15 —
  two-component exciton model; A. J. Koning et al., *TALYS: modeling of
  nuclear reactions*, Eur. Phys. J. A **59** (2023) 131.
- ✓ U. Tippawan et al., Phys. Rev. C **73** (2006) 034611
  (arXiv:nucl-ex/0501014) — 96 MeV n + ¹⁶O light-ion DDX (verified scope);
  and Phys. Rev. C **79** (2009) 064611 (arXiv:0812.0701) — 96 MeV n + ¹²C.
- ✓ R. Bevilacqua et al., *Light-ion production from O, Si, Fe and Bi induced
  by 175 MeV quasi-monoenergetic neutrons*, arXiv:1303.4637 — verified:
  cross-model comparison, pre-equilibrium composite emission called the most
  crucial aspect.
- E. Benck et al., At. Data Nucl. Data Tables — lcp production, n + ¹⁶O
  (25–65 MeV) and n + ¹²C (25–75 MeV).
- F. E. Bertrand, R. W. Peelle, Phys. Rev. C **8** (1973) 1045 — p-induced
  H/He spectra, A = 12–209, 29–62 MeV.
- G. Battistoni et al., *The FLUKA code: an accurate simulation tool for
  particle therapy*, Front. Oncol. **6** (2016) 116 — PEANUT chain
  (GINC + coalescence along the chain + exciton preeq + de-excitation).
- M. Fippel, M. Soukup, Med. Phys. **31** (2004) 2263 — simplified nuclear
  interactions for fast proton MC (VMCpro).
- K. Souris, J. A. Lee, E. Sterpin, Med. Phys. **43** (2016) 1700 — MCsquare:
  inelastic final states sampled from ICRU 63 data, alphas transported.
- C. Grassberger, H. Paganetti, *Elevated LET components in clinical proton
  beams*, Phys. Med. Biol. **56** (2011) 6677 — secondary-particle LET_d
  context.
- ICRU Report 63, *Nuclear Data for Neutron and Proton Radiotherapy and for
  Radiation Protection* (2000).
- A. J. Koning, D. Rochman et al., *TENDL: Complete Nuclear Data Library for
  Innovative Nuclear Science and Technology*, Nucl. Data Sheets **155**
  (2019) 1 — proton sub-library: explicit channels ≤ 30 MeV, MF6/MT5
  continuum spectra ≤ 200 MeV.
- J. P. Bondorf et al., Phys. Rep. **257** (1995) 133 — SMM (future, A > 16).
- Geant4 sources used for sizing/extraction: `pre_equilibrium/exciton_model`
  (transition rates, cluster emission widths, condensation factors),
  `inclxx` (clustering algorithm), `cascade/cascade` (Bertini).
