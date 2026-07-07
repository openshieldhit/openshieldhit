# Nuclear inelastic reactions reference

This page documents the proton-induced nuclear reaction chain in
OpenShieldHIT: what happens when the inelastic hazard fires on a transport
step, which models produce the secondary particles, and where the calibration
knobs live.

The implementation lives in `src/physics/nuclear/`:

```text
osh_nuclear_handler        per-step hazard evaluation, struck-element selection, chain dispatch
osh_nuclear_tripathi       total reaction cross section σ_R (Tripathi 1999)
osh_nuclear_pp             p+p (hydrogen) elastic channel
osh_nuclear_elastic        p+A elastic channel (recoil nucleus)
osh_nuclear_abrasion       fast stage: mini intranuclear cascade + knockout retention
osh_nuclear_preeq          pre-equilibrium: single-component exciton model
osh_nuclear_fermi_breakup  equilibrium stage: microcanonical Fermi break-up (A ≤ 16)
osh_nuclear_compound       compound-nucleus router (neutron path; future SMM seam)
```

and is invoked from the ion stepper (`src/transport/osh_transport_ion_step.c`),
which owns the event struct and dispatches the produced secondaries to the
neutron/ion pools.

## The `NUCRE` switch

`NUCRE` in `beam.dat` selects the nuclear channels. Integer values match
SHIELD-HIT12A where they exist.

| mode | channels | note |
|------|----------|------|
| 0 | off | no nuclear interactions |
| 1 | inelastic + elastic | the physical "on" setting |
| 2 | elastic only | OpenShieldHIT-only decomposition of mode 1 |
| 3 | inelastic only | OpenShieldHIT-only decomposition of mode 1 |

## The inelastic chain

An inelastic event runs a staged pipeline (the architecture of issue
[#221](https://github.com/openshieldhit/openshieldhit/issues/221); stage
boundaries are fixed interfaces so that a future real intranuclear cascade
and a coalescence stage for ion projectiles can slot in without rework):

```text
σ_R hazard (Tripathi)  →  struck-element selection (per-nuclide rates)
  →  1. FAST STAGE      abrasion: knockout nucleons + excited prefragment (A, Z, E*, p, h)
  →  2. (coalescence)   reserved slot — no-op for proton projectiles
  →  3. PRE-EQUILIBRIUM exciton model: fast {n, p, d, t, ³He, α} emission
  →  4. EQUILIBRIUM     Fermi break-up (A ≤ 16); heavier residues: sink (future SMM/evaporation)
  →  fragment pool → ion/neutron transport or point deposit
```

Only protons run this pipeline today; ion projectiles resolve as `ABSORB`
(energy deposited locally) until the ion-projectile fast stage lands.
Hydrogen targets are excluded from the inelastic hazard for protons below
280 MeV (no open channel below the pion threshold); the p+p elastic channel
covers them.

### Fast stage — abrasion with knockout retention

`osh_nuclear_abrasion.c` implements a Bowman-Swiatecki-Tsang abrasion in an
intranuclear-cascade picture: the proton undergoes ν quasi-elastic collisions
(⟨ν⟩ = σ_pN·A/σ_pA, zero-truncated Poisson), each knocking out one nucleon
(isotropic CM, full relativistic equal-mass kinematics) and deflecting the
continuing cascade proton. Each hole charges
`OSH_ABRASION_EXCITATION_PER_HOLE_MEV` (13.3 MeV, Gaimard–Schmidt) to the
cascade proton, booked as prefragment excitation, so the kinetic budget
`T_in = Σ KE_escaped + E*` is exact.

Low-energy knockouts are **re-absorbed** ("back to spectator", INCL4.6,
Boudard et al. PRC 87 (2013) 014606 Sec. II D1): a smooth Fermi-function
retention around `OSH_ABRASION_RETENTION_THRESHOLD_MEV` (ξ = 7 MeV; protons
add ⅔ of the Coulomb barrier). Retained kinetic energy funds E\*; the
retained nucleon becomes a **particle exciton**. Every collision leaves a
**hole exciton**; the exciton configuration (p, h) is exported in
`struct osh_nuclear_fragment` as the input of the pre-equilibrium stage.
The threshold default (ξ = 15 MeV) was calibrated at transport level in the
[#225](https://github.com/openshieldhit/openshieldhit/issues/225) study:
retaining the 10–15 MeV knockouts broadens the prefragment E\* tail and
hardens the break-up alpha sector into its acceptance windows, at a measured
cost confined to the already-low sub-20 MeV proton bands (the validated
20–50 MeV band is untouched) — the quantitative form of INCL4.6's caution
about their earlier ξ = 18 MeV variant.

| knob (compile-time) | default | meaning |
|---|---|---|
| `OSH_ABRASION_SIGMA_PN_MB` | 30.0 | p+nucleon σ for ⟨ν⟩ |
| `OSH_ABRASION_EXCITATION_PER_HOLE_MEV` | 13.3 | hole excitation (do not retune — #260) |
| `OSH_ABRASION_RETENTION_THRESHOLD_MEV` | 15.0 | back-to-spectator ξ (0 disables); calibrated in #225 (ξ ∈ {7,12,15,18} transport scan) |
| `OSH_ABRASION_RETENTION_WIDTH_MEV` | 2.0 | smooth turn-on width |

### Pre-equilibrium — single-component exciton model

`osh_nuclear_preeq.c` (issue
[#225](https://github.com/openshieldhit/openshieldhit/issues/225)) de-excites
the prefragment before equilibrium break-up. Starting from (p, h, E\*), each
step either emits an ejectile or thermalizes one Δn = +2 collision, until the
equilibrium exciton number n_eq = √(2 g E\*) is reached:

- **State density**: Williams, ω(p, h, E) = g (gE)^(p+h−1) / (p! h! (p+h−1)!),
  with g = `OSH_PREEQ_LEVEL_DENSITY_PER_A`·A (0.075/MeV per nucleon).
- **Emission widths** (detailed balance, Betak-Dobes form) for
  {n, p, d, t, ³He, α}: exact mass-table separation energies,
  Dostrovsky-style inverse cross sections, a composition factor from the
  residue N/Z, and for clusters the condensation probability γ_j
  (γ_d = 16/A, γ_t = γ_³He = 243/A², γ_α = 4096/A³ — the Geant4 closed
  forms), each times a per-species calibration scale.
- **Transition rate** λ₊ (Δn = +2 only): CEM form — in-medium NN cross
  section at the mean exciton relative energy, Pauli-blocking factor,
  interaction volume (Gudima–Mashnik–Toneev).
- **Contracts**: a fragment with (p, h) = (0, 0) is thermalized and passes
  through untouched (compound nuclei, break-up residues); per-emission
  bookkeeping is exact (E\* before = E\* after + B_j + T); only the
  whitelist {n, p, d, t, ³He, α} is emitted. Emission angles are isotropic
  (Kalbach continuum systematics are a planned refinement).

| knob (compile-time) | default | meaning |
|---|---|---|
| `OSH_PREEQ_LEVEL_DENSITY_PER_A` | 0.075 | single-particle level density / A [1/MeV] |
| `OSH_PREEQ_FERMI_ENERGY_MEV` | 35.0 | Fermi energy in λ₊ |
| `OSH_PREEQ_TRANSITIONS_R0_FM` | 0.6 | interaction-volume radius in λ₊ |
| `OSH_PREEQ_EMISSION_R0_FM` | 1.25 | inverse-σ radius parameter |
| `OSH_PREEQ_GAMMA_SCALE_{D,T,HE3,ALPHA}` | 1.0 | cluster condensation calibration scales |

**Calibration status.** The γ scales ship neutral (1.0) and are *expected*
to be fitted (CEM03.03 multiplies its equivalent factors by empirical fits).
The #225 measurement showed that on abrasion-fed exciton configurations the
deuteron channel responds strongly to γ_d while t/³He/α are exciton-starved
(they need p ≥ 3–4 particle excitons), and that cluster emission is
E\*-redistributive against the break-up stage — the hard α component needs a
direct knockout channel (Stage 3 of #221) or a richer cascade. The joint fit
is deferred until the isotope range-table defect
([#267](https://github.com/openshieldhit/openshieldhit/issues/267)) is fixed,
because it distorts the transported d/t/³He fluences the fit targets.

### Equilibrium — Fermi break-up

`osh_nuclear_fermi_breakup.c` de-excites light residues (A ≤ 16) by the
microcanonical statistical Fermi break-up: all N = 2..6 ground-state
partitions over the isotope database, weighted by the canonical phase-space
probability with a free-volume coefficient (single knob
`OSH_FERMI_BREAKUP_R0_FM` = 0.50 fm, calibrated against G4FermiBreakUp
multiplicities). Particle-unstable intermediates (Be-8, He-5, Li-5) decay on
a work stack; only the whitelist {n, p, d, t, ³He, α} is ever transported;
bound non-whitelist nuclides go to the fragment pool. Ground states only —
systematic excited-state level tables are issue
[#196](https://github.com/openshieldhit/openshieldhit/issues/196).

### Compound routing and heavier residues

`osh_nuclear_compound.c` routes an equilibrated compound nucleus: A ≤ 16 to
Fermi break-up, A > 16 to a warn-once sink that deposits E\* locally — the
seam where evaporation
([#175](https://github.com/openshieldhit/openshieldhit/issues/175)) and SMM
([#176](https://github.com/openshieldhit/openshieldhit/issues/176)) will plug
in. The neutron-capture path (`src/physics/neutron/osh_neutron_reaction.c`)
enters here with a thermalized (p, h) = (0, 0) configuration.

## Event contract and transport handoff

The chain writes one `struct osh_nuclear_event`: up to 32 secondaries
(direction, kinetic energy, species pointer) and up to 4 residual fragments
(A, Z, E\*, lab momentum, exciton configuration). Transport injects
secondaries into the neutron/ion pools; fragments above the per-nucleon
threshold enter the ion pool as recoil species, below it they are
point-deposited. Kinetic energy is conserved exactly through abrasion, and
total energy (with masses) through pre-equilibrium and break-up; unit tests
enforce both (`tests/unit/test_osh_nuclear_*.c`).

## Validation tooling

- `tests/reference/idd_water_200mev_nucre{0..3}/` — NUCRE isolation decks
  (200 MeV p → water, MSCAT/STRAGG off) with species-resolved dose/fluence,
  plateau spectra, and LET; SH12A mirror fixtures under
  `tests/reference/shieldhit/`.
- `tools/plot_nucre.py` — OpenShieldHIT vs SH12A overlay report.
- `examples/06_nucre_validation/nucre_scan` — transport-free production
  spectra of the abrasion → preeq → break-up chain (yields per event,
  per-species spectra, prefragment E\* and exciton statistics).

## Provenance and references

The models are implemented from the published literature; Geant4 sources
were consulted as a formula cross-reference only (no code reuse — license
differs from OpenShieldHIT's MIT). SHIELD-HIT12A shares the input format and
application domain, not code; its reference spectra are used as validation
fixtures only.

- Tripathi, Cucinotta, Wilson, NASA/TP-1999-209726 — reaction cross section.
- Bowman, Swiatecki, Tsang (1973); Gaimard & Schmidt, Nucl. Phys. A 531
  (1991) 709 — abrasion, hole excitation.
- Boudard, Cugnon, David, Leray, Mancusi, Phys. Rev. C 87 (2013) 014606 —
  INCL4.6; knockout retention ("back to spectator").
- Gudima, Mashnik, Toneev, Nucl. Phys. A 401 (1983) 329; Mashnik et al.,
  arXiv:0805.0751 (CEM03.03) — exciton model, transition rate, n_eq exit.
- Betak & Dobes, Z. Phys. A 279 (1976) 319 — cluster emission widths.
- Dostrovsky, Fraenkel, Friedlander, Phys. Rev. 116 (1959) 683 — inverse
  cross sections.
- Fermi, Prog. Theor. Phys. 5 (1950) 570 — break-up; validated against
  G4FermiBreakUp (Pshenichnov test data, used with permission).
