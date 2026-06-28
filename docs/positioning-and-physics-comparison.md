# OpenShieldHIT — Positioning in the Monte Carlo Transport Landscape

*A physics-capability comparison for the first public release.*

> **Read this first.** OpenShieldHIT (OSH) is a young, lean, from-scratch Monte
> Carlo transport code. This document is deliberately honest about what is and
> is not implemented, so that nobody arrives expecting a drop-in replacement for
> FLUKA or Geant4. It also points out the (real) cases where the implemented
> physics is already *good enough* to be useful today.

---

## TL;DR — where OpenShieldHIT sits

OpenShieldHIT is a **high-energy charged-particle transport engine for the
therapy / shielding energy window** (roughly a few MeV/u up to ~GeV/u), with a
working ion → secondary → neutron chain, true 3D combinatorial + DICOM-CT
geometry, and clinical dose/LET/fluence scoring.

- It is **not** a low-energy ion-implantation / radiation-damage tool — that is
  SRIM/TRIM territory, and OSH does not try to compete there.
- It is **not yet** a general-purpose multi-particle code — there is no
  electron, positron, or photon transport, and the hadronic models are
  deliberately compact rather than full intranuclear-cascade.
- It **is** a transparent, fast, modern, MIT-licensed engine that already
  reproduces the dominant physics of a proton/ion Bragg curve in water and on a
  patient CT, end to end, in a code base you can actually read.

Think of it as occupying the niche of **"a small SHIELD-HIT12A core for ions in
the therapeutic range"** — strong where that domain is well-defined, explicitly
incomplete elsewhere.

---

## Capability snapshot

| Area | Status in OSH today |
|---|---|
| **Primary ions** (p, d, t, ³He, α, ¹²C, heavier) | ✅ Full CSDA wavefront transport |
| **Neutrons (fast)** | ✅ Condensed model, JEFF-4.0 tables (35 nuclides) + Tripathi fallback |
| **Neutrons (thermal, < 1 eV)** | ❌ No separate thermal regime / S(α,β) kernels |
| **Photons** | ❌ Stub only (returns "not implemented") |
| **Electrons / positrons** | ❌ Not modelled (no δ-ray transport) |
| **Electronic stopping power** | ✅ Bethe-Bloch + Lindhard-Scharff low-E join + Hubert effective charge + Sternheimer density correction |
| **Multiple Coulomb scattering** | ✅ Highland (Gaussian-Molière) + random-hinge; Molière-specific tail reserved but not yet distinct |
| **Energy straggling** | ✅ Bohr Gaussian; Vavilov reserved but not yet distinct |
| **Nuclear inelastic (attenuation)** | ✅ Tripathi σ_R parametrisation |
| **Nuclear elastic** | ⚠️ pp elastic only (recoil protons tracked); pn / pα / … not yet |
| **Fragmentation / de-excitation** | ⚠️ Abrasion + sequential-binary Fermi break-up for residuals A ≤ 16; A > 16 counted but not de-excited; no full SMM/INC |
| **Geometry — combinatorial (CSG)** | ✅ SPH, BOX, RPP, RCC, REC, TRC, ELL, WED, ARB, planes + Boolean zones |
| **Geometry — voxel / DICOM CT** | ✅ Schneider & Permatassari HU calibration, IEC 61217 gantry/couch |
| **Scoring meshes** | ✅ Cartesian (X,Y,Z) and cylindrical (R,Z) |
| **Scored quantities** | ✅ Dose, DoseGy, Fluence, Energy, dose-/track-avg LET, dose-/track-avg (z_eff/β)² |
| **Spectral / differential scoring** | ✅ Single & double differential in E, E/u, LET, (z_eff/β)² |
| **Output formats** | ✅ BDO2019, TEXT, DICOM RTDOSE round-trip |
| **Beam delivery** | ✅ Pencil/Gaussian/uniform spots, spot scanning, SAD fan-out, divergence, external spectrum |
| **Range modulator / ridge filter** | ❌ Not yet |
| **Phase-space (MCPL) I/O** | ❌ Not yet |
| **RTSTRUCT / structure scoring** | ❌ Not yet |
| **Parallelism** | ⚠️ Groundwork laid (per-worker accumulators, RNG keyed by history); thread pool in progress |
| **Platforms** | ✅ Linux · macOS · Windows (CI-tested every commit) |

✅ available · ⚠️ partial / approximate · ❌ not yet

---

# Part 1 — OpenShieldHIT vs SRIM / TRIM

The user request expected SRIM/TRIM to be "on a similar level." In terms of
*ambition and footprint* (a focused, lean tool rather than a giant framework)
that is fair. But it is important to be precise: **SRIM/TRIM and OpenShieldHIT
solve largely different problems, in different energy regimes.** They overlap
only at the edges.

### What SRIM/TRIM is built for

SRIM (Stopping and Range of Ions in Matter) and its Monte Carlo companion TRIM
(Transport of Ions in Matter) by Ziegler, Biersack & Littmark are the de-facto
standard for **low-energy ion–solid interaction**: implantation, radiation
damage, and sputtering, typically from ~10 eV/u up to ~MeV/u.

TRIM's strengths — none of which OSH targets:

- **Best-in-class stopping tables** at low energy, including the regime where
  shell corrections, the Barkas term, and charge-state effects dominate.
- **Nuclear (elastic) stopping** against lattice atoms — the binary-collision
  recoil cascade. This is *the* point of TRIM.
- **Damage metrics**: vacancies, replacement collisions, displacements per atom
  (dpa, NRT), and the ionization/phonon energy partition.
- **Sputtering yields** and backscattering for surface physics.
- Detailed range and lateral/longitudinal **straggling distributions** for
  implant profiles.

### What SRIM/TRIM cannot do (and OSH can)

- **No nuclear reactions / fragmentation.** TRIM has no inelastic hadronic
  channel — no proton-induced spallation, no projectile fragmentation, no
  secondary neutrons. OSH models exactly these (Tripathi attenuation, pp
  elastic, abrasion + Fermi break-up, a fast-neutron pool).
- **No true 3D geometry.** TRIM targets are stacked planar layers; there is no
  combinatorial CSG or patient CT. OSH has both.
- **Therapeutic / shielding energies.** TRIM is impractical and physically
  off-design for the hundreds-of-MeV/u therapy window and GeV/u shielding
  problems that are OSH's home turf.
- **Clinical workflow.** Spot scanning, SAD fan-out, dose-to-water, RTDOSE
  output, gantry/couch geometry — out of scope for TRIM, native to OSH.
- **Speed at high energy / thick targets.** TRIM tracks every recoil in the
  cascade and is slow for thick, high-energy problems; OSH uses a
  condensed-history CSDA wavefront tuned for exactly that case.

### Where OSH is currently *weaker* than SRIM/TRIM — be honest

- **Low-energy stopping accuracy.** OSH uses Bethe-Bloch with a Lindhard-Scharff
  low-energy join and Hubert effective charge. This is excellent in the
  therapeutic window but **less accurate than SRIM below ~1–2 MeV/u**, where
  SRIM's semi-empirical fits are the reference.
- **No nuclear stopping / recoil cascades.** OSH treats target-atom elastic
  recoil only through MCS deflection of the projectile; it does not track
  lattice-atom cascades, so it produces **no damage, dpa, or sputtering data**.
- **No implantation profiles below the therapy regime** — OSH is not designed to
  give you an implant depth distribution at keV energies.

### One-line positioning vs SRIM/TRIM

> **Different energy regimes, different jobs.** Use SRIM/TRIM for low-energy
> implantation, material damage, and sputtering; use OpenShieldHIT for
> high-energy ion-beam dose, range, LET, and secondary production in 3D and on
> CT. OSH is not a SRIM replacement, and SRIM cannot do what OSH does at therapy
> energies.

---

# Part 2 — OpenShieldHIT vs full-scale MC codes

Here the framing flips: against Geant4, FLUKA, SHIELD-HIT12A, PHITS, and MCNP6/X,
OpenShieldHIT is the **small, young, focused newcomer.** It does *a slice* of
what those codes do, and does it transparently and fast — but the slice is
genuinely smaller, and we should say so plainly.

### Feature comparison

| Capability | OpenShieldHIT | SHIELD-HIT12A | FLUKA | Geant4 | PHITS | MCNP6/X |
|---|---|---|---|---|---|---|
| **License / openness** | ✅ MIT, fully open, readable | Closed source | Closed (free binary) | ✅ Open (toolkit) | Closed (licensed) | Export-controlled |
| **Ion transport (therapy E)** | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| **Electron/photon EM (EGS-class)** | ❌ | ⚠️ limited | ✅ | ✅ | ✅ | ✅ |
| **Neutron transport (full)** | ⚠️ fast, condensed | ✅ | ✅ | ✅ | ✅ | ✅ (reference-grade) |
| **Thermal neutrons / S(α,β)** | ❌ | ✅ | ✅ | ✅ | ✅ | ✅ |
| **Nuclear inelastic model** | ⚠️ Tripathi σ_R + abrasion/Fermi | ✅ (full) | ✅ PEANUT/INC | ✅ (BIC/INCL++) | ✅ (INCL/JQMD) | ✅ (CEM/LAQGSM) |
| **Projectile fragmentation** | ⚠️ light residues only | ✅ | ✅ | ✅ | ✅ | ✅ |
| **Full INC / QMD** | ❌ | ✅ | ✅ | ✅ | ✅ | ✅ |
| **Combinatorial geometry** | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ (gold standard) |
| **DICOM CT / clinical workflow** | ✅ native | ✅ | via FLAIR/tools | via GATE/TOPAS | ✅ | via wrappers |
| **Magnetic-field tracking** | ❌ | ⚠️ | ✅ | ✅ | ✅ | ⚠️ |
| **Variance reduction** | ❌ | ⚠️ | ✅ extensive | ✅ | ✅ | ✅ extensive |
| **Decades of validation** | ❌ (new) | ✅ | ✅ | ✅ | ✅ | ✅ |
| **Lines of code to understand it** | ~tens of k, traceable | large | very large | millions | large | very large |
| **Build & run footprint** | Tiny, single binary | Moderate | Moderate | Heavy | Moderate | Moderate |

✅ full · ⚠️ partial/approximate · ❌ absent

### Honest read of this table

- **Physics breadth:** the full-scale codes win decisively. They transport
  electrons, positrons, and photons; they have full intranuclear-cascade /
  QMD / evaporation chains; they handle thermal neutrons, magnetic fields, and
  variance reduction. OSH has none of these yet.
- **Validation maturity:** Geant4, FLUKA, PHITS, MCNP, and SHIELD-HIT12A carry
  decades of benchmarking and published validation. OSH is brand new and
  unvalidated by comparison — this is the single most important caveat for a
  user deciding whether to trust a number.
- **Hadronic depth:** where a full code de-excites every residual nucleus
  through evaporation/SMM and tracks the whole fragment spectrum, OSH currently
  de-excites only light residues (A ≤ 16) via a sequential-binary Fermi
  break-up placeholder and counts the rest.

### Where OpenShieldHIT genuinely competes

It is not all caveats. For a focused set of problems OSH is already attractive:

- **Transparency & pedagogy.** The entire transport chain — parsing, geometry
  compilation, stepping, scoring — is traceable through a small number of
  readable C files. No other code on this list can claim that. For teaching what
  a therapy MC actually does, OSH is arguably *better* than the giants precisely
  because it is small.
- **Speed & footprint.** A single small binary, no heavy framework, AVX2 geometry
  kernels, no hot-path allocations. For sweeping a clean proton/ion Bragg-curve
  or simple-geometry study, it starts instantly and runs lean.
- **Modern engineering.** Cross-platform CI (Linux/macOS/Windows on every
  commit), explicit memory layout, MIT license, embeddable C API. Easy to fork,
  extend, or wire into a pipeline.
- **Clinical I/O out of the box.** DICOM CT in, RTDOSE out, IEC 61217
  positioning, SH12A-compatible BDO and dose-to-water — without bolting on a
  separate framework.

### Where the implemented physics is already "good enough"

These are legitimate strengths to advertise, with the boundary stated:

1. **Pristine proton/ion Bragg peaks in water and homogeneous media.** Bethe-Bloch
   + Lindhard + Hubert + Sternheimer reproduces range and the dominant
   energy-loss curve well in the therapeutic window. For range/dose verification
   with nuclear reactions off, results track analytic and reference expectations.
2. **Lateral penumbra from MCS.** Highland + random-hinge gives a sound Gaussian
   spreading of the beam — good enough for penumbra and field-width studies.
3. **Primary-fluence attenuation from nuclear reactions.** Tripathi σ_R captures
   the loss of primaries with depth (the characteristic pre-peak dose build-down
   and tail), which matters for SOBP and carbon studies.
4. **LET and spectral scoring.** Dose- and track-averaged LET, plus differential
   dΦ/dE and dΦ/dLET (including LET in a chosen material such as silicon),
   are implemented and SH12A-layout-compatible — directly useful for
   microdosimetry-adjacent and detector-response studies.
5. **CT-based dose on a patient phantom.** The full DICOM → transport → RTDOSE
   path works end to end with correct HU calibration and gantry/couch geometry.

### Where it is explicitly *not* good enough yet — set expectations

- **No electrons/photons** → cannot do photon teletherapy, brachytherapy,
  bremsstrahlung, or anything needing EM showers.
- **Approximate nuclear secondaries** → fragment spectra and out-of-field /
  neutron dose are not yet quantitative; heavy-residue de-excitation is missing.
- **No thermal neutrons, no variance reduction, no magnetic fields.**
- **Not validated** → no published benchmark suite yet; treat all absolute
  numbers as provisional until a validation campaign is published.
- **Straggling/MCS model selection is on/off in practice** → the `Vavilov` and
  distinct `Molière` switches are reserved in the input but currently resolve to
  the Bohr-Gaussian and Highland models respectively.

---

## Suggested release-announcement blurb

> **OpenShieldHIT 0.x** is a lean, open-source (MIT) Monte Carlo engine for
> **high-energy ion-beam transport** — protons and heavier ions in the
> therapy/shielding energy range. It already delivers end-to-end ion transport
> with Bethe-Bloch stopping, Highland multiple scattering, Bohr straggling,
> Tripathi nuclear attenuation, pp elastic recoils, light-fragment break-up, a
> fast-neutron pool, true 3D combinatorial **and** DICOM-CT geometry, and
> clinical dose/LET/fluence scoring with RTDOSE output — all in a code base
> small enough to read in an afternoon.
>
> **What it is:** a fast, transparent, modern, embeddable ion-transport core and
> a teaching reference for how a therapy MC works under the hood.
>
> **What it is not (yet):** it does **not** transport electrons or photons, does
> **not** include a full intranuclear cascade or thermal-neutron physics, and
> has **not** undergone the decades of validation behind FLUKA, Geant4,
> SHIELD-HIT12A, PHITS, or MCNP. It is also **not** a low-energy
> implantation/damage tool — for that, use SRIM/TRIM.
>
> If you need a small, hackable, open ion-transport engine for the therapeutic
> window — or a transparent platform to learn and extend — OpenShieldHIT is for
> you. If you need turnkey, fully validated, all-particle physics today, reach
> for one of the established codes. We are building toward that, in the open.

---

*Generated for the OpenShieldHIT first-release positioning discussion. Reflects
the state of the `main`/comparison branch at the time of writing; see
[`TODO.md`](../TODO.md) and [`src/physics/README.md`](../src/physics/README.md)
for the live roadmap.*
