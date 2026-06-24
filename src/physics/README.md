# Physics modules

This directory contains the physics models used during ion transport.
Each subdirectory is a self-contained static library with no dependency on
the transport loop itself — the models are pure functions that take numbers
in and return numbers out.

```
physics/
    atomic/    electromagnetic (EM) interactions with target electrons
    nuclear/   hadronic interactions with target nuclei
    neutron/   neutron cross sections and reaction final states
```

---

## atomic/ — Electromagnetic physics

Three processes govern how a charged ion loses energy and changes direction
as it passes through matter.  All three act continuously along every step.

### Stopping power — `osh_physics_bethe`

The **Bethe-Bloch formula** gives the mean energy loss per unit areal density
(mass stopping power, dE/dx [MeV cm²/g]).  It is the dominant process for
therapeutic ion beams and determines the Bragg peak position.

At low energies, where Bethe gives unphysical results, the **Lindhard-Scharff**
extension takes over.  The two are joined tangentially (the "sewing point") at
compile time per (projectile, material) pair.

The **Hubert effective-charge parametrisation** accounts for partial electron
capture by heavy ions at lower energies, reducing their effective charge and
therefore their stopping power.

Key references:
- Bethe, Ann. Phys. 5, 325 (1930); Bloch, ibid. 16, 285 (1933)
- Sternheimer & Peierls, Phys. Rev. B 3, 3681 (1971) — density correction
- Hubert, Bimbot & Gauvin, NIMB 36, 357 (1989) — effective charge
- Lindhard, Scharff & Schiott, Mat. Fys. Medd. 33:14 (1963) — low-energy extension

### Multiple Coulomb scattering — `osh_physics_moliere`

Every Coulomb interaction with a target nucleus deflects the ion slightly.
The cumulative effect over many small deflections is described by the
**Highland formula** (a Gaussian approximation to Molière theory), giving an
RMS projected scattering angle per step.  The transport loop uses the
**random-hinge method** (Fippel & Soukup 2004) to apply this as a single
angular deflection at a randomly sampled point along the step.

Key references:
- Highland, NIMB 129, 497 (1975) — Gaussian approximation
- Molière, Z. Naturforsch. 3a, 78 (1948) — full theory
- Fippel & Soukup, Med. Phys. 31, 2263 (2004) — random-hinge method

### Energy straggling — `osh_physics_straggling`

The Bethe formula gives only the *mean* energy loss.  In reality, statistical
fluctuations in the number of collisions broaden the energy distribution —
this is **straggling**.  The **Bohr variance** gives the width of this
distribution per step; the transport loop samples it as a Gaussian correction
to the exit energy.

Key references:
- Bohr, Mat. Fys. Medd. 18:8 (1948)

---

## nuclear/ — Nuclear (hadronic) physics

Nuclear interactions occur when the projectile comes close enough to a target
nucleus to interact via the strong force rather than the electromagnetic force.
Unlike the continuous EM processes above, nuclear interactions are discrete
stochastic events.

### Total nuclear reaction cross section — `osh_nuclear_tripathi`

The **Tripathi formula** (NASA/TP-1999-209726) gives the total nuclear reaction
cross section σ_R for any projectile–target pair as a function of energy.
At each transport step, the survival probability `exp(−ds/λ)` is evaluated
(where λ is the mean free path derived from σ_R) and compared to a uniform
random draw.  If the primary is killed, its remaining kinetic energy escapes
with the (untracked) nuclear fragments and is not scored locally; only the
ionisation energy deposited along the step before the reaction is recorded.

The formula is purely parametric (no nuclear data tables) and includes a
natural Coulomb threshold: the cross section goes to zero below the Coulomb
barrier, which for protons on oxygen is around 3–4 MeV.

The energy argument is T/A [MeV/nucleon] — kinetic energy divided by the
integer mass number — not MeV/u (which differs by ~0.7 % for protons due
to the proton mass exceeding 1 u).

Key references:
- Tripathi, Cucinotta & Wilson, NASA/TP-1999-209726 (1999)

### Abrasion + Fermi break-up — `osh_nuclear_abrasion`, `osh_nuclear_fermi_breakup`

When the inelastic channel fires, the **abrasion** stage (wounded-nucleon
picture, Bowman–Swiatecki–Tsang 1973) emits fast knock-out nucleons and leaves
an excited prefragment with E* from a per-hole estimate (Gaimard & Schmidt
1991, clamped to the absorbed cascade energy) and a momentum from the event
balance.  Light prefragments (A ≤ 16) are then de-excited by the **Fermi
break-up** stage: sequential binary splits weighted by two-body phase space
(g₁g₂ μ^(3/2) √E_kin), resolved over a startup-compiled channel table built
from the isotope database.  Only n, p, d, t, ³He and α are emitted as
transportable products; other residues are counted in the fragment pool.
The sequential-binary scheme is a semiphysical placeholder for the full
simultaneous n-body Fermi model.

Key references:
- Fermi, Prog. Theor. Phys. 5 (1950) 570
- Gaimard & Schmidt, Nucl. Phys. A 531 (1991) 709

### Coming next

- **SMM** (Statistical Multifragmentation Model, Bondorf et al. 1995):
  de-excitation of heavy residues (A > 16) that the Fermi break-up stage
  leaves unprocessed, and a full simultaneous n-body break-up.
- **BNAB-26**: 26-group neutron cross sections for simple neutron transport.

---

## neutron/ — Neutron physics

The neutron module is the physics layer for neutron transport.  The current
implementation starts with the minimal above-thermal transport path, but the
module boundary is intended to also cover thermal-neutron treatment once those
models are added.  It deliberately sits between the transport loop and the
nuclear back-ends:

- `osh_neutron_xsec` owns neutron cross-section lookup.  It interpolates the
  condensed JEFF-derived tables for the supported Tier-1 nuclides and falls
  back to a Tripathi/geometric approximation for missing nuclides.
- `osh_neutron_reaction` owns target selection, channel sampling, and final
  state construction for one neutron interaction.  The transport loop samples
  the flight distance and manages pools; this layer decides whether an event is
  elastic, capture, charge exchange, compound decay, or a local-deposit sink.

Compound neutron reactions are routed through `osh_nuclear_compound_step()`,
which is the adapter from neutron-induced compound nuclei to the existing
Fermi break-up back-end or to the heavy-residue sink.  This keeps
`osh_transport_neutron.c` free of Fermi break-up details.

The current reaction model is intentionally minimal:

- elastic scattering uses an isotropic CM approximation; H-1 recoils are
  returned as proton secondaries, while heavier recoils deposit their recoil
  energy locally.
- `(n,gamma)` capture kills the neutron and deposits the neutron kinetic energy
  locally.
- `(n,p)` and `(n,alpha)` use two-body relativistic kinematics when the mass
  lookup and threshold allow the channel; otherwise they fall back to local
  deposit.
- `(n,n')`, `(n,2n)`, and generic non-elastic remainder build a compound
  nucleus `(Z,A+1,E*)` and delegate de-excitation to the compound adapter.

Thermal-neutron physics is not yet modeled separately; current lookups simply
use the available cross-section tables at the requested energy.  A later
thermal layer can live under the same `neutron/` ownership boundary.

`osh_transport_neutron_run()` currently drains the neutron pool with a
wavefront loop that mirrors the ion driver at a high level: batch zone lookup,
batch boundary-distance lookup, per-slot cutoff/escape/material handling,
exponential free-path sampling from the macroscopic total cross section, and
reaction sampling through `osh_neutron_reaction_sample()`.  Elastic events keep
the neutron alive with updated direction and energy; compound events may push
neutron secondaries back into the neutron pool for the next wavefront pass; all
other current channels kill the neutron.  Local energy deposits and charged
ion feedback are marked at the event boundary but are not fully scored or
fed back into ion transport yet.

---

## How the modules connect to transport

```
osh_transport_ion_step()
    Phase 2b  ion_step_length()              — CSDA range (Bethe tables)
    Phase 3   ion_step_hinge_and_scatter()   — Highland/Molière angle
    Phase 4   ion_step_energy_and_straggling() — exit energy + Bohr straggling
    Phase 5   ion_step_nuclear()             — Tripathi survival probability
    Phase 6   ion_step_commit()              — score + update pool
```

The EM processes (phases 2–4) act on every material step.  The nuclear
process (phase 5) is only active when `beam->nuclear` is set in `beam.dat`.

For neutron transport, `osh_transport_neutron.c` calls the neutron module only
after it has sampled an interaction point in the current material.  It then
applies the returned `osh_neutron_reaction_event`: update the neutron for
elastic scatter, push neutron secondaries back into the neutron pool, and mark
charged secondaries/local deposits for the ion-feedback and point-scoring paths
that will be wired later.
