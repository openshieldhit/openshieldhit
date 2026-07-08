# Neutron transport reference

This page documents the current neutron-transport model in OpenShieldHIT.  It
is a work-in-progress reference for developers and reviewers: the implementation
now has a minimal neutron pool drain and reaction layer, while several physics
and scoring paths remain explicit TODOs.

## Scope

The current neutron transport path handles neutrons produced by nuclear
interactions in ion transport.  Those neutrons are banked in `struct
osh_neutron_pool` and later drained by `osh_transport_neutron_run()`.

The model currently transports neutrons down to the configured `NEUTRLCUT`
cutoff, using the cross-section table floor `1e-9 MeV` (`1 meV`) when the input
leaves `NEUTRLCUT` at `0.0`.  Elastic collisions that would down-scatter below
room-temperature thermal energy are clamped to `2.53e-8 MeV` (`25.3 meV`), so
default runs keep thermal neutrons alive until capture or escape.  Cross-section
lookup uses the available neutron tables at the requested energy.

## Data Flow

```text
ion nuclear event
  -> neutron secondary pushed to osh_neutron_pool
  -> osh_transport_neutron_run()
      -> geometry zone and boundary-distance batches
      -> macroscopic total cross section
      -> sampled free path
      -> osh_neutron_reaction_sample()
      -> apply event to neutron pool / future ion feedback / future scoring
```

The neutron pool is a structure-of-arrays store with per-slot position,
direction, kinetic energy, statistical weight, ancestor primary index,
generation number, and RNG stream.  Dead entries are marked by setting
`e <= 0` and removed by `osh_neutron_pool_compact()` at the end of each
wavefront pass.

## Flight Sampling

For a neutron with kinetic energy `E` in material `m`, the transport loop
computes the macroscopic total cross section

```text
Sigma_tot(E, m) = sum_i n_i sigma_tot_i(E)
```

where `n_i` is the element number density in `cm^-3` and `sigma_tot_i` is the
microscopic total cross section.  The implementation stores microscopic
cross sections in millibarns and converts the final sum to `cm^-1`.

The sampled interaction distance is

```text
l = -log(U) / Sigma_tot
```

with `U` drawn from the neutron slot's RNG stream.  If `l` is beyond the current
geometry boundary, the neutron is advanced to the boundary plus a small nudge
and survives to the next wavefront pass.  Otherwise it is advanced to the
interaction point and handed to the neutron reaction layer.

Finite neutron transport steps are passed to the existing step-scoring path with
zero continuous energy deposition, enabling neutron fluence maps.  A neutron
that is killed by the configured energy cutoff inside positive-density material
point-deposits its remaining kinetic energy locally, so the cutoff residual is
booked as `Energy` / `Dose` / `DoseGy`.  Other point-like local deposits from
neutron reactions remain separate from this and require additional wiring.

## Boundary And Material Cases

Per neutron slot, the current loop handles:

- outside geometry: kill the neutron
- blackhole material: kill the neutron
- energy at or below the effective neutron cutoff in positive-density material:
  point-deposit the remaining kinetic energy, then kill the neutron
- vacuum or zero-density material: advance to the next boundary
- zero macroscopic cross section: advance to the next boundary
- finite macroscopic cross section: sample a free path and possibly interact

## Reaction Layer

`physics/neutron/osh_neutron_reaction.{h,c}` owns target selection, channel
sampling, and final-state construction for one interaction.  Transport owns
where the interaction happens and how returned products are pushed into pools.

The current channels are:

| Channel | Current behavior |
|---|---|
| elastic | update neutron direction and energy; H-1 recoil proton is returned as a secondary; heavier recoil energy is marked as local deposit |
| `(n,gamma)` | kill neutron and mark local deposit |
| `(n,p)` / `(n,alpha)` | create charged secondary with two-body relativistic kinematics when mass lookup and threshold allow it; otherwise local deposit |
| `(n,n')`, `(n,2n)`, remainder | build compound nucleus `(Z,A+1,E*)` and route through `osh_nuclear_compound_step()` |

Compound events may produce neutron secondaries through Fermi break-up; those
are pushed back into the neutron pool and processed in a later wavefront pass.

## Current Limitations

These are intentional boundaries of the current implementation:

- local neutron-reaction energy deposits are marked but not yet scored
- charged secondaries from neutron reactions ((n,p), (n,α)) are not yet fed
  back into the ion transport family; their energy is deposited locally instead
- thermal-neutron models are not implemented separately (separate issue #178)
- elastic scattering uses an isotropic center-of-mass approximation

## Cross-Section Data

The Tier-1 neutron cross sections are condensed from JEFF-4.0 PENDF0K at 0 K.
They use a 31-point irregular energy grid, giving 30 interpolation intervals,
from `1e-9 MeV` (`1 meV`) to `20 MeV`.  The lower part of the grid keeps
thermal and epithermal absorber behavior available to the current default
transport loop.

See [Neutron cross-section tables](neutron-cross-sections.md) for the stored
channels, the full nuclide list, and the Tier-2 fallback rule.

## Code Map

| File | Responsibility |
|---|---|
| `src/transport/osh_neutron_pool.{h,c}` | neutron pool allocation, reset, compact |
| `src/transport/osh_transport_neutron.c` | neutron wavefront transport loop |
| `src/physics/neutron/osh_neutron_xsec.{h,c}` | microscopic cross-section lookup |
| `src/physics/neutron/osh_neutron_reaction.{h,c}` | target/channel sampling and final states |
| `src/physics/nuclear/osh_nuclear_compound.{h,c}` | compound adapter to Fermi break-up or heavy-residue sink |
