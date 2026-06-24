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
cutoff, using `1e-3 MeV` (`1 keV`) when the input leaves `NEUTRLCUT` at `0.0`.
The embedded cross-section tables already extend below that to `1e-9 MeV`
(`1 meV`) for future thermal-neutron transport, but
thermal-neutron physics is not modeled as a separate regime yet.
Cross-section lookup uses the available neutron tables at the requested energy.

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

## Boundary And Material Cases

Per neutron slot, the current loop handles:

- energy at or below the effective neutron cutoff: kill the neutron
- outside geometry: kill the neutron
- blackhole material: kill the neutron
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

- local neutron energy deposits are marked but not yet scored
- charged secondaries from neutron reactions are not yet fed back into the ion
  transport family
- thermal-neutron models are not implemented separately
- angular distributions are simplified; elastic scattering currently uses an
  isotropic center-of-mass approximation
- the transport scheduler still needs full cross-family orchestration around
  the ion and neutron families

## Cross-Section Data

The Tier-1 neutron cross sections are condensed from JEFF-4.0 PENDF0K at 0 K.
They use a 31-point irregular energy grid, giving 30 interpolation intervals,
from `1e-9 MeV` (`1 meV`) to `20 MeV`.  The lower part of the grid keeps
thermal and epithermal absorber behavior available even though the current
transport loop cuts off at `1 keV`.

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
