# Neutron cross-section tables

OpenShieldHIT embeds a condensed neutron cross-section table for the Tier-1
nuclides listed below.  The table is generated from JEFF-4.0 PENDF0K at 0 K
with `tools/condense_neutron_xsec.py` and is stored in
`src/physics/neutron/osh_neutron_xsec_data.h`.

## Energy grid

The current table uses a 31-point irregular energy grid from `1e-9 MeV`
to `20 MeV`.

| Region | Grid density | Purpose |
|---|---:|---|
| `1e-9` to `1 MeV` | 3 points per decade | thermal, epithermal, resonance, and low-energy fast region |
| `1 MeV` to `20 MeV` | 2 points per decade | fast-neutron region where the stored channels are smoother |

The lower edge, `1e-9 MeV`, is `1 meV`.  This preserves the low-energy `1/v`
behavior of important absorber channels such as `B-10(n,alpha)` and
`Li-6(n,alpha)`.  The current neutron transport loop uses `1e-3 MeV`
(`1 keV`) as the default cutoff when `NEUTRLCUT <= 0`, so the sub-keV table
entries are present for future thermal-neutron transport work and for runs that
lower the cutoff explicitly.

The upper edge, `20 MeV`, follows the JEFF-4.0 evaluation range used here.  For
nuclides outside the table, `osh_neutron_xsec_lookup()` falls back to the
Tier-2 optical approximation described below.

## Stored channels

For each listed Tier-1 nuclide, the generated table stores:

| Channel | ENDF MT | Meaning |
|---|---:|---|
| `tot` | 1 | total cross section |
| `el` | 2 | elastic cross section |
| `nn` | 4 | `(n,n')` inelastic cross section |
| `n2n` | 16 | `(n,2n)` cross section |
| `ng` | 102 | `(n,gamma)` capture cross section |
| `np` | 103 | `(n,p)` cross section |
| `na` | 107 | `(n,alpha)` cross section |

Missing MT channels in JEFF are stored as zero arrays and marked in the
generated header.  Lookup uses log-log interpolation on the irregular grid,
with linear interpolation at threshold crossings where one endpoint is zero.
Energies outside the stored range are clamped to the nearest grid edge.

Material compositions may describe natural elements with `A=0`.  Neutron
lookup resolves those entries to a representative mass number before computing
number density or selecting a Tier-1 table row, for example H-1, C-12, N-14,
O-16, Cl-35, and Ca-40 for the elements used in the water-equivalent phantom
cases.  This keeps natural material definitions finite and prevents them from
falling through to the Tier-2 path solely because the material card omitted an
isotope number.

## Tier-1 nuclides

| Nuclide | Z | A | Main reason for inclusion |
|---|---:|---:|---|
| H-1 | 1 | 1 | water, tissue, dominant elastic target |
| H-2 | 1 | 2 | light water, 2.2 MeV capture line |
| He-3 | 2 | 3 | He-3 proportional counter |
| Li-6 | 3 | 6 | LiF TLD, CLYC detector |
| Li-7 | 3 | 7 | LiF TLD, Li-glass |
| Be-9 | 4 | 9 | beam line window |
| B-10 | 5 | 10 | BF3 detector, borated polyethylene |
| B-11 | 5 | 11 | boron shielding |
| C-12 | 6 | 12 | tissue, carbon target, graphite |
| N-14 | 7 | 14 | tissue, air |
| O-16 | 8 | 16 | water, tissue, air |
| F-19 | 9 | 19 | fluoride compounds |
| Na-23 | 11 | 23 | soft tissue |
| Mg-24 | 12 | 24 | bone mineral |
| Al-27 | 13 | 27 | structural aluminium |
| Si-28 | 14 | 28 | silicon detector |
| P-31 | 15 | 31 | bone mineral |
| S-32 | 16 | 32 | soft tissue |
| Cl-35 | 17 | 35 | soft tissue |
| Ar-40 | 18 | 40 | air |
| K-39 | 19 | 39 | soft tissue |
| Ca-40 | 20 | 40 | bone mineral |
| Fe-56 | 26 | 56 | iron and steel shielding |
| Cu-63 | 29 | 63 | copper nozzle, coils, brass |
| Zn-64 | 30 | 64 | brass, natural zinc isotope |
| Zn-66 | 30 | 66 | brass, natural zinc isotope |
| Zn-68 | 30 | 68 | brass, natural zinc isotope |
| Cd-113 | 48 | 113 | cadmium cover, thermal absorber |
| Cd-114 | 48 | 114 | cadmium foil, abundant Cd isotope |
| W-182 | 74 | 182 | tungsten collimator |
| W-183 | 74 | 183 | tungsten collimator |
| W-184 | 74 | 184 | tungsten collimator |
| W-186 | 74 | 186 | tungsten collimator |
| Au-197 | 79 | 197 | activation foil |
| Pb-208 | 82 | 208 | lead shielding |

## Tier-2 fallback

Nuclides not listed above use a simple optical fallback: Tripathi total nuclear
reaction cross section for the non-elastic part and a geometric elastic
approximation.  Sub-channels are set to zero, and the reaction layer treats the
non-elastic remainder as a generic compound event.  The first lookup for each
missing `(Z,A)` pair emits a diagnostic message, capped by the warning tracker
in `struct osh_neutron_xsec`.
