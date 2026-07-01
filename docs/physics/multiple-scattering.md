# Multiple Coulomb scattering reference

This page documents the multiple Coulomb scattering (MCS) model in
OpenSHIELDHIT: how a charged ion's direction is changed by the many small-angle
Coulomb deflections it undergoes while traversing matter.

The implementation lives in `src/physics/atomic/`:

```text
osh_physics_scat            dispatcher: enum osh_mcs_model + osh_physics_mcs_scatter()
osh_physics_scat_highland   Highland θ₀ (the scattering width / magnitude)
osh_physics_scat_moliere    full Bethe-Molière angular distribution (shape: core + tail)
```

and it is invoked from the random-hinge phase of the ion stepper
(`ion_step_hinge_and_scatter()` in `src/transport/osh_transport_ion_step.c`).

## Models (the `MSCAT` switch)

`MSCAT` in `beam.dat` selects the angular model. Integer values match
SHIELD-HIT12A.

| mode | name | description |
|------|------|-------------|
| 0 | off | no scattering |
| 1 | Gaussian | Highland θ₀, sampled as two Gaussian projected angles — Gaussian *core only*, no tail |
| 2 | Molière | full Bethe-Molière: Gaussian-like core **plus** the Rutherford single-scattering tail |
| 3 | Wentzel | **reserved** for a future Geant4 Wentzel-VI/Urban screened model; rejected at validation today |

"Gaussian" and "Highland" are not two different models — Highland *is* the formula
for the width of the Gaussian. The meaningful distinction between mode 1 and mode 2
is whether the heavy single-scattering tail is present.

## Strategy: Highland sets the width, Molière sets the shape

The central design decision is to separate the **magnitude** of the scattering
from its **angular distribution shape**.

### Width — Highland θ₀

The Highland formula

```
θ₀ = (13.6 MeV / βcp) · z · √(d/X₀) · [1 + 0.038 · ln(z²·s/X₀)]
```

gives the RMS projected scattering angle and is the validated scattering *power*:
it agrees with experiment (e.g. Gottschalk's proton lateral-spread data) across
materials. **Both** mode 1 and mode 2 use θ₀ for their width.

Two length scales appear deliberately:

- the **substep thickness** `d` in the leading `√(d/X₀)`, so that independent
  per-substep Gaussian variances add correctly across the hundreds of
  energy-loss substeps in a track;
- the **macroscopic path scale** `s` (the residual CSDA range at the substep
  entry) in the logarithmic correction, so that the summed per-substep variance
  reproduces the *full-path* Highland value rather than under-counting it
  (Lynch & Dahl). Passing `s ≤ 0` falls back to `d` (legacy single-step
  behaviour).

This path-scale split is verified by a unit test: splitting a step into N
substeps and summing the variances reproduces the single full-step θ₀².

### Shape — Bethe-Molière

Mode 2 adds the genuine angular distribution. In reduced angle ϑ (with physical
space angle θ = scale·ϑ), Bethe's expansion is

```
f(ϑ) = f⁽⁰⁾(ϑ) + f⁽¹⁾(ϑ)/B + f⁽²⁾(ϑ)/B²
```

where `f⁽⁰⁾(ϑ) = 2·e^(−ϑ²)` is the Gaussian core and the f⁽¹⁾, f⁽²⁾ corrections
add the Rutherford tail. The B parameter solves Molière's equation
`B − ln B = ln Ω`, with `Ω = χ_c²/(1.167·χ_a²)` the effective number of
scatters; χ_c² (the characteristic angle) and χ_a² (the screening angle) follow
the standard Bethe / Lynch-Dahl forms.

### Anchoring the Molière width to Highland

A subtlety: the Molière reduced distribution carries `⟨ϑ²⟩ ≈ 1.2–1.4` (the
second moment is inflated by the Rutherford tail), not 1. If one set the
Molière *core* width equal to Highland θ₀ and then sampled the full distribution,
the tail would push the overall RMS ~15 % above Highland — double-counting the
tail's contribution to the width.

Instead, the sampled angle is **rescaled so its RMS equals Highland θ₀**:

```
θ = θ₀ · √(2 / ⟨ϑ²⟩(B)) · ϑ
```

(the factor 2 converts the projected θ₀ to a space-angle variance). Highland
fixes the magnitude; Molière supplies only the shape (sharp core + heavy tail).
`⟨ϑ²⟩(B)` is precomputed alongside the inverse-CDF table (below).

The payoff is internal consistency: mode 1 and mode 2 have the **same width**,
differing only in that mode 2 has a sharper core and a Rutherford tail. Both
match experiment. (For comparison, SH12A's mode-1 Gaussian is noticeably
narrower than its own mode-2 Molière; OpenSHIELDHIT's two modes agree.)

## Performance: precomputed inverse-CDF

The reduced-angle density depends on B only through `f⁽⁰⁾ + f⁽¹⁾/B + f⁽²⁾/B²`, so
the inverse CDF ϑ(u) is a smooth function of (B, u). At startup
(`osh_physics_moliere_init`):

1. the f⁽⁰⁾/f⁽¹⁾/f⁽²⁾ grids are evaluated from their Bethe integral definitions
   (numerical Hankel/Bessel quadrature — our own code, no external tables);
2. for a linear grid of B nodes, the clamped cumulative of `ϑ·f(ϑ;B)` is built
   and inverted onto a uniform u grid, giving a 2-D support table
   `invcdf[B][u]`; `⟨ϑ²⟩(B)` is accumulated at the same time.

On the hot path, `osh_physics_moliere_sample_reduced` is then a single **O(1)
bilinear lookup** — no per-scatter CDF reconstruction. The B grid covers
`[1.6, 40]`; computed B values stay in `~9–16` for materials from water to
tungsten and 1–200 MeV, so clamping is never reached in practice.

## Material independence

Nothing is hardcoded to a specific material. The two per-medium constants

- `moliere_chic2 = 0.157 · Σᵢ wᵢ Zᵢ(Zᵢ+1)/Aᵢ`  (the χ_c² coefficient), and
- `moliere_screen_z` (an effective screening Z for χ_a),

are derived from the element composition at material-compile time
(`compute_moliere_constants()` in `src/material/runtime/osh_material_compile.c`)
and stored in `struct osh_material_runtime`, next to the radiation length X₀.
The f⁽ⁿ⁾ tables and the inverse-CDF are universal (parameterised only by B).

Because the width is anchored to Highland θ₀ (which uses the per-material X₀),
material-level approximations — notably the single effective screening Z for
compounds spanning very different Z — only shift the tail/core split via B; they
do not change the overall scattering width.

## Random hinge

The chosen deflection is applied by the **random-hinge method** (Fippel &
Soukup 2004): a single angular kick placed at a uniformly random point along the
step. This both distributes the deflection along the step and reproduces the
correct lateral displacement, so the angular models above only need to supply the
deflection — not the offset.

## Provenance (licensing)

This is a **clean-room** implementation built solely from the published
references below. No source code or tabulated numbers were taken from Geant3's
`gmolie.F`, SHIELD-HIT12A, or third-party reimplementations (Geant 3.21 is
GPLv3; OpenSHIELDHIT is MIT). The reduced-angle functions are computed from
their integral definitions, and the sampling/interpolation routines are our own.

## Validation

For 200 MeV protons in water (NUCRE off, STRAGG off — MCS in isolation), lateral
RMS radius near the Bragg peak:

| | lateral RMS |
|---|---|
| OpenSHIELDHIT Gaussian (mode 1) | 0.531 cm |
| OpenSHIELDHIT Molière (mode 2) | 0.528 cm |
| SHIELD-HIT12A Molière (mode 2) | 0.533 cm |
| Gottschalk rule of thumb (≈2 %·range) | ~0.52 cm |

The benchmark decks are `tests/reference/idd_water_200mev_scat{0,1,2}` (and the
SH12A-syntax mirrors under `tests/reference/shieldhit/`); the comparison plots
are produced by `tools/plot_mcs_scat.py`. Unit tests in
`tests/unit/test_osh_physics_moliere.c` cover the Highland path-scale additivity,
the B equation, the analytic f⁽⁰⁾ and the f⁽ⁿ⁾ normalisation moments, the
inverse-CDF monotonicity, and the sampled-distribution moments.

## References

- V. L. Highland, *NIM* **129**, 497 (1975) — Gaussian width
- G. Molière, *Z. Naturforsch.* **2a**, 133 (1947); **3a**, 78 (1948) — theory
- H. A. Bethe, *Phys. Rev.* **89**, 1256 (1953) — reduced-angle functions, B equation
- W. T. Scott, *Rev. Mod. Phys.* **35**, 231 (1963) — review
- G. R. Lynch, O. I. Dahl, *NIM B* **58**, 6 (1991) — χ_c, screening, path-scale additivity
- B. Gottschalk, *Radiation Dosimetry: Proton Therapy* (2010) — proton lateral spread
- M. Fippel, M. Soukup, *Med. Phys.* **31**, 2263 (2004) — random-hinge method
