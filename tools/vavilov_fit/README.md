# vavilov_fit — offline straggling sampler tooling

Offline pipeline that produces the **reference fixtures** and (later) the fitted
inverse-CDF **coefficients** for the Vavilov + Landau energy-straggling samplers.
Nothing here is compiled into the OpenShieldHIT C build.

## Provenance / licensing

Everything is **clean-room from published theory** — no source code or tabulated
numbers from GEANT3 (GPL) or from SHIELD-HIT12A / Thomsen's `libvav`:

- `vavilov_exact.py` — exact Vavilov density from its published definition
  (Vavilov, Sov. Phys. JETP 5, 749 (1957); B. Schorr, Comput. Phys. Commun. 7,
  215 (1974)), evaluated via the imaginary-axis Bromwich (Fourier) integral with
  `scipy.special.sici`. Self-validated to machine precision against the
  closed-form moments `mean = γ_E − 1 − β² − ln κ`, `var = (1 − β²/2)/κ`.
- `gen_fixtures.py` — writes `tests/fixtures/vavilov/*.csv`: golden
  `λ(u; κ, β²)` inverse-CDF tables the C unit tests validate the runtime sampler
  against. Vavilov regime κ ∈ [0.01, 10]; u ≤ 0.999.

The Landau regime (κ < 0.01) is cross-checked against `scipy.stats.landau` after
mapping its parametrization to the physics λ convention (affine shift/scale).

## Requirements (dev only)

    python3 -m pip install numpy scipy

## Usage

    python3 tools/vavilov_fit/vavilov_exact.py     # moment self-test
    python3 tools/vavilov_fit/gen_fixtures.py       # (re)generate fixtures
