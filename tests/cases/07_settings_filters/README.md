# 07_settings_filters

150 MeV proton pencil beam in a PMMA cylinder (R=10 cm, L=20 cm).
100 primaries, Molière scattering, Vavilov straggling, nuclear reactions on.

Scoring: single mesh geometry (21×21×40 cm, 800 z-bins).

## Filters

| Name    | Description                    |
|---------|--------------------------------|
| Primary | Generation-0 protons only      |
| Protons | All protons (primary + secondary) |

## Settings

Material overrides for stopping-power-based corrections:

| Name    | Description              |
|---------|--------------------------|
| inWater | Water stopping power     |
| inSi    | Silicon stopping power   |

## Scored quantities

| Quantity        | Description                                          |
|-----------------|------------------------------------------------------|
| Dose            | Dose in transport medium (PMMA)                      |
| Dose inWater    | Dose-to-water (SPR correction applied) [verified]    |
| Dose inSi       | Dose-to-silicon (SPR correction applied) [verified]  |
| Fluence         | All particles                                        |
| Fluence Primary | Primary protons only (filter applied)                |
| DLET            | Dose-averaged LET in transport medium (PMMA)         |
| DLET inWater    | Dose-averaged LET in water                           |
| DLET inSi       | Dose-averaged LET in silicon                         |
| TLET            | Track-averaged LET in transport medium (PMMA)        |
| TLET inWater    | Track-averaged LET in water                          |
| TLET inSi       | Track-averaged LET in silicon                        |
| DQEFF           | Dose-averaged (z\_eff/β)²                            |
| TQEFF           | Track-averaged (z\_eff/β)²                           |

Output in both TEXT and BDO formats.
