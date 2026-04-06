## Runtime / API Status

The repo has moved past the original bootstrap phase. We now have:

- [x] Internal CLI parser in `src/cli/`
- [x] Thin `src/main.c` that configures a public API context and calls the facade
- [x] Public header `include/openshieldhit/openshieldhit.h`
- [x] Opaque `openshieldhit_context_t`
- [x] Context lifecycle, deep-copying setters, version helpers, and `openshieldhit_last_error()`
- [x] `OPENSHIELDHIT_RUN_VALIDATE` mode for parse/validate-only execution
- [x] Geometry loading wired into `openshieldhit_run()`
- [x] Startup / validation summary output
- [x] Basic tests for the public API (`tests/unit/test_osh_main.c`)
- [x] Basic tests for the internal CLI parser (`tests/unit/test_osh_cli.c`)
- [x] `osh_beam.h` struct layout finalised (beam_spot, beam_shared, beam_phsp, beam_workspace)
- [x] `osh_beam_setup_from_path(path, lg, wb_out)` constructor pattern established
- [x] Path normalisation (`osh_path_normalize`) and dirname (`osh_path_dirname`) in `osh_file.h`
- [x] Windows path separator handling centralised in `main.c` (no per-subsystem `#ifdef`)

## Beam loader — next steps

The struct layout and constructor API for the beam subsystem are finalised.
The parser exists but is not yet connected to the new API. Work to do in order:

- [ ] Fix `osh_beam_setup_from_path()` body — currently still calls `osh_beam_shared_init(wb->shared)` on a NULL pointer (shared is a pointer, not embedded value); decide: embed `beam_shared` by value in `beam_workspace` or allocate it
- [ ] Fix parser references to `beam->spot0` (field does not exist) → use `wb->spots[0]`
- [ ] Fix parser references to fields that moved (e.g. `sad`, `focus`, `use_div`, `use_sad` now in `beam_shared`; `tcut`, `pcut` moved to `beam_workspace`)
- [ ] Implement post-parse init step in `osh_beam_setup_from_path()`:
  - Derive `wdir` from path (done) and store `fname`
  - Convert p0/psigma from t0/tsigma using particle mass (relativistic)
  - Compute `shared->emax` and `shared->pmax` over all spots
  - Build `_tm[16]` rotation+translation matrix per spot from `shared->theta/phi` and `spot->p[]`
- [ ] Declare `osh_relative_path_to_file()` in `osh_file.h` (currently used in parser but undeclared)
- [ ] Wire `osh_beam_setup_from_path()` into `openshieldhit_run()` in `src/openshieldhit.c`
- [ ] Apply CLI override: after beam load, overwrite `wb->nstat` if `cfg->has_nstat` is set

## Material Runtime Design

The material parser should remain a raw input layer. Before transport starts,
the parsed `material_workspace` needs a material assembly step that builds a
cache-friendly runtime material database.

- [x] Parse named `MATERIAL` blocks with dense internal material indices
- [x] Reserve material index 0 for `blackhole` and index 1 for `vacuum`
- [x] Parse explicit composition cards and ICRU references as alternative composition sources
- [x] Store scalar user overrides: `RHO`, `STATE`, material/element mean excitation, `LOADDEDX`
- [ ] Import or regenerate the ICRU material database from `_temp_shieldhit/material`
- [ ] Add a material assembly layer after all input files are parsed
- [ ] Expand ICRU references into composition, density, state, and mean excitation defaults
- [ ] Preserve explicit scalar user overrides when filling unset ICRU/default properties
- [ ] Derive complementary composition fields: atom counts, mass fractions, number densities, electron densities
- [ ] Resolve GEMCA `zone->material_name` to dense `zone->material_idx`

Runtime representation should separate atomic transport data from nuclear target
sampling data. Atomic transport can then work before the fragmentation generator
exists.

- [ ] Build atomic transport tables for each material/projectile pair:
  - stopping power / dE/dx
  - range
  - optical depth
- [ ] Keep transport tables in dense, cache-friendly arrays indexed by material, projectile, and energy grid
- [ ] Build a separate nuclear target-sampling table per material from elemental/isotopic composition
- [ ] Use the nuclear table only when a nuclear interaction is sampled, then pass the sampled target to the future fragmentation generator
- [ ] Add a projectile table registry with a default ion set and dense projectile indices
- [ ] Support lazy extension of material/projectile tables at safe setup or batch boundaries, not inside the hot stepping loop

CT/voxel geometry needs an extra calibration layer. `RHO` in `MATERIAL` should
be treated as the default/reference density, while voxel geometry may provide a
local density per step from the CT image.

- [ ] Keep `HU -> material_idx` segmentation in the voxel/CT geometry layer, not in `struct material`
- [ ] Keep `HU -> rho` density calibration in the voxel/CT geometry layer, not in `struct material`
- [ ] Represent `HU -> rho` and `HU -> WEPL` calibrations as piecewise-linear tables with precomputed coefficients
- [ ] Represent `HU -> material_idx` as piecewise-constant HU bins mapped to dense material indices
- [ ] Prefer mass-normalized transport tables so local density scaling is cheap: mass stopping power, mass range, and mass interaction coefficients where applicable
- [ ] Let the transport step carry or query both `material_idx` and local `rho`; constant-density geometry uses the material default

## Current High-Priority TODO

- [ ] Wire beam loading into `openshieldhit_run()`
- [ ] Wire material loading into `openshieldhit_run()`
- [ ] Wire detect/scoring loading into `openshieldhit_run()`
- [ ] Implement `OPENSHIELDHIT_RUN_NORMAL`
  - Current state: returns `OPENSHIELDHIT_STATUS_NOT_SUPPORTED`.
- [ ] Add a validate-mode integration test using `tests/fixtures/test01/`
- [ ] Add tests for `openshieldhit_last_error()`
  - [x] Empty on fresh context
  - [x] Empty on NULL context (no crash)
  - [x] Set on unsupported run mode
  - [ ] Error reset on a subsequent successful call
  - [ ] Set on missing geometry file
- [ ] Add tests for context configure ownership semantics

## Packaging / Install TODO

- [ ] Decide whether the public facade should ship as static-only, shared-only, or both
- [ ] Install the public library target and headers under `include/openshieldhit/`
- [ ] Keep internal headers and the internal CLI parser non-installed

## Cleanup / Follow-up

- [ ] Make gemca parser library-safe (non-terminating)
  - Current state: `osh_gemca_parse()` calls `osh_error()` / `exit()` on failures, bypassing the `OPENSHIELDHIT_STATUS_PARSE_ERROR` return path. The `03_malformed_geometry` test exits with code 78 by accident (via `exit(EX_CONFIG)` inside the parser). See comment in `src/gemca/osh_gemca2.c`.
- [ ] Make beam parser library-safe (non-terminating)
  - Same pattern as gemca: all `osh_err()` call sites in `src/beam/osh_beam_parse.c` call `exit()`.
- [ ] Improve validation diagnostics — include file and line info from parsers
- [ ] Wire `osh_beam_print` / `osh_beam_print_spot` to logger instead of stdout
- [ ] Revisit public API once beam/material/detect loading are wired

## Notes

- The public API is context-based and opaque; the CLI parser remains internal in `src/cli/`.
- `--dry-run` CLI maps to `OPENSHIELDHIT_RUN_VALIDATE` in the public API.
- `osh_path_normalize()` must be called in `main.c` on all CLI-supplied paths before any library call — library code assumes forward slashes throughout.
- Beam energies/momenta are always TOTAL (MeV, MeV/c) — never per nucleon. Transport engine divides by `part->a` at table lookup time.
- nstat override pattern: beam file sets default, CLI `--nstat` wins if `cfg->has_nstat` is set.
