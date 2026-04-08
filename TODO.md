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

The beam parser is wired into validate mode. Remaining items are follow-ups,
not blockers for the current material-runtime work:

- [x] Fix `osh_beam_setup_from_path()` shared initialisation for embedded `beam_shared`
- [x] Fix parser references to `beam->spot0`; parser writes through `wb->spots[0]`
- [x] Fix parser references to moved fields such as `sad`, `focus`, `use_div`, `use_sad`, `tcut`, and `pcut`
- [x] Derive `wdir` from the beam file path
- [ ] Store the parsed beam file path in `beam_workspace::fname`
- [x] Convert p0/psigma from t0/tsigma using particle mass (relativistic)
- [x] Compute `shared.emax` and `shared.pmax` over all spots
- [x] Build `_tm[16]` rotation+translation matrix per spot from `shared.theta/phi` and `spot->p[]`
- [x] Declare `osh_relative_path_to_file()` in `osh_file.h`
- [x] Wire `osh_beam_setup_from_path()` into `openshieldhit_run()` in `src/openshieldhit.c`
- [x] Apply CLI `--nstat` override after beam load when `cfg->has_nstat` is set
- [ ] Add ridge-modulator / ripple-filter support; both should share the same implementation path
- [ ] Add MCPL phase-space import

## Material Runtime Design

The material parser should remain a raw input layer. Before transport starts,
the parsed `material_workspace` needs a material assembly step that builds a
cache-friendly runtime material database.

- [x] Parse named `MATERIAL` blocks with dense internal material indices
- [x] Reserve material index 0 for `blackhole` and index 1 for `vacuum`
- [x] Parse explicit composition cards and ICRU references as alternative composition sources
- [x] Store scalar user overrides: `RHO`, `STATE`, material/element mean excitation, `LOADDEDX`
- [x] Import or regenerate the ICRU material database from `_temp_shieldhit/material` / `_temp_libdedx`
- [ ] Refresh material data with ICRU90 density, stopping-power, and mean-excitation updates from `_temp_libdedx`
- [x] Add a material assembly layer after material parsing
- [x] Expand ICRU references into composition, density, state, and mean excitation defaults
- [x] Preserve explicit scalar user overrides when filling unset ICRU/default properties
- [x] Derive complementary composition fields: atom counts and mass fractions
- [ ] Derive number densities and electron densities
- [x] Resolve GEMCA `zone->material_name` to dense `zone->material_idx` in validate-mode setup

Runtime representation should separate atomic transport data from nuclear target
sampling data. Atomic transport can then work before the fragmentation generator
exists.

- [ ] Add `src/physics/` for Bethe fallback and transport-table generation
- [ ] Finalise the runtime table contract in `src/transport/prepare/osh_transport_runtime.h`
  - dense row-major ordering: `[material][projectile][energy]`
  - keep `energy_idx` innermost for interpolation-friendly hot-path access
  - store mass stopping power and mass range first; leave nuclear/optical tables separate
- [ ] Build atomic transport tables for each material/projectile pair:
  - stopping power / dE/dx
  - range
  - optical depth
- [x] Add initial runtime table layout sketch in `src/transport/prepare/osh_transport_runtime.h`
- [ ] Keep transport tables in dense, cache-friendly arrays indexed by material, projectile, and energy grid
- [ ] Build a separate nuclear target-sampling table per material from elemental/isotopic composition
- [ ] Use the nuclear table only when a nuclear interaction is sampled, then pass the sampled target to the future fragmentation generator
- [ ] Add a projectile table registry with a default ion set and dense projectile indices
- [ ] Support lazy extension of material/projectile tables at safe setup or batch boundaries, not inside the hot stepping loop

CT/voxel geometry needs an extra calibration layer. `RHO` in `MATERIAL` should
be treated as the default/reference density, while voxel geometry may provide a
local density per step from the CT image.

- [ ] Add GEMCA voxel geometry support
- [ ] Keep `HU -> material_idx` segmentation in the voxel/CT geometry layer, not in `struct material`
- [ ] Keep `HU -> rho` density calibration in the voxel/CT geometry layer, not in `struct material`
- [ ] Represent `HU -> rho` and `HU -> WEPL` calibrations as piecewise-linear tables with precomputed coefficients
- [ ] Represent `HU -> material_idx` as piecewise-constant HU bins mapped to dense material indices
- [ ] Prefer mass-normalized transport tables so local density scaling is cheap: mass stopping power, mass range, and mass interaction coefficients where applicable
- [ ] Let the transport step carry or query both `material_idx` and local `rho`; constant-density geometry uses the material default

## Current High-Priority TODO

- [x] Wire beam loading into `openshieldhit_run()`
- [x] Wire material loading into `openshieldhit_run()`
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

- [ ] Revisit data-module structure (`material/`, `particle/`, embedded tables)
  - Current state: works, but atomic/isotopic/material reference data is split in a way that makes ownership and discovery harder than it should be.
  - Think through whether a future `src/data/atomic/` and `src/data/material/` split would make lookup helpers and embedded datasets easier to find.
- [ ] Make gemca parser library-safe (non-terminating)
  - Current state: `osh_gemca_parse()` calls `osh_error()` / `exit()` on failures, bypassing the `OPENSHIELDHIT_STATUS_PARSE_ERROR` return path. The `03_malformed_geometry` test exits with code 78 by accident (via `exit(EX_CONFIG)` inside the parser). See comment in `src/gemca/osh_gemca2.c`.
- [x] Make beam parser library-safe (non-terminating)
  - `src/beam/osh_beam_parse.c` now reports parse/allocation failures via return codes instead of process exit.
- [ ] Improve validation diagnostics — include file and line info from parsers
- [x] Wire `osh_beam_print` / `osh_beam_print_spot` to logger instead of stdout
- [ ] Revisit public API once beam/material/detect loading are wired

## Notes

- The public API is context-based and opaque; the CLI parser remains internal in `src/cli/`.
- `--dry-run` CLI maps to `OPENSHIELDHIT_RUN_VALIDATE` in the public API.
- `osh_path_normalize()` must be called in `main.c` on all CLI-supplied paths before any library call — library code assumes forward slashes throughout.
- Beam energies/momenta are always TOTAL (MeV, MeV/c) — never per nucleon. Transport engine divides by `part->a` at table lookup time.
- nstat override pattern: beam file sets default, CLI `--nstat` wins if `cfg->has_nstat` is set.
