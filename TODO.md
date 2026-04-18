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

## Beam loader

- [x] Fix `osh_beam_setup_from_path()` shared initialisation for embedded `beam_shared`
- [x] Fix parser references to moved fields (`sad`, `focus`, `use_div`, `use_sad`, `tcut`, `pcut`)
- [x] Derive `wdir` from the beam file path
- [x] Convert p0/psigma from t0/tsigma using particle mass (relativistic)
- [x] Compute `shared.emax` and `shared.pmax` over all spots
- [x] Build `_tm[16]` rotation+translation matrix per spot from `shared.theta/phi` and `spot->p[]`
- [x] Declare `osh_relative_path_to_file()` in `osh_file.h`
- [x] Wire `osh_beam_setup_from_path()` into `openshieldhit_run()` in `src/openshieldhit.c`
- [x] Apply CLI `--nstat` override after beam load when `cfg->has_nstat` is set
- [ ] Store the parsed beam file path in `beam_workspace::fname`
- [ ] Add ridge-modulator / ripple-filter support
- [ ] Add MCPL phase-space import

## Material runtime

- [x] Parse named `MATERIAL` blocks with dense internal material indices
- [x] Reserve material index 0 for `blackhole` and index 1 for `vacuum`
- [x] Parse explicit composition cards and ICRU references
- [x] Store scalar user overrides: `RHO`, `STATE`, mean excitation, `LOADDEDX`
- [x] Import ICRU material database from `_temp_shieldhit/material` / `_temp_libdedx`
- [x] Add material assembly layer after parsing
- [x] Expand ICRU references into composition, density, state, and mean excitation defaults
- [x] Preserve explicit scalar user overrides over ICRU defaults
- [x] Derive complementary composition fields: atom counts and mass fractions
- [x] Resolve GEMCA `zone->material_name` to dense `zone->material_idx` in validate-mode setup
- [x] Build atomic transport tables: stopping power / dE/dx, CSDA range
- [x] Add a projectile table registry with a default ion set and dense projectile indices
- [ ] Refresh material data with ICRU90 density, stopping-power, and mean-excitation updates
- [ ] Derive number densities and electron densities
- [ ] Build optical depth tables
- [ ] Build nuclear target-sampling table per material
- [ ] Support lazy extension of material/projectile tables at batch boundaries
- [ ] Add batch/SIMD runtime lookup helpers for stopping power and CSDA range
- [ ] Add batch inverse-range-to-energy lookup for residual-range transport steps

## Transport

- [x] Straight-line CSDA transport loop (`src/transport/osh_transport.c`)
- [x] DELTAE fractional energy-loss step criterion
- [x] Boundary-limited steps with exit energy from residual CSDA range
- [x] Blackhole and vacuum zone handling
- [x] Transport wired into `openshieldhit_run()` (end-to-end Bragg peak produces correct output)
- [x] Multiple Coulomb scattering (Moliere / random hinge)
- [x] Gaussian energy straggling
- [ ] Gaussian MCS mode
- [ ] Vavilov energy straggling
- [ ] Nuclear interactions and secondary particles
- [ ] Batch ion step phases around runtime lookup hot spots (`range_lookup`, inverse residual-range lookup)
- [ ] Measure SIMD benefit for mixed-material pools vs same-material/species micro-batches

## Scoring

- [x] Detect file parser (`src/scoring/parse/`)
- [x] Filter rules (Z, A, E, GEN, ID, NPRIM) with compiled comparison operators
- [x] Geometry types: Mesh, Cyl, Zone, Voxel, All
- [x] Quantity types: ENERGY, FLUENCE, DOSE, LETFLU, DLET, TLET, COUNT, NORMCOUNT, NKERMA, ALANINE, MCPL
- [x] Scoring runtime compilation (`src/scoring/runtime/osh_scoring_compile.c`)
- [x] Per-step scoring wired to transport (`osh_scoring_score_step`)
- [x] ASCII text output
- [x] BDO2019 binary output
- [x] Scoring wired into `openshieldhit_run()`
- [ ] Detect loading wired into `openshieldhit_run()` (currently bypassed; scoring uses inline setup)
- [ ] Cylindrical (R,Z) mesh scoring
- [ ] Zone scoring (GEMCA zone membership)
- [ ] Voxel scoring (CT grid via raytrace)
- [ ] LET scoring (DLET, TLET algorithms)
- [ ] Alanine detector response
- [ ] MCPL phase-space output

## Current High-Priority TODO

- [x] Wire beam loading into `openshieldhit_run()`
- [x] Wire material loading into `openshieldhit_run()`
- [x] Wire scoring into `openshieldhit_run()`
- [x] Implement `OPENSHIELDHIT_RUN_NORMAL` — end-to-end transport runs
- [ ] Wire detect/scoring loading via public API (currently bypassed in `openshieldhit_run()`)
- [ ] Add a validate-mode integration test using `tests/fixtures/test01/`
- [ ] Add tests for `openshieldhit_last_error()`
  - [x] Empty on fresh context
  - [x] Empty on NULL context (no crash)
  - [x] Set on unsupported run mode
  - [ ] Error reset on a subsequent successful call
  - [ ] Set on missing geometry file
- [ ] Add tests for context configure ownership semantics

## CT / Voxel geometry

- [ ] Add GEMCA voxel geometry support
- [ ] `HU -> material_idx` segmentation in the voxel/CT geometry layer
- [ ] `HU -> rho` density calibration as piecewise-linear table
- [ ] `HU -> WEPL` calibration table
- [ ] Prefer mass-normalized transport tables so local density scaling is cheap

## Packaging / Install

- [ ] Decide whether the public facade should ship as static-only, shared-only, or both
- [ ] Install the public library target and headers under `include/openshieldhit/`
- [ ] Keep internal headers and the internal CLI parser non-installed

## Geometry Runtime (GEMCA)

- [x] Flat runtime structs: surfaces[], bodies[], zones[] replacing pointer-linked AST
- [x] RPN instruction evaluator replacing recursive AST traversal (zone membership + distance)
- [x] GUARD_BODY fast-reject prepended to each zone at setup time
- [x] AVX2+FMA zone-batch path (`osh_gemca_runtime_get_zone_batch_avx2`) with runtime dispatch
- [ ] AVX2 batch path for `eval_distance` / `osh_gemca_runtime_get_distance_batch`: vectorise
      inner body-distance loop with `_mm256_min_pd` (Roth minpos) and `_mm256_blendv_pd`
      (CSG combine); advance 4 rays simultaneously in the step-loop.  See TODO in eval_distance
      and runtime/README.md.
- [ ] Flat `insns_flat[]` + `insn_begin[]` layout in `gemca_runtime` to close the GPU-migration
      gap (current `zones[j].insns` is a host heap pointer; see runtime/README.md).
- [ ] Jacobs voxel-traversal dispatch in `eval_distance` (current fallback: RPP distance)

## Cleanup / Follow-up

- [ ] Revisit data-module structure (`material/`, `particle/`, embedded tables)
- [ ] Make gemca parser library-safe (non-terminating)
  - Current state: `osh_gemca_parse()` calls `osh_error()` / `exit()` on failures. The `03_malformed_geometry` test exits with code 78 via `exit(EX_CONFIG)` inside the parser.
- [x] Make beam parser library-safe (non-terminating)
- [x] Remove POSIX-only APIs (`<strings.h>`, `mkdtemp`, `getpid`) — see `DEVELOPER.md` portability section
- [x] Dispatch tables in `src/scoring/parse/` — fixed redefinition pattern (forward decl + initialiser)
- [ ] Improve validation diagnostics — include file and line info from parsers
- [x] Wire `osh_beam_print` / `osh_beam_print_spot` to logger instead of stdout
- [ ] Revisit public API once beam/material/detect loading are fully wired
- [ ] Rename `struct ray` to `struct ray_v` throughout (see TODO in `osh_coord.h`)

## Notes

- The public API is context-based and opaque; the CLI parser remains internal in `src/cli/`.
- `--dry-run` CLI maps to `OPENSHIELDHIT_RUN_VALIDATE` in the public API.
- `osh_path_normalize()` must be called in `main.c` on all CLI-supplied paths before any library call — library code assumes forward slashes throughout.
- Beam energies/momenta are always TOTAL (MeV, MeV/c) — never per nucleon. Transport engine divides by `part->a` at table lookup time.
- nstat override pattern: beam file sets default, CLI `--nstat` wins if `cfg->has_nstat` is set.
- Material names are case-sensitive; `osh_material_by_name` uses `strcmp`. Use lowercase `"blackhole"` and `"vacuum"` in geometry files.
- Keywords stored by parsers (geometry kind, quantity, fileformat) are always lowercased at storage time via `osh_lower_inplace`; compare with lowercase literals.
