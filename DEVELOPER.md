# Developer Notes for OpenShieldHIT

This document is the authoritative source for internal coding rules and design
rationale in OpenShieldHIT. It is intentionally more than a checklist: many
contributors use it to learn the project style, so each rule should explain the
reason behind the choice when the reason is not obvious.

When reviewing code, prefer referring to the numbered rules below, for example
"violates §1" or "see §11". Tool-specific guidance files should point here
instead of repeating these rules.

## How to Read This Document

Each numbered section is a rule group. A group may contain:

- **Rule:** the required project convention.
- **Reason:** the design motivation, written in prose when that helps learning.
- **Exceptions:** allowed cases where the rule would otherwise be too blunt.
- **Examples:** small examples of preferred and avoided forms.

The rules are grouped by topic rather than by severity. A later section is not
less important than an earlier one.

Unless a rule explicitly says otherwise, apply these rules to new and touched
code. Prefer focused fixes over mechanical tree-wide churn: cleanup is welcome
when it is near the change being made or removes real confusion, but style-only
sweeps should be deliberate PRs.

## §1 C Dialect, Declarations, and Initialization

### §1.1 C Dialect and Declaration Placement

**Rule:** OpenShieldHIT is C11, but uses a K&R-inspired declaration style. Declare
variables at the beginning of the block. Do not declare variables inside `for`,
`if`, `while`, or other control statements.

### §1.2 Declare Early, Define Late

Prefer **declare early, define late**:

- Put declarations at the top of the block.
- Assign a value only when the value is logically known; do not initialize
  merely to make an uninitialized-use warning disappear.
- Avoid placeholder initializers such as `0` or `NULL` merely to silence warnings.
- Assign loop counters and temporary variables close to the loop or computation
  that uses them.

**Reason:** The project chooses this style because it keeps local state visible
when scanning the start of a block. That makes review easier and keeps the code
visually consistent with the plain-C style used throughout the codebase. This is
not a claim that modern C cannot support mixed declarations; it is a deliberate
project convention.

Avoiding unnecessary initialization is a bug-finding tool, not just style. A
variable that is accidentally used before it is assigned should look suspicious
to reviewers and should remain visible to compiler diagnostics; it must not
silently contain a placeholder value that lets the bug pass through review.

### §1.3 Initialization Exceptions

**Exceptions:** Some objects must be initialized at declaration for correctness
or clarity:

- `const` locals must receive their value when declared.
- An API may require a pointer or handle to start as `NULL` because cleanup code
  may free it or because an output function may overwrite it.
- Defensive initialization is acceptable when it clearly improves robustness and
  does not hide control-flow mistakes.

Use exceptions because they express correctness, not as a way to reintroduce
mixed declaration style.

**Examples:**

**Preferred:**

```c
size_t i;
double sum;

sum = 0.0;
for (i = 0u; i < n; ++i) {
    sum += values[i];
}
```

**Allowed exception:**

```c
double const scale = 1.0 / norm;
```

**Avoid:**

```c
double sum = 0.0;
for (size_t i = 0u; i < n; ++i) {
    sum += values[i];
}
```

## §2 Type Spelling, `const`, and Boolean Storage

### §2.1 Struct and Enum Spelling

**Rule:** Use `struct foo` and `enum bar` explicitly for project-defined types.
Do not hide structs behind `typedef struct`.

In general, avoid `typedef` for project-defined types. Prefer spelling out the
kind of type in declarations.

**Reason:** This is a readability choice. The code should show when a value is a
struct instead of hiding that fact behind an alias. Seeing
`struct osh_beam_workspace` or `enum osh_status` at a call site tells the reader
something useful about ownership, layout, and API boundaries. It also keeps
project-defined types from looking like built-in scalar types.

### §2.2 Pointer `const` Spelling

**Rule:** Use `const` on pointers when the function does not modify the
referenced data. Prefer:

```c
double const *p;
```

over:

```c
const double *p;
```

**Reason:** Both spellings are valid C, but the project standardizes on
right-binding `const` because it scales more predictably to complex pointer
declarations.

### §2.3 Boolean Storage and Exact-Width Integers

**Rule:** Prefer `char` for boolean flag fields inside structs rather than
`<stdbool.h>`. Use `int` for predicate return values.

**Reason:** `char` makes storage size explicit in structs and keeps layout
predictable. Predicate functions use `int` because that is the idiomatic C return
type for yes/no checks and interacts cleanly with `if (predicate(...))`.

Use `stdint.h` types such as `uint32_t` only when the exact width matters.

## §3 Struct Layout and Data Locality

### §3.1 Field Ordering

**Rule:** Lay out structs deliberately:

- Put largest or most-aligned fields first.
- Keep fields that are commonly read together near each other.
- Keep rarely touched caches, diagnostics, or optional state separate from hot
  fields when that improves locality.

**Reason:** Struct layout is part of performance and readability. A reader should
be able to see which fields are central to the object's identity and which are
supporting caches or diagnostics. Good layout also avoids accidental padding and
helps future parallel code keep hot state compact.

### §3.2 Cache-Line Alignment

Use `_Alignas(64)` on per-thread or per-worker data structures when false sharing
across cache lines is a real concern.

## §4 Comments and Doxygen Documentation

### §4.1 Production Comments

**Rule:** Production code uses block comments, not `//` comments. Single-line
`//` comments are reserved for temporary notes or work-in-progress markers and
must not remain in production code.

**Reason:** This makes temporary code visually distinct during review. If `//`
appears in a diff, it should attract attention.

### §4.2 Inline Declaration Comments

**Rule:** Encourage short inline comments on variable declarations when the name
or type alone does not make the role obvious.

Use this especially for local state in larger functions, scratch buffers, cached
values, ownership-sensitive pointers, and variables whose unit or lifetime is not
immediately clear.

**Reason:** Declaration comments make `.c` files much easier to scan. A reader
can understand the local working set of a function at the declaration block
before reading the control flow. This supports the same goal as §1: the beginning
of a block should give a useful overview of the state the code will manipulate.

**Preferred:**

```c
double *dist_batch;              /* Caller-owned boundary-distance scratch. */
size_t primaries_completed;      /* Histories finished after pool compaction. */
struct osh_diag_sink const *diag; /* Borrowed; NULL means silent. */
```

**Avoid:**

```c
size_t i; /* i */
double sum; /* double sum */
```

### §4.3 Function Documentation

**Rule:** Use Doxygen-style `/** ... */` block comments on non-trivial functions.
Place the comment immediately before the function signature.

Use this structure:

```c
/**
 * @brief One-sentence summary of what the function does.
 *
 * @details
 * Expanded explanation: sign conventions, units, algorithm choice,
 * design rationale, or anything a reader cannot infer from the code alone.
 * Use @details whenever the brief is not self-contained.
 *
 * @param[in]     foo  Description.
 * @param[in,out] bar  Description.
 * @param[out]    baz  Description.
 *
 * @returns Description of return value, or "Does not return" if the
 *          function aborts via a noreturn fatal helper on all error paths.
 */
```

**Details:**

- `@brief` is a single line: the summary a reader sees in an overview.
- `@details` holds explanations that need more than one sentence.
- Use `@param[in]`, `@param[out]`, and `@param[in,out]` accurately.
- Omit `@param` and `@returns` only for trivial flag-setter functions where the
  signature is completely self-explanatory.
- Inline `/* */` comments inside a function body are for non-obvious
  implementation details, such as unit conversions or indexing choices. They
  should not repeat what the Doxygen block already says.

**Reason:** OpenShieldHIT is also a pedagogical codebase. Documentation should
help readers understand not only what a function does, but why the implementation
chooses a particular convention, algorithm, unit, or ownership model.

### §4.4 Long-Form Documentation Drift

**Rule:** Before merging a PR to `main`, update applicable long-form
documentation when the change affects invariants, user-visible behavior, module
ownership, file formats, workflows, or non-obvious design decisions.

Check `docs/`, `src/<module>/README.md`, and other subdirectory READMEs that
own the changed area. If none of them are affected, no documentation edit is
needed, but the review should make that choice deliberately.

**Reason:** Code comments explain local implementation choices; long-form docs
explain the model future contributors rely on. Keeping them current at merge
time prevents `docs/` and module READMEs from drifting into historical notes.

## §5 Return Value Conventions

### §5.1 Predicates

**Rule:** Predicate functions return `int`:

- `1` means true or condition holds.
- `0` means false or condition does not hold.
- Use `if (foo(...))` at call sites.
- Predicate names should read like checks or questions where practical, for
  example `_in_body()`, `_in_zone()`, or `osh_gemca2_check_surface_side()`.

### §5.2 Operations

**Rule:** Operation functions return `enum osh_status`:

- `OSH_OK = 0` means success.
- Any other value is a specific failure reason.
- Use `if (foo(...) != OSH_OK)` at call sites.

**Reason:** These two conventions intentionally differ. Predicate truth uses
non-zero as true, while operation success uses zero through `OSH_OK`. Mixing them
is a common C bug, so call sites should make the convention obvious.

Functions that return computed values, such as zone IDs, distances, or counts,
use their natural type (`size_t`, `double`, and so on) and document any sentinel
or error value in the Doxygen `@returns` field.

### §5.3 Infallible Functions

If a function is trivially infallible, such as a pure copy helper, prefer `void`
over a dummy success return.

## §6 Naming and File Organization

### §6.1 Function Naming

**Rule:** Public functions use the `osh_` prefix:

```c
osh_foobar()
```

Private file-local functions are `static` and do not use an `_osh` prefix:

```c
foobar()
```

Private `static inline` helpers use a leading underscore for now:

```c
_foobar()
```

This `static inline` convention may change later if the project settles on a
better scheme.

**Reason:** The prefix marks API surface. File-local helpers should stay local
and readable rather than pretending to be public API.

### §6.2 File Layout

**Rule:** Prefer top-down file layout:

- In public headers, place forward declarations and public API near the top so
  opening the file gives an immediate module overview.
- In `.c` files, place the exported or top-level entry points near the top, in
  the order a reader should understand the module. Put lower-level auxiliary
  helpers closer to the bottom. Avoid bottom-up files where the main behavior is
  hidden after pages of support code.
- Use `static` prototypes, handler typedefs, and dispatch-table declarations
  near the top when needed to make top-down function order compile cleanly.
- Order public data structures from high-level concepts to supporting detail
  where type dependencies allow it.

**Reason:** Top-down layout lets a reader understand the module's public shape
and main behavior before diving into local mechanics.

## §7 Helpers, Constants, and Magic Numbers

### §7.1 Helper Functions

**Rule:** Do not add helper functions merely to shorten one call site. First look
for an existing function in the module that owns the concept, and scan the shared
utility modules before writing local arithmetic or geometry helpers. Examples:

- Use `src/common/osh_vect.h` for vector operations, affine transforms, dot/cross
  products, normalization, and rotations.
- Use `src/physics/osh_kinematics.h` for generic relativistic kinematics and
  direction-rotation utilities.
- Look in the owning module before adding a file-local helper; for example,
  parsing helpers belong near parsing code, scoring helpers near scoring code,
  and transport helpers near transport code.

Add a new helper only when it has a clear owner, captures reusable non-obvious
logic, or prevents a real consistency bug across call sites.

**Exception:** Hot-path helpers may be kept local, `static`, or `static inline`
when inlining is important for performance or when keeping the helper next to
the hot loop makes ownership clearer. This exception is for measured or
well-motivated hot paths, not for avoiding a search for existing utilities.

**Reason:** Helpers are abstraction, not decoration. A helper that does not own a
concept can make the call graph harder to understand. Reusing the existing
utility functions also keeps numerics and edge-case behavior consistent across
modules.

### §7.2 Conditional Operator

**Rule:** Do not introduce new conditional operators (`?:`) in production C.
Prefer explicit `if` / `else` control flow.

**Exception:** A `?:` may be acceptable inside a macro definition when the macro
genuinely must expand to a single expression, such as a local MIN/MAX-style
helper. Keep both arms simple, parenthesize macro arguments correctly, and avoid
side effects in the condition or either arm. Do not use nested conditional
operators.

Existing ternaries do not need a mechanical tree-wide cleanup. When touching
nearby code, replace them opportunistically, and treat new non-macro ternaries
or nested ternaries as review findings.

**Reason:** Readability comes first: explicit `if` / `else` blocks are easier to
scan, review, and debug than compact conditional expressions. Ternaries are also
easy to pack into larger expressions where evaluation order, units, casts, or
ownership become harder to see. OpenShieldHIT favors plain control flow for
maintainability and review.

### §7.3 Numbered Lists and Constants

**Rule:** Prefer `enum` instead of `#define` for numbered lists.

**Rule:** Avoid magic numbers in implementation code. Use named constants
(`enum`, `static const`, or `#define` where appropriate) for option IDs, buffer
sizes, exit codes, default filenames, and similar values.

Before introducing a new mathematical or physical constant, check
`include/openshieldhit/const.h` (and the internal compatibility include
`src/common/osh_const.h`). Use project constants such as `OSH_M_PI` and
`OSH_M_PI_180`; do not rely on non-standard platform macros such as `M_PI` from
`math.h`.

Particle identities and properties must also have one source of truth. Do not
redefine proton, neutron, electron, ion, or PDG constants locally. Use
`src/particle/osh_particle_const.h` for particle masses, `src/particle/osh_particle_pdg.h`
for PDG codes, and the accessors in `src/particle/osh_particle.h` for derived
particle properties such as nuclear or atomic mass. For example, use
`OSH_PART_MASS_PROTON` or `osh_particle_nuclear_mass_mev_from_za()` rather than a
local `938.272...` literal.

**Exceptions:** Small arithmetic literals directly tied to local logic, such as
`+1` for NUL terminator sizing, are acceptable when obvious and local.

## §8 Public Headers and C++ Consumers

### §8.1 Include Guards

**Rule:** Public headers use include guards of the form:

```c
#ifndef OSH_FOO_H
#define OSH_FOO_H
...
#endif /* OSH_FOO_H */
```

### §8.2 C++ Linkage

**Rule:** Public headers must be consumable from C++. Wrap public declarations in:

```c
#ifdef __cplusplus
extern "C" {
#endif
```

and close with:

```c
#ifdef __cplusplus
}
#endif
```

**Reason:** OpenShieldHIT provides a C API that should be usable from C++ code
without name-mangling surprises. Internal headers under `src/` should still keep
clean include boundaries, but `include/openshieldhit/` is the public contract.

## §9 Module Layout and Dependency Direction

### §9.1 Header Locations

**Rule:** Internal headers live alongside `.c` files under `src/`. Public headers
live under `include/openshieldhit/`. Do not create `include/internal/` or similar
parallel internal include trees.

### §9.2 Module Directory Shape

**Rule:** Prefer this per-module layout when it fits the module lifecycle:

```text
src/<module>/
src/<module>/parse/
src/<module>/runtime/
src/<module>/README.md
```

Meaning:

- `src/<module>/` holds the module's main API, shared domain types, and
  setup-facing entry points.
- `src/<module>/parse/` holds raw input parsing and parse-only helpers.
- `src/<module>/runtime/` holds the compiled simulation-ready representation of
  parsed or setup data.
- `src/<module>/README.md` should briefly describe what the module owns, what
  `parse/` means there, and what `runtime/` means there.

**Exceptions:** Do not create `parse/` or `runtime/` mechanically for every
module. Use those layers only when the module actually has those phases.

### §9.3 Dependency Direction

**Dependency direction rules:**

- `parse/` may depend on `common/` and its own module headers, but should avoid
  transport/runtime-specific coupling.
- `runtime/` may depend on `common/` and other modules' public headers when
  needed for simulation-time integration.
- `transport/` should depend on module `runtime/` layers rather than owning other
  modules' preparation code.
- `runtime/` layers must not depend on `transport/`.
- One module's `runtime/` layer should not reach into another module's runtime
  internals unless that relationship is explicitly intended and documented.

**Reason:** The project separates cold input/setup state from hot
simulation-ready state. The directory structure should make that lifecycle clear
instead of hiding parsing, compilation, and transport responsibilities in one
place.

## §10 No Allocations on the Simulation Hot Path

### §10.1 Allocation Rule

**Rule:** `malloc`, `calloc`, `realloc`, and `free` must not be called during
simulation hot paths, meaning inside `osh_scoring_score_step()`,
`osh_transport_step()`, or anything they call.

All scratch buffers must be pre-allocated at compile/setup time and stored where
hot-path code can reuse them.

**Reason:** Per-step allocation is expensive enough to dominate runtime. Earlier
scratch-buffer allocation showed zeroing via `memset` accounting for roughly
60 percent of all instructions in a CT transport run. Hot code should spend time
on transport and scoring, not allocator and zero-fill overhead.

### §10.2 Caller-Owned Scratch and Counters

Buffers that traversal mutates per step are caller-owned so they can become
per-worker without churn. The voxel-crossing scratch lives in
`struct osh_scoring_scratch`, is passed into `osh_scoring_score_step()`, and the
serial driver owns a long-lived copy pre-sized at compile time in
`osh_scoring_runtime.master_scratch` (available through
`osh_scoring_runtime_master_scratch()`). A parallel worker owns its own scratch.

The same caller-owned discipline applies to mutable counters on the hot path. The
transport profile (`struct osh_transport_profile`) is per-worker: a worker
accumulates phase timers and event counters into its own profile carried on
`osh_worker_context.profile` and threaded into `osh_transport_ion_step()`. Workers
must not race on a shared profile. The serial worker points its profile straight
at the run master, so serial values stay unchanged. A parallel driver gives each
worker a private profile and folds them afterwards with
`osh_transport_profile_merge()`. Counters and per-phase `*_s` values sum as
aggregate work; `total_s` combines by maximum because concurrent workers overlap
in wall-clock time.

### §10.3 Expensive Math on Hot Paths

**Rule:** Avoid repeated calls to expensive transcendental functions such as
`sin()`, `cos()`, `tan()`, `acos()`, `asin()`, `atan2()`, `log()`, `exp()`, and
`pow()` on hot paths when the value can be precomputed, tabulated, cached, or
rewritten with simpler arithmetic.

**Reason:** These functions are often much more expensive than additions,
multiplications, fused multiply-adds, or table lookups. In per-step transport,
geometry, scoring, or inner physics loops, even a single avoidable
transcendental call can become visible in profiles.

**Exceptions:** Use the mathematically correct function when the physics or
geometry requires it. The rule is to avoid unnecessary repeated evaluation, not
to replace clear and correct physics with fragile approximations. Good patterns
include computing `sin`/`cos` once at parse/setup time for fixed angles, carrying
both `cos_phi` and `sin_phi` from a sampling helper, using squared lengths when
only comparisons are needed, and documenting any approximation used in hot code.

**Example:** Avoid generic `pow()` for simple fixed exponents in hot code.

```c
y = x * sqrt(x); /* preferred for x^(3/2) */
```

instead of:

```c
y = pow(x, 1.5);
```

The `pow()` form is general and may cost substantially more CPU cycles than the
specialized expression.

## §11 Portability and Banned APIs

### §11.1 Banned APIs

**Rule:** Windows is a target platform. POSIX-only APIs and APIs missing from
MSVC must not appear in source or test files.

| Banned | Reason | Replacement |
|---|---|---|
| `<strings.h>` / `strcasecmp` | Not available on MSVC | Lowercase strings at storage time with `osh_lower_inplace()`, then use `strcmp` |
| `mkdtemp`, `mkstemp` | POSIX temp-file helpers | Write output to the current working directory with a fixed filename |
| `getpid()` | POSIX process ID | Not needed; use fixed filenames in tests |
| `<unistd.h>` | POSIX umbrella header | Use C standard library equivalents (`<stdio.h>`, `<stdlib.h>`, `<string.h>`) |
| `<sys/stat.h>`, `<sys/types.h>` | POSIX filesystem headers | Avoid; write to the current directory instead |
| `<threads.h>` | C11 threads, not implemented in MSVC | Use pthreads or OpenMP on non-Windows; guard with `#if !defined(_WIN32)` if needed |

**Reason:** Portability should be visible during implementation and review, not
discovered later in Windows CI.

### §11.2 Case-Insensitive Keyword Matching

**Pattern for case-insensitive keyword matching:** Lowercase the string once at
parse/storage time using `osh_lower_inplace()` from `common/osh_readline.h`, then
compare with lowercase string literals using plain `strcmp`. Do not use
`strcasecmp` at comparison time.

### §11.3 Test Output Files

**Pattern for test output files:** Write generated files to `.` (the CTest
working directory) with descriptive fixed names. Clean them up with `remove()` at
the end of the test. Do not use `/tmp`, `mkdtemp`, or `getpid`-based names.

### §11.4 Atomics and Threads

**Rule:** Use `<stdatomic.h>` for atomic operations on shared tallies and counters
when atomics are needed, for example `_Atomic` and `atomic_fetch_add`. Prefer
lock-free atomics over mutexes for hot paths.

**Rule:** Avoid `<threads.h>` for portability. Use pthreads or OpenMP for thread
pools where appropriate, with platform guards where needed.

## §12 API Stability and Pedagogical Comments

### §12.1 API Stability

**Rule:** The public API is still early and not yet frozen. During this stage,
prefer the cleaner interface over preserving temporary compatibility glue. Small
or medium API changes are acceptable when they improve ownership, readability,
naming, or module boundaries.

Once the core architecture settles, API changes should become much more
deliberate.

**Reason:** A young public API should not preserve mistakes too early. It is
better to converge on clear ownership and boundaries before the interface becomes
hard to change.

### §12.2 Pedagogical Comments

**Rule:** Explain non-obvious design choices near the decision when practical:
struct layout, ownership model, unit conventions, algorithm selection, dependency
direction, and similar.

**Reason:** OpenShieldHIT is intended as a pedagogical reference. A reader should
be able to understand not only what the code does, but why, without digging
through git history.

## §13 Source-Code Provenance

**Rule:** OpenShieldHIT shares no source code with SHIELD-HIT or SHIELD-HIT12A.
Do not describe it as a source-code port or derivative of those projects.

**Reason:** The connection is the application domain, compatible input concepts,
and scientific/algorithmic inspiration. Source-code provenance matters legally
and scientifically, so wording should not imply copied or derived code.

## §14 Build, Test, and Presets

### §14.1 Presets

**Rule:** The project uses CMake with a modular subdirectory layout. Named
configure and build presets are defined in `CMakePresets.json`.

| Preset | Binary dir | Flags | Use for |
|---|---|---|---|
| `debug` | `build_debug/` | `-Og -g3 -fno-omit-frame-pointer` | Day-to-day development |
| `asan` | `build_asan/` | `-Og -g3 -fsanitize=address,undefined` | Instrumented `ctest` (ASan + UBSan) |
| `release` | `build/` | `-O3` | Fastest baseline benchmarking |
| `relwithdebinfo` | `build-rel/` | `-O3 -g` | Optimized benchmarking with symbols |
| `prof` | `build_prof/` | `-O3 -g -fno-omit-frame-pointer` | `perf record` / flamegraphs |

> The `debug` preset does **not** enable sanitizers — it is a plain `-Og` build.
> Use the `asan` preset (or `-DOSH_ENABLE_SANITIZERS=ON` on any GCC/Clang
> configure) to build with AddressSanitizer + UndefinedBehaviorSanitizer; CI runs
> the full test suite under it on three toolchains — Linux GCC, Linux Clang, and
> macOS Apple Clang — so implementation-specific reports surface on each. Run
> leak checks separately: the `asan` CI job sets
> `ASAN_OPTIONS=detect_leaks=0` while the codebase's benign at-exit leaks are
> worked through, so LeakSanitizer does not fire there by default.

### §14.2 Common Commands

Configure, build, and test the debug build:

```bash
cmake --preset debug
```

```bash
cmake --build --preset debug --parallel
```

```bash
ctest --test-dir build_debug --output-on-failure
```

Configure and build the release build:

```bash
cmake --preset release
```

```bash
cmake --build --preset release --parallel
```

Build one target:

```bash
cmake --build --preset release --parallel --target gemca_raycast_bench
```

Configure, build, and run the test suite under AddressSanitizer + UBSan (matches
the `sanitizers` CI job; LeakSanitizer is disabled while benign at-exit leaks are
worked through):

```bash
cmake --preset asan
```

```bash
cmake --build --preset asan --parallel
```

```bash
ASAN_OPTIONS=detect_leaks=0 ctest --test-dir build_asan --output-on-failure
```

### §14.3 Formatting and Static Analysis

**Rule:** Respect the repository's local formatting and static-analysis
configuration. Use the checked-in clang-format and clang-tidy settings; do not
apply a personal formatter style or editor defaults that disagree with the repo.

Run formatting before submitting C or header changes:

```bash
./tools/clang-format-all.sh
```

**Reason:** Formatting and static-analysis settings are part of the project
contract. Keeping them local and version-controlled avoids style drift between
contributors, editors, operating systems, and LLM-generated patches.

## §15 Hot-Path Data Layout, SIMD, and Offload Readiness

### §15.1 SoA-Ready Hot Paths

**Rule:** Design simulation hot paths so they remain friendly to
structure-of-arrays (SoA) storage, vectorization, SIMD kernels, and future
offload backends where applicable.

Prefer APIs and data flow that operate on contiguous arrays or explicit batches
of the same field, rather than interleaving unrelated per-history state in a
single large object. Keep cold metadata separate from per-step mutable state.
When adding hot-path state, ask whether a future worker, SIMD lane, or GPU kernel
could consume it without unpacking an array-of-structs layout first.

**Reason:** The current serial code should not block future performance work.
OpenShieldHIT already uses pool-style transport state and batched geometry/scoring
operations; new hot code should preserve that direction. SoA-friendly design
makes cache behavior easier to reason about, enables explicit SIMD kernels, and
keeps a path open for thread, accelerator, or GPU offload work if a module later
justifies it.

**Examples:**

- Prefer separate arrays such as `x[]`, `y[]`, `z[]`, `ux[]`, `uy[]`, `uz[]`,
  `e[]` for live particle state over an array of large per-particle structs.
- Pass caller-owned scratch arrays into hot kernels instead of allocating or
  hiding scratch inside the callee.
- Keep immutable species/material metadata shared and separate from per-history
  mutable state.

### §15.2 Current AVX2 Dispatch

**Rule:** Do not add a separate preset just for AVX2. AVX2+FMA support is detected
automatically at configure time for every preset via `check_c_compiler_flag`.

When the compiler supports `-mavx2 -mfma`,
`src/gemca/runtime/osh_gemca_runtime_avx2.c` is compiled with those flags and
linked in. Runtime dispatch through `__builtin_cpu_supports("avx2")` keeps the
binary usable on CPUs without AVX2.

The configure output says which path was selected:

```text
-- gemca_runtime: AVX2+FMA zone batch enabled
-- gemca_runtime: AVX2+FMA not available, using scalar fallback
```

Use the `release` preset for SIMD benchmarking.
