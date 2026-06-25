# Developer Notes for OpenShieldHIT

This document collects internal development guidelines and decisions for OSH.
Not intended for end-users.

## Style

- Use `struct foo` explicitly; **do not use `typedef struct`** to hide the fact that it is a struct.
- In general, discourage `typedef` for project-defined types. Prefer the explicit type kind in declarations (`struct foo`, `enum bar`) because it makes ownership, layout, and API boundaries easier to read when scanning code.
- For structs:
  - Largest / most-aligned fields first
  - Keep fields that are commonly read together near each other
  - Keep “rarely touched” caches separate

- Compilation follows the **C11 standard**, but with some **C90-inspired** stylistic ideas:
  - `//` single-line comments are permitted **only for temporary code**, such as quick notes or `TODO`s. They must not remain in production code. This helps spotting temporary and WIP code.
  - Use `stdint.h` types (e.g., `uint32_t`) only when needed.
  - Prefer `char` for boolean flag fields in structs rather than `<stdbool.h>` — it is explicit about storage size and keeps struct layout predictable. Use `int` for function return codes.
  - Use `<stdatomic.h>` for atomic operations on shared tallies/counters (e.g., `_Atomic`, `atomic_fetch_add`). Prefer lock-free atomics over mutexes for hot paths.
  - Use `_Alignas(64)` on per-thread data structures to avoid false sharing across cache lines.
  - Avoid `<threads.h>` for portability (MSVC does not implement it); use pthreads or OpenMP for the thread pool instead.

- Use `const` on pointers when the function does not modify the referenced data. This clarifies ownership and protects against accidental changes.
- Prefer `double const *p` rather than `const double *p`.

- Do **not** declare variables inside `for`, `if`, or other control blocks.
  Declare them at the beginning of the block (old K&R style).

### Initialization and Lifetime

- Prefer **declare early, define late**:
  - Variables should be declared at the beginning of the block, but assigned or initialized only at the point where their value first becomes logically known.
  - Avoid assigning placeholder values such as `0` or `NULL` merely to silence warnings; this can mask logic errors where a variable is used before being properly set.
  - Loop counters and temporary variables should typically be assigned immediately before the loop or computation that uses them, not at function entry.
  - The goal is that uninitialized use is detectable during review, testing, or with compiler warnings, rather than silently producing incorrect results.

- Do **not** initialize variables or structs at declaration unless required for correctness or safety.
  Avoiding unnecessary initialization helps the compiler or developer detect logic errors where values are used before being properly assigned.

- Acceptable reasons to initialize at declaration include:
  - A variable must have a defined value before first conditional use.
  - An API requires initialization (for example, pointers passed to functions that may free or overwrite them).
  - Defensive initialization clearly improves robustness without hiding logic errors.

### Documentation

Use Doxygen-style `/** ... */` block comments on all non-trivial functions,
placed immediately before the function signature (not inside the body).

Follow this structure:

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

Rules:
- `@brief` is a single line — the summary a reader sees in an overview.
- `@details` holds everything that needs more than one sentence.
- Use `@param[in]`, `@param[out]`, `@param[in,out]` as appropriate.
- Omit `@param` and `@returns` only for trivial flag-setter functions
  where the signature is completely self-explanatory.
- Inline `/* */` comments inside the function body are for
  non-obvious *implementation* details (unit conversions, index tricks),
  not for repeating what the Doxygen block already says.

### Design Decisions

#### No Allocations on the Simulation Hot Path

`malloc`, `calloc`, `realloc`, and `free` must not be called during simulation,
i.e. inside `osh_scoring_score_step()`, `osh_transport_step()`, or anything they
call. Per-step `calloc`/`free` of scratch buffers was found to dominate CPU time
(zeroing via `memset` accounted for ~60 % of all instructions in a CT transport
run). All scratch buffers must be pre-allocated at compile/setup time and stored
in the corresponding runtime struct (e.g. `osh_scoring_runtime.crossing_buf`).

#### API Stability

The public API is still early and **not yet frozen**.

During this stage, prefer the cleaner interface over preserving temporary
compatibility glue. Small or medium API changes are acceptable when they
improve ownership, readability, naming, or module boundaries. Once the core
architecture settles, API changes should become much more deliberate.

This project is also intended as a pedagogical reference. Where practical,
non-obvious design choices should be briefly explained in a comment near the
decision — struct layout, ownership model, unit conventions, algorithm
selection, and similar. The goal is that a reader can understand not just
*what* the code does but *why*, without having to dig through git history.

### General Coding Conventions
#### Return Value Conventions

Two distinct return value conventions are used, and they must not be mixed:

**Predicate functions** return `int` with the meaning "yes/no":
- `1` = true / condition holds
- `0` = false / condition does not hold
- Use `if (foo(...))` at call sites.
- Typically named as questions: `_in_body()`, `_in_zone()`, `osh_gemca2_check_surface_side()`.

**Operation functions** return `enum osh_status` (defined in `common/osh_rc.h`):
- `OSH_OK = 0` = success
- Any other value = specific failure reason (see `osh_rc.h`)
- Use `if (foo(...) != OSH_OK)` at call sites.
- Note that `OSH_OK = 0` means success, which is the **opposite** of predicate convention — do not confuse the two.

Functions that return actual computed values (zone IDs, distances, counts) use their natural type (`size_t`, `double`, etc.) and document their sentinel/error value in the Doxygen `@returns` field.

If a function is trivially infallible (e.g. a pure memcpy wrapper), prefer `void` over a dummy `return 1`.

---

- Public functions have `osh` prefix : `osh_foobar()`
- Private/internal functions use `static` and no `_osh` prefix : `foobar()`
- Private/inline functions use `static inline` and `_` prefix: `_foobar()`. This naming may be subject to change later when we settle on a better convention.
- Prefer a top-down file layout:
  - In public headers, place forward declarations and the public API near the top so opening the file gives an immediate module overview.
  - In `.c` files, place exported/top-level functions before lower-level helpers. Use `static` prototypes, handler typedefs, and dispatch-table declarations near the top when needed to support this layout.
  - Order public data structures from high-level concepts to supporting detail where the type dependencies allow it.
- Do not add helper functions merely to shorten one call site. First look for an
  existing function in the module that owns the concept. Add a new helper only
  when it has a clear owner, captures reusable non-obvious logic, or prevents a
  real consistency bug across call sites.
- Prefer `enum` instead of `#define` numbered lists.
- Avoid magic numbers in implementation code:
  - Use named constants (`enum`, `static const`, or `#define` where appropriate) for option IDs, buffer sizes, exit codes, default filenames, and similar values.
  - Small arithmetic literals directly tied to local logic (e.g. `+1` for NUL terminator sizing) are acceptable when obvious and local.

- Public headers must use include guards of the form:

```c
#ifndef OSH_FOO_H
#define OSH_FOO_H
...
#endif /* OSH_FOO_H */
```

- The idea is to provide a public API as well, allowing these header files to be linked against C++, so encapsulate headers in:

```c
#ifdef __cplusplus
extern "C" {
#endif
```
and the matching closing brace:
```c
#ifdef __cplusplus
}
#endif
```


## Layout

- All internal headers are located alongside `.c` files in `src/*`.
- `include/` is reserved for public headers only, can be added later.
- No use of `include/internal/` or similar.
- Prefer a consistent per-module layout where it fits the module's lifecycle:

```text
src/<module>/
src/<module>/parse/
src/<module>/runtime/
src/<module>/README.md
```

- Intended meaning of these layers:
  - `src/<module>/` holds the module's main API, shared domain types, and setup-facing entry points.
  - `src/<module>/parse/` holds raw input parsing and parse-only helpers.
  - `src/<module>/runtime/` holds the compiled simulation-ready representation of parsed/setup data.
  - `src/<module>/README.md` should briefly describe what the module owns, what `parse/` means there, and what `runtime/` means there.

- Dependency direction rules:
  - `parse/` may depend on `common/` and its own module headers, but should avoid transport/runtime-specific coupling.
  - `runtime/` may depend on `common/` and other modules' public headers when needed for simulation-time integration.
  - `transport/` should depend on module `runtime/` layers rather than owning other modules' preparation code.
  - `runtime/` layers must not depend on `transport/`.
  - One module's `runtime/` layer should not reach into another module's `runtime/` internals unless that relationship is explicitly intended and documented.

- Do not create `parse/` or `runtime/` mechanically for every module. Use the structure when the module actually has those phases.

## Portability — banned POSIX-only APIs

Windows is a target platform. The following POSIX-only APIs must **not** appear
in any source or test file:

| Banned | Reason | Replacement |
|---|---|---|
| `<strings.h>` / `strcasecmp` | Not available on MSVC | Lowercase strings at storage time with `osh_lower_inplace()`, then use `strcmp` |
| `mkdtemp`, `mkstemp` | POSIX temp-file helpers | Write output to the current working directory with a fixed filename |
| `getpid()` | POSIX process ID | Not needed; use fixed filenames in tests |
| `<unistd.h>` | POSIX umbrella header | Use C standard library equivalents (`<stdio.h>`, `<stdlib.h>`, `<string.h>`) |
| `<sys/stat.h>`, `<sys/types.h>` | POSIX filesystem headers | Avoid; write to the current directory instead |
| `<threads.h>` | C11 threads (not implemented in MSVC) | Use pthreads or OpenMP on non-Windows; guard with `#if !defined(_WIN32)` if needed |

**Pattern for case-insensitive keyword matching:**
Lowercase the string once at parse/storage time using `osh_lower_inplace()` (from
`common/osh_readline.h`), then compare with lowercase string literals using plain
`strcmp`. Do not use `strcasecmp` at comparison time.

**Pattern for test output files:**
Write generated files to `.` (the CTest working directory) with descriptive
fixed names. Clean them up with `remove()` at the end of the test.
Do not use `/tmp`, `mkdtemp`, or `getpid`-based naming.

## Porting Notes

- This is a full reimplementation. No direct code copied from SHIELD-HIT12A.

## Build

- Uses CMake with a modular subdirectory layout.
- Unit tests are built using `CTest`.
- Named presets are defined in `CMakePresets.json`.

### Presets

| Preset | Binary dir | Flags | Use for |
|---|---|---|---|
| `debug` | `build_debug/` | `-O0 -g` | Day-to-day development, sanitisers |
| `release` | `build/` | `-O3` | Fastest baseline benchmarking |
| `relwithdebinfo` | `build-rel/` | `-O3 -g` | Optimised benchmarking with symbols |
| `prof` | `build_prof/` | `-O3 -g -fno-omit-frame-pointer` | `perf record` / flamegraphs |

Configure (first time, or after `CMakeLists` changes):

```bash
cmake --preset release
```

Build everything:

```bash
cmake --build --preset release --parallel
```

Build one target:

```bash
cmake --build --preset release --parallel --target gemca_raycast_bench
```

Run tests:

```bash
ctest --preset debug
```

### AVX2 / SIMD acceleration

AVX2+FMA is detected automatically at configure time for every preset via
`check_c_compiler_flag`.  When the compiler supports `-mavx2 -mfma`,
`src/gemca/runtime/osh_gemca_runtime_avx2.c` is compiled with those flags and
linked in.  A runtime dispatch (`__builtin_cpu_supports("avx2")`) ensures the
binary still runs on CPUs without AVX2.

The configure output will tell you whether it was detected:

```
-- gemca_runtime: AVX2+FMA zone batch enabled
-- gemca_runtime: AVX2+FMA not available, using scalar fallback
```

You do **not** need a separate preset for AVX2.  The `release` preset is the
one to use for SIMD benchmarking.
