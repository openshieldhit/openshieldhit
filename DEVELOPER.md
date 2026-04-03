# Developer Notes for OpenShieldHIT

This document collects internal development guidelines and decisions for OSH.
Not intended for end-users.

## Style

- Use `struct foo` explicitly; **do not use `typedef struct`** to hide the fact that it is a struct.
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

### Design Decisions

This project is also intended as a pedagogical reference. Where practical,
non-obvious design choices should be briefly explained in a comment near the
decision — struct layout, ownership model, unit conventions, algorithm
selection, and similar. The goal is that a reader can understand not just
*what* the code does but *why*, without having to dig through git history.

### General Coding Conventions
- Public functions have `osh` prefix : `osh_foobar()`
- Private/internal functions use `static` and no `_osh` prefix : `foobar()`
- Private/inline functions use `static inline` and `_` prefix: `_foobar()`. This naming may be subject to change later when we settle on a better convention.
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

## Porting Notes

- This is a full reimplementation. No direct code copied from SHIELD-HIT12A.

## Build

- Uses CMake with a modular subdirectory layout.
- Unit tests are built using `CTest`.
