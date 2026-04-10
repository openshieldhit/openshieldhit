# Contributing to OpenShieldHIT

Thank you for your interest in contributing.
This file covers the practical steps to get the code building and a patch submitted.
For architecture decisions and internal design rationale see [DEVELOPER.md](DEVELOPER.md).

## Prerequisites

| Tool | Minimum version | Notes |
|------|-----------------|-------|
| C compiler | C11 | GCC, Clang, or MSVC |
| CMake | 3.15 | |
| clang-format | 18+ | for formatting checks |
| clang-tidy | 18+ | for static analysis |
| libSDL2 | any | optional, only needed for example programs |

## Build

```sh
cmake -B build
cmake --build build
ctest --test-dir build
```

All three steps must pass before submitting a PR.  The CI matrix runs on
Ubuntu, macOS, and Windows — if you only have one platform available, the
others are covered by CI.

## Code style

Formatting is enforced by clang-format.  Run it before committing:

```sh
clang-format -i $(git diff --name-only HEAD | grep '\.[ch]$')
```

Or check without modifying:

```sh
clang-format --dry-run --Werror src/**/*.c src/**/*.h
```

Beyond formatting, the style rules that clang-format cannot enforce are in
[DEVELOPER.md](DEVELOPER.md).  The short version:

- C11, but declare variables at the top of each block (K&R style)
- No `typedef struct` — use `struct foo` explicitly
- `double const *p` not `const double *p`
- Block (`/* */`) comments only; `//` comments are reserved for temporary notes

## Submitting a patch

1. Fork the repository and create a feature branch off `main`.
2. Keep commits focused — one logical change per commit.
3. Make sure `ctest` passes locally.
4. Open a pull request against `main`.  The PR description should explain
   *why* the change is needed, not just what it does.
5. CI runs tests, clang-format, and clang-tidy automatically.  Fix any
   failures before requesting review.

## Reporting bugs

Open a GitHub issue.  Include the OS, compiler version, and a minimal
reproducer (ideally one of the `tests/cases/` input sets with a modified
parameter that triggers the bug).
