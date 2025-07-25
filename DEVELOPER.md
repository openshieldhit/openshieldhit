# Developer Notes for OpenShieldHIT

This document collects internal development guidelines and decisions for OSH.
Not intended for end-users.

## Style

- Use `struct foo` explicitly; **do not use `typedef struct`** to hide the fact that it is a struct.
- All code must adhere to **C90**, with the following *exceptions*:
  - `//` single-line comments are permitted **only for temporary code**, such as quick notes or `TODO`s. They must not remain in production code.
  - `stdint.h` may be used for fixed-width integers (e.g., `uint32_t`)
- Use `const` on pointers when the function does not modify the referenced data.
  This clarifies ownership and protects against accidental changes.
- Do **not** declare variables inside `for`, `if`, or other control blocks — declare them at the beginning of the block.
- Private/internal functions use `static` and a `_name()` prefix.
- Public headers must use include guards of the form:
  ```c
  #ifndef OSH_FOO_H
  #define OSH_FOO_H
  ...
  #endif /* OSH_FOO_H */
  ```
  - Do **not** by default initialize all variable upon declaration. This helps the compiler or developer detect logic errors where values are used before being properly assigned.

## Layout

- All internal headers are located alongside `.c` files in `src/*`.
- `include/` is reserved for public headers only.
- No use of `include/internal/` or similar.

## Porting Notes

- This is a full reimplementation. No direct code copied from SH12A.

## Build

- Uses CMake with a modular subdirectory layout.
- Unit tests are built using `CTest`.


