# Contributing

See `CONTRIBUTING.md` in the repository root for the full contributor guide.

Key points:

- All C code follows the style defined in `DEVELOPER.md` (clang-format enforced in CI)
- Every public function must have a Doxygen doc block (`@brief`, `@param`, `@returns`)
- New features need a unit test under `tests/unit/` and ideally an end-to-end case under `tests/cases/`
- Open an issue first for non-trivial changes; PRs without a linked issue may wait longer for review
