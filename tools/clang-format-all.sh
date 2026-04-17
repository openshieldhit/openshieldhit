#!/usr/bin/env bash
# Run this script from the project root (or any subdirectory inside this repo).
set -euo pipefail

if ! command -v clang-format >/dev/null 2>&1; then
    echo "error: clang-format not found in PATH" >&2
    exit 1
fi

repo_root="$(git rev-parse --show-toplevel 2>/dev/null || true)"
if [[ -z "${repo_root}" ]]; then
    echo "error: not inside a git repository" >&2
    exit 1
fi

cd "${repo_root}"

before="$(mktemp)"
after="$(mktemp)"
files_list="$(mktemp)"
trap 'rm -f "${before}" "${after}" "${files_list}"' EXIT

# Snapshot currently modified tracked C/H files before formatting.
git diff --name-only -- '*.c' '*.h' | sort > "${before}"

find . -type f \( -name '*.c' -o -name '*.h' \) \
    -not -path './build/*' \
    -not -path './build-*/*' \
    -not -path './_temp_shieldhit/*' \
    -print0 > "${files_list}"

if [[ ! -s "${files_list}" ]]; then
    echo "No .c/.h files found."
    exit 0
fi

echo "Running clang-format using .clang-format..."
xargs -0 clang-format -i -style=file < "${files_list}"

# Snapshot currently modified tracked C/H files after formatting.
git diff --name-only -- '*.c' '*.h' | sort > "${after}"

echo
echo "Files changed by this clang-format run:"
comm -13 "${before}" "${after}" | sed '/^$/d'

if [[ "$(comm -13 "${before}" "${after}" | wc -l | tr -d ' ')" -eq 0 ]]; then
    echo "(none)"
fi
