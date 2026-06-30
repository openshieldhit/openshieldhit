# fetch_fixtures.cmake — download the DICOM test fixtures from a GitHub Release.
#
# The DCPT head-phantom CT dataset (182 *.dcm files, ~267 MB uncompressed) used
# to live in Git LFS, but the repository ran out of LFS bandwidth (see #187).
# The dataset now ships as a release asset and is fetched on demand.
#
# Usage (standalone — works on Linux, macOS and Windows):
#
#     cmake -P cmake/fetch_fixtures.cmake
#
# Optionally override the destination directory (defaults to
# <repo>/tests/fixtures/dicom):
#
#     cmake -DFIXTURES_DIR=/some/path -P cmake/fetch_fixtures.cmake
#
# This script relies only on CMake built-ins — file(DOWNLOAD) and
# `cmake -E tar` (libarchive is bundled with CMake) — so it behaves identically
# on every CI platform and needs neither curl/wget nor a system tar.

cmake_minimum_required(VERSION 3.10)

# ---- Release coordinates (single source of truth) -------------------------
set(FIXTURE_TAG     "test-fixtures-v1")
set(FIXTURE_ARCHIVE "DCPT_headphantom.tar.gz")
set(FIXTURE_URL
    "https://github.com/openshieldhit/openshieldhit/releases/download/${FIXTURE_TAG}/${FIXTURE_ARCHIVE}")
set(FIXTURE_SHA256  "5ab7ac365ccbbf4e5d3eb915b30f8a123fce3d460b05c7a17ce4ce9d64707095")

# ---- Resolve paths relative to this script, not the caller's CWD ----------
get_filename_component(_repo_root "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)
if(NOT DEFINED FIXTURES_DIR)
    set(FIXTURES_DIR "${_repo_root}/tests/fixtures/dicom")
endif()

set(_dataset_dir "${FIXTURES_DIR}/DCPT_headphantom")
set(_marker      "${_dataset_dir}/.fixtures-${FIXTURE_TAG}.stamp")

# ---- Idempotent: skip if this exact version is already unpacked -----------
if(EXISTS "${_marker}")
    message(STATUS "DICOM fixtures already present (${FIXTURE_TAG}); skipping download.")
    return()
endif()

file(MAKE_DIRECTORY "${FIXTURES_DIR}")
set(_archive "${FIXTURES_DIR}/${FIXTURE_ARCHIVE}")

# ---- Download (hash-verified) ---------------------------------------------
message(STATUS "Downloading DICOM fixtures: ${FIXTURE_URL}")
file(DOWNLOAD "${FIXTURE_URL}" "${_archive}"
     EXPECTED_HASH SHA256=${FIXTURE_SHA256}
     TLS_VERIFY ON
     SHOW_PROGRESS
     STATUS _dl_status)
list(GET _dl_status 0 _dl_code)
list(GET _dl_status 1 _dl_msg)
if(NOT _dl_code EQUAL 0)
    file(REMOVE "${_archive}")
    message(FATAL_ERROR
        "Failed to download DICOM fixtures from ${FIXTURE_URL}\n  ${_dl_msg}")
endif()

# ---- Extract --------------------------------------------------------------
# Drop any stale extraction first so a re-fetch starts clean.
file(REMOVE_RECURSE "${_dataset_dir}")
message(STATUS "Extracting ${FIXTURE_ARCHIVE} into ${FIXTURES_DIR}")
execute_process(
    COMMAND "${CMAKE_COMMAND}" -E tar xzf "${_archive}"
    WORKING_DIRECTORY "${FIXTURES_DIR}"
    RESULT_VARIABLE _tar_result)
file(REMOVE "${_archive}")
if(NOT _tar_result EQUAL 0)
    message(FATAL_ERROR "Failed to extract ${_archive} (tar exit ${_tar_result})")
endif()

if(NOT IS_DIRECTORY "${_dataset_dir}")
    message(FATAL_ERROR
        "Expected ${_dataset_dir} after extraction, but it is missing — "
        "the archive layout may have changed.")
endif()

# ---- Stamp so subsequent configures/fetches are no-ops --------------------
file(WRITE "${_marker}"
    "DCPT_headphantom DICOM fixtures\n"
    "tag:    ${FIXTURE_TAG}\n"
    "source: ${FIXTURE_URL}\n"
    "sha256: ${FIXTURE_SHA256}\n")
message(STATUS "DICOM fixtures ready in ${_dataset_dir}")
