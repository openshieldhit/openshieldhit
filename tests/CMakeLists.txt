# ---- DICOM test fixtures (fetched from a GitHub release, not stored in git) ----
# The DCPT head-phantom dataset used by the DICOM tests is downloaded on demand
# (see cmake/fetch_fixtures.cmake and tests/fixtures/dicom/README.md).
#
# The ~27 MB download is deferred to *test* time — it is NOT done at configure or
# build time. It runs as a CTest setup fixture (named DICOM_FIXTURES) that CTest
# executes only when a test which requires it is about to run. So:
#   * `cmake --build ...` (building openshieldhit to just run it) never downloads.
#   * `ctest -R something_non_dicom` / `ctest -LE dicom` never downloads.
#   * Running the DICOM tests triggers the fetch once, beforehand.
# The fetch script is idempotent (a stamp file records a completed download), so
# repeat `ctest` runs re-check the stamp and skip re-downloading.
add_test(NAME fetch_dicom_fixtures
         COMMAND "${CMAKE_COMMAND}" -P "${PROJECT_SOURCE_DIR}/cmake/fetch_fixtures.cmake")
set_tests_properties(fetch_dicom_fixtures PROPERTIES
    FIXTURES_SETUP DICOM_FIXTURES
    LABELS "fixtures")

add_subdirectory(unit)
add_subdirectory(cases)
add_subdirectory(reference)
add_subdirectory(bench)
