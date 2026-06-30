# ---- DICOM test fixtures (fetched from a GitHub release, not stored in git) ----
# The DCPT head-phantom dataset used by the DICOM tests is downloaded on demand
# (see cmake/fetch_fixtures.cmake and tests/fixtures/dicom/README.md). Fetching
# at configure time keeps it automatic for developers and CI alike; the script
# is idempotent so repeat configures don't re-download. Turn OFF to build/run
# the non-DICOM tests offline.
option(OSH_FETCH_DICOM_FIXTURES
    "Download the DCPT head-phantom DICOM test fixtures at configure time" ON)

if(OSH_FETCH_DICOM_FIXTURES)
    execute_process(
        COMMAND "${CMAKE_COMMAND}" -P "${PROJECT_SOURCE_DIR}/cmake/fetch_fixtures.cmake"
        RESULT_VARIABLE _osh_fetch_fixtures_result)
    if(NOT _osh_fetch_fixtures_result EQUAL 0)
        message(WARNING
            "Could not fetch DICOM test fixtures; DICOM tests will fail until they "
            "are present. Fetch manually with:\n"
            "    cmake -P cmake/fetch_fixtures.cmake")
    endif()
endif()

add_subdirectory(unit)
add_subdirectory(cases)
add_subdirectory(reference)
add_subdirectory(bench)
