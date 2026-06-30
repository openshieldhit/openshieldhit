# ---- Integration cases: one CTest entry per case directory ----
#
# Each subdirectory is a self-contained test case:
#
#   <name>/
#     geo.dat, beam.dat, mat.dat, detect.dat   (inputs)
#     expected/
#       exit_code.txt   (optional; default 0)
#       stdout.txt      (optional; if present, output is compared exactly)
#       stderr.txt      (optional; if present, output is compared exactly)
#     args.cmake        (optional; overrides default CLI args for this case)
#     README            (optional; human-readable description)
#
# Binary reference files (*.bdo) are ignored by the driver — they will be
# used for numeric output validation once transport is implemented.
#
# To add a new case: create a new subdirectory with the required input files.
# No CMake changes are needed.

file(GLOB CASE_DIRS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*/")

foreach(case_dir ${CASE_DIRS})
    if(NOT IS_DIRECTORY "${case_dir}")
        continue()
    endif()

    get_filename_component(case_name "${case_dir}" NAME)

    add_test(
        NAME "cases::${case_name}"
        COMMAND ${CMAKE_COMMAND}
            -DOSH_EXECUTABLE=$<TARGET_FILE:openshieldhit>
            -DCASE_DIR=${case_dir}
            -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/${case_name}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/run_case.cmake
    )

    set_tests_properties("cases::${case_name}" PROPERTIES
        LABELS "cases"
    )
endforeach()

# The DICOM integration case needs the fetched head-phantom dataset; depend on
# the DICOM_FIXTURES setup fixture so ctest downloads it on demand, only when
# this case runs (see tests/CMakeLists.txt). The "dicom" label lets callers
# select/exclude it, e.g. `ctest -LE dicom`.
if(TEST "cases::05_dicom_simple")
    set_tests_properties("cases::05_dicom_simple" PROPERTIES
        FIXTURES_REQUIRED DICOM_FIXTURES)
    set_property(TEST "cases::05_dicom_simple" APPEND PROPERTY LABELS "dicom")
endif()
