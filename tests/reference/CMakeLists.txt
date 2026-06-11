# ---- Reference benchmarks: IDD comparison against external codes ----
#
# Each idd_* subdirectory is a self-contained benchmark case:
#
#   idd_<name>/
#     geo.dat, beam.dat, mat.dat, detect.dat   (inputs, full transport)
#     reference/
#       idd_<code>.dat      (external reference curves: SH12A, FLUKA, ...)
#     compare_args.cmake    (optional; per-case tolerance overrides)
#     README                (setup description for reproducing the reference)
#
# Tests carry the label "reference" (excluded from the unit/cases labels;
# run with `ctest -L reference`).  A case with no reference curve yet is
# reported as SKIPPED, so the infrastructure can precede the data.
#
# The shieldhit/ and topas/ subdirectories hold mirrored input decks for the
# external codes; they are not registered as tests.

file(GLOB IDD_CASE_DIRS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/idd_*/")

foreach(case_dir ${IDD_CASE_DIRS})
    if(NOT IS_DIRECTORY "${case_dir}")
        continue()
    endif()

    get_filename_component(case_name "${case_dir}" NAME)

    add_test(
        NAME "reference::${case_name}"
        COMMAND ${CMAKE_COMMAND}
            -DOSH_EXECUTABLE=$<TARGET_FILE:openshieldhit>
            -DCASE_DIR=${case_dir}
            -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/${case_name}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/run_idd.cmake
    )

    set_tests_properties("reference::${case_name}" PROPERTIES
        LABELS "reference"
        SKIP_REGULAR_EXPRESSION "SKIPPED"
        TIMEOUT 1800
    )
endforeach()
