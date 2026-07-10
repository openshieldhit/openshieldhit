# ---- Reference benchmarks: IDD comparison against external codes ----
#
# Each idd_* subdirectory is a self-contained benchmark case:
#
#   idd_<name>/
#     geo.dat, beam.dat, mat.dat, detect.dat   (inputs, full transport)
#     reference/
#       idd_<code>.dat      (external reference curves: SH12A, FLUKA, ...)
#     compare_args.cmake    (optional; per-case tolerance overrides)
#     README.md             (setup description for reproducing the reference)
#
# Tests carry the label "reference" (excluded from the unit/cases labels;
# run with `ctest -L reference`).  A case with no reference curve yet is
# reported as SKIPPED, so the infrastructure can precede the data.
#
# The shieldhit/ and topas/ subdirectories hold mirrored input decks for the
# external codes; they are not registered as tests.

file(GLOB IDD_CASE_DIRS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/idd_*/")

# A real benchmark carries all four documented input files (see the header
# above). Register a directory only when the complete set is present, so a
# half-checked-out idd_*/ directory from a stashed/rebased branch is not turned
# into a phantom failing test (issue #296). A missing reference/ curve is a
# separate, expected state: run_idd.cmake reports such a (complete) case as
# SKIPPED at run time, letting the infrastructure precede the reference data.
set(IDD_REQUIRED_INPUTS geo.dat beam.dat mat.dat detect.dat)

foreach(case_dir ${IDD_CASE_DIRS})
    if(NOT IS_DIRECTORY "${case_dir}")
        continue()
    endif()

    get_filename_component(case_name "${case_dir}" NAME)

    set(_missing_inputs "")
    foreach(_input ${IDD_REQUIRED_INPUTS})
        if(NOT EXISTS "${case_dir}/${_input}")
            list(APPEND _missing_inputs "${_input}")
        endif()
    endforeach()

    # Every idd_*/ directory is meant to be a benchmark, so any missing input is
    # an incomplete case worth surfacing rather than dropping silently.
    if(_missing_inputs)
        string(REPLACE ";" ", " _missing_str "${_missing_inputs}")
        message(WARNING
            "tests/reference: skipping incomplete benchmark '${case_name}' "
            "(missing input file(s): ${_missing_str})")
        continue()
    endif()

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
