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
#     README.md         (optional; human-readable description)
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

    # A real case always has a detect.dat; skip any other directory the glob may
    # catch (e.g. a stray __pycache__/ left by running the Python test helpers),
    # so it is not registered as a broken case.
    if(NOT EXISTS "${case_dir}/detect.dat")
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

# ---- --score-replicas reproducibility contract (issue #230) ------------------
# Real-transport runs of an EM-only and a nuclear case, each transported serially
# and via 1 and N sequential private-accumulator replicas, then cross-compared:
#   serial == --score-replicas 1  (byte-identical) ; == --score-replicas N (tol) ;
#   --score-replicas N run twice   (byte-identical).
# Driven by run_score_replicas.cmake (separate from the plain dry-run case entry).
foreach(repl_case "00_minimal" "13_minimal_nucre")
    add_test(
        NAME "cases::score_replicas::${repl_case}"
        COMMAND ${CMAKE_COMMAND}
            -DOSH_EXECUTABLE=$<TARGET_FILE:openshieldhit>
            -DCASE_DIR=${CMAKE_CURRENT_SOURCE_DIR}/${repl_case}
            -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/score_replicas_${repl_case}
            -DNSTAT=200
            -DREPLICAS=4
            -P ${CMAKE_CURRENT_SOURCE_DIR}/run_score_replicas.cmake
    )
    set_tests_properties("cases::score_replicas::${repl_case}" PROPERTIES LABELS "cases")
endforeach()

# ---- VARIANCE standard-error contract (issue #209) ---------------------------
# Real-transport runs of an EM-only and a nuclear case: a plain run and a
# "VARIANCE N" run are cross-checked (value columns agree; a sane non-negative
# standard-error column is added per quantity), and the VARIANCE run is asserted
# byte-identical across two invocations (the deterministic count-cadence guard of
# #168).  REPLICAS additionally runs a --score-replicas variant to confirm the
# replica split doubles as the batch-means source.  Driven by run_variance.cmake.
foreach(var_case "00_minimal" "13_minimal_nucre")
    add_test(
        NAME "cases::variance::${var_case}"
        COMMAND ${CMAKE_COMMAND}
            -DOSH_EXECUTABLE=$<TARGET_FILE:openshieldhit>
            -DCASE_DIR=${CMAKE_CURRENT_SOURCE_DIR}/${var_case}
            -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/variance_${var_case}
            -DNSTAT=200
            -DNBATCHES=5
            -DREPLICAS=5
            -P ${CMAKE_CURRENT_SOURCE_DIR}/run_variance.cmake
    )
    set_tests_properties("cases::variance::${var_case}" PROPERTIES LABELS "cases")
endforeach()
