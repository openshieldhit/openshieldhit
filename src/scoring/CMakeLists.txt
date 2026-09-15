add_library(osh_scoring
    osh_scoring.c
    runtime/osh_scoring_accumulator.c
    runtime/osh_scoring_compile.c
    runtime/osh_scoring_estimator.c
    runtime/osh_scoring_estimator_common.c
    runtime/osh_scoring_estimator_point.c
    runtime/osh_scoring_estimator_step.c
    runtime/osh_scoring_filter_runtime.c
    runtime/osh_scoring_geometry_runtime.c
    runtime/osh_scoring_postprocess.c
    runtime/osh_scoring_shadow.c
    runtime/osh_scoring_snapshot.c
    runtime/osh_scoring_step.c
    runtime/osh_scoring_point.c
    save/osh_scoring_save.c
    save/osh_scoring_sink.c
    save/osh_scoring_save_ascii.c
    save/osh_scoring_save_bdo2019_raw.c
    save/osh_scoring_save_bdo2019.c
    save/osh_scoring_save_rtdose.c
    save/osh_scoring_save_plot.c
    save/osh_scoring_plot_spec.c
    save/osh_scoring_svg.c
)

target_include_directories(osh_scoring
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
        ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(osh_scoring
    PUBLIC
        osh_common
        osh_dicom
        osh_raytrace
)

# MCPL phase-space I/O (-DOSH_ENABLE_MCPL, ON by default; see the root
# CMakeLists.txt for why it is switchable).  Everything that includes <mcpl.h>
# lives in this one file, so switching the feature off drops the save layer and
# the vendored Apache-2.0 target together and leaves no third-party code in the
# link.  The osh_scoring_mcpl_record buffer that score_step appends to is plain
# openshieldhit data with no MCPL dependency, so it stays compiled either way;
# osh_scoring_compile() refuses an MCPL detect.dat when OSH_WITH_MCPL is unset,
# which keeps that buffer unreachable at runtime.  PRIVATE: only this library's
# own sources test the macro.
if(OSH_ENABLE_MCPL)
    target_sources(osh_scoring PRIVATE save/osh_scoring_save_mcpl.c)
    target_link_libraries(osh_scoring PRIVATE osh_thirdparty_mcpl)
    target_compile_definitions(osh_scoring PRIVATE OSH_WITH_MCPL)
endif()
