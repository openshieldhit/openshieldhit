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
)

# Native SVG plot writer (issue #238 prototype): compiled in only when the
# opt-in OSH_ENABLE_PLOT switch is set, so the stock library carries no
# plotting code.  The dispatcher in osh_scoring_save.c is #ifdef-gated to match.
if(OSH_ENABLE_PLOT)
    target_sources(osh_scoring PRIVATE save/osh_scoring_save_plot.c)
endif()

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
