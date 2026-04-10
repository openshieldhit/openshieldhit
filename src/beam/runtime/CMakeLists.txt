add_library(osh_beam_runtime
    osh_beam_runtime.c
)

target_include_directories(osh_beam_runtime
    PUBLIC
    ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(osh_beam_runtime
    PRIVATE
        osh_beam
        osh_common
)
