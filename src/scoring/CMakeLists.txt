add_library(osh_scoring
    osh_scoring.c
    parse/osh_scoring_parse.c
)

target_include_directories(osh_scoring
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(osh_scoring
    PRIVATE
        osh_common
)
