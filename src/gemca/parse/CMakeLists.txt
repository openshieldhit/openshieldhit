# Gemca2 library
add_library(osh_gemca2_parse
    parse/osh_gemca2_parse.c
    parse/osh_gemca2_parse_body.c
    parse/osh_gemca2_parse_zone.c
    parse/osh_gemca2_parse_medium.c
    parse/osh_gemca2_parse_stack.c
)

# Allow includes like "gemca/..." and "common/..."
target_include_directories(osh_gemca2_parse
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

# common headers/functions, link it
target_link_libraries(osh_gemca2
    PRIVATE
        osh_common
        m
)
