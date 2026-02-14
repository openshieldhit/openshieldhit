# Gemca2 library
add_library(osh_gemca2
    osh_gemca2.c
    osh_gemca2_calc_body.c
    osh_gemca2_calc_surface.c
    osh_gemca2_calc_zone.c
    osh_gemca2_dist.c

    parse/osh_gemca2_parse.c
    parse/osh_gemca2_parse_body.c
    parse/osh_gemca2_parse_zone.c
    parse/osh_gemca2_parse_medium.c
    parse/osh_gemca2_parse_stack.c
)

# Allow includes like "gemca/..." and "common/..."
target_include_directories(osh_gemca2
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

# common headers/functions, link it
target_link_libraries(osh_gemca2
    PRIVATE
        osh_common
)

# Link math library (Unix/Linux only - Windows has it in the C runtime)
if(UNIX)
    target_link_libraries(osh_gemca2 PRIVATE m)
endif()
