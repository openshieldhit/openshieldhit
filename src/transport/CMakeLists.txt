add_library(osh_transport
    osh_transport.c
)

# Make sure consumers of the library see the headers
target_include_directories(osh_transport
    PUBLIC
    ${PROJECT_SOURCE_DIR}/src
)

# Link math library (Unix/Linux only - Windows has it in the C runtime)
target_link_libraries(osh_transport
    PRIVATE
        osh_common
        osh_random
        osh_beam
        osh_material
        osh_gemca2
        osh_scoring
        $<$<BOOL:${UNIX}>:m>
)
