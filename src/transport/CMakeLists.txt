add_library(osh_transport
    osh_transport.c
    prepare/osh_transport_material_prepare.c
)

# Make sure consumers of the library see the headers
target_include_directories(osh_transport
    PUBLIC
    ${PROJECT_SOURCE_DIR}/src
)

# Link math library (Unix/Linux only - Windows has it in the C runtime)
target_link_libraries(osh_transport
    PRIVATE
        osh_physics
        osh_material
        osh_particle
        $<$<BOOL:${UNIX}>:m>
)
