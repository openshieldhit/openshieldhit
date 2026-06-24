add_library(osh_transport
    osh_transport.c
    osh_transport_scheduler.c
    osh_neutron_pool.c
    osh_transport_ion.c
    osh_transport_ion_step.c
    osh_transport_neutron.c
    osh_transport_photon.c
)

# Make sure consumers of the library see the headers
target_include_directories(osh_transport
    PUBLIC
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/include
)

# Link math library (Unix/Linux only - Windows has it in the C runtime)
target_link_libraries(osh_transport
    PRIVATE
        osh_common
        osh_physics
        osh_nuclear
        osh_random
        osh_beam
        osh_beam_runtime
        osh_material
        osh_gemca2
        osh_gemca_runtime
        osh_scoring
        $<$<BOOL:${UNIX}>:m>
)
