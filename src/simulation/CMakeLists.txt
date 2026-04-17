add_library(osh_simulation STATIC osh_simulation.c)

target_include_directories(osh_simulation
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include
    PRIVATE
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(osh_simulation
    PUBLIC
        osh_common
    PRIVATE
        osh_beam
        osh_beam_runtime
        osh_gemca2
        osh_gemca_runtime
        osh_material
        osh_scoring
        osh_transport
)

if(UNIX)
    target_link_libraries(osh_simulation PRIVATE m)
endif()
