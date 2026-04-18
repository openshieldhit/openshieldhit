add_library(osh_app_osh
    osh_beam_parse.c
    osh_beam_spotlist.c
    osh_material_loaddedx.c
    osh_material_parse.c
    osh_scoring_parse.c
    osh_scoring_parse_filter.c
    osh_scoring_parse_settings.c
    osh_scoring_parse_geometry.c
    osh_scoring_parse_output.c
    osh_geometry_parse.c
    osh_app_osh.c
    osh_run.c
)

target_include_directories(osh_app_osh
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include
    PRIVATE
        ${PROJECT_SOURCE_DIR}/src
        ${PROJECT_BINARY_DIR}/generated
)

target_link_libraries(osh_app_osh
    PUBLIC
        osh_common
    PRIVATE
        osh_beam
        osh_gemca2
        osh_material
        osh_particle
        osh_scoring
        osh_simulation
)

if(UNIX)
    target_link_libraries(osh_app_osh PRIVATE m)
endif()
