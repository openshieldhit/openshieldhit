add_library(osh_app_osh
    osh_beam_parse.c
    osh_app_osh.c
)

target_include_directories(osh_app_osh
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include
    PRIVATE
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(osh_app_osh
    PRIVATE
        osh_beam
        osh_common
        osh_particle
)

if(UNIX)
    target_link_libraries(osh_app_osh PRIVATE m)
endif()
